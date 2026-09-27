"""OpenMM Simulation Engine for Self-Avoiding Walk (SAW) Polymers.

Replaces LAMMPS with native OpenMM Python APIs for:
1. Topology and System construction from SAW configurations or PDB/data files.
2. Energy minimization (LocalEnergyMinimizer).
3. Equilibrium relaxation under NPT ensemble (LangevinMiddleIntegrator + MonteCarloBarostat).
4. Uniaxial tensile deformation & Young's Modulus calculation (MonteCarloAnisotropicBarostat).
5. Trajectory and structure export in PDB and DCD formats.
"""

import os
import sys
import json
import numpy as np
import openmm
from openmm import app, unit


def make_periodic_box_vectors(lx, ly, lz):
    """Return a single Quantity of 3 Vec3s in nanometers, as expected by Topology.setPeriodicBoxVectors."""
    def to_nm(v):
        if hasattr(v, "value_in_unit"):
            return v.value_in_unit(unit.nanometer)
        return float(v) * 0.1  # assume Angstroms if dimensionless float

    return (
        openmm.Vec3(to_nm(lx), 0.0, 0.0),
        openmm.Vec3(0.0, to_nm(ly), 0.0),
        openmm.Vec3(0.0, 0.0, to_nm(lz)),
    ) * unit.nanometer


def set_topology_box(topology, box_vectors):
    """Safely apply periodic box vectors to an OpenMM Topology."""
    if box_vectors is None:
        return
    if isinstance(box_vectors, (list, tuple)) and len(box_vectors) == 3:
        def to_nm(v):
            if hasattr(v, "value_in_unit"):
                return v.value_in_unit(unit.nanometer)
            elif isinstance(v, (list, tuple, openmm.Vec3)):
                return [to_nm(x) for x in v]
            return float(v) * 0.1

        v0 = to_nm(box_vectors[0])
        v1 = to_nm(box_vectors[1])
        v2 = to_nm(box_vectors[2])
        topology.setPeriodicBoxVectors(
            (openmm.Vec3(*v0), openmm.Vec3(*v1), openmm.Vec3(*v2)) * unit.nanometer
        )


class OpenMMPolymerSimulation:
    def __init__(self, device="cpu"):
        """Initialize OpenMM simulation platform.

        Args:
            device (str): "cuda", "cpu", or "reference".
        """
        self.device = device.lower()
        self.platform = self._select_platform(self.device)

    def _select_platform(self, device):
        available = [
            openmm.Platform.getPlatform(i).getName()
            for i in range(openmm.Platform.getNumPlatforms())
        ]
        target = "CUDA" if device == "cuda" else "CPU"
        if target in available:
            return openmm.Platform.getPlatformByName(target)
        elif "CPU" in available:
            print(f"[Warning] Platform '{target}' not available, falling back to 'CPU'. Available: {available}")
            return openmm.Platform.getPlatformByName("CPU")
        else:
            print(f"[Warning] Falling back to '{available[0]}'.")
            return openmm.Platform.getPlatformByName(available[0])

    def build_system_from_saw(self, saw, parameter_package):
        """Construct an OpenMM Topology, System, and Positions from a SAW instance.

        Args:
            saw (SAW): Generated SAW instance with particule_table and dimensions.
            parameter_package (list): [bond_k, bond_l, angle_k, angle_l, epsilon, sigma]
                - bond_k: kcal/mol/A^2 (LAMMPS harmonic bond constant)
                - bond_l: Angstrom
                - angle_k: kcal/mol/rad^2
                - angle_l: degrees or radians
                - epsilon: kcal/mol
                - sigma: Angstrom

        Returns:
            tuple: (topology, system, positions)
        """
        bond_k, bond_l, angle_k, angle_l, epsilon, sigma = parameter_package

        # System setup
        system = openmm.System()
        box_vectors = make_periodic_box_vectors(saw.init_box_size, saw.init_box_size, saw.init_box_size)
        system.setDefaultPeriodicBoxVectors(*box_vectors)

        topology = app.Topology()
        topology.setPeriodicBoxVectors(box_vectors)

        # Determine chain monomer lengths
        if hasattr(saw, "number_unit_list") and saw.number_unit_list is not None:
            chain_lengths = saw.number_unit_list.tolist()
        else:
            chain_lengths = [saw.res] * saw.chains

        total_atoms = saw.particule_table.shape[0]
        for _ in range(total_atoms):
            system.addParticle(saw.mass * unit.dalton)

        # Create topology atoms and chains
        atom_list = []
        atom_idx = 0
        for chain_idx, length in enumerate(chain_lengths):
            chain = topology.addChain(id=str(chain_idx + 1))
            for _ in range(length):
                res = topology.addResidue("POL", chain)
                atom = topology.addAtom("C", app.Element.getBySymbol("C"), res)
                atom_list.append(atom)
                atom_idx += 1

        # Forces
        # OpenMM HarmonicBondForce uses E = 0.5 * k * (r - r0)^2
        # LAMMPS uses E = K * (r - r0)^2 => k_openmm = 2 * K_lammps
        bonds_force = openmm.HarmonicBondForce()
        bonds_force.setUsesPeriodicBoundaryConditions(True)
        k_bond_openmm = 2.0 * bond_k * unit.kilocalories_per_mole / (unit.angstrom**2)
        r0_bond_openmm = bond_l * unit.angstrom

        # HarmonicAngleForce uses E = 0.5 * k * (theta - theta0)^2
        angles_force = openmm.HarmonicAngleForce()
        angles_force.setUsesPeriodicBoundaryConditions(True)
        k_angle_openmm = 2.0 * angle_k * unit.kilocalories_per_mole / (unit.radian**2)
        if angle_l > np.pi:
            theta0_openmm = angle_l * unit.degree
        else:
            theta0_openmm = angle_l * unit.radian

        # NonbondedForce (Lennard-Jones)
        nb_force = openmm.NonbondedForce()
        nb_force.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
        nb_force.setCutoffDistance(20.0 * unit.angstrom)

        for _ in range(total_atoms):
            nb_force.addParticle(
                0.0 * unit.elementary_charge,
                sigma * unit.angstrom,
                epsilon * unit.kilocalories_per_mole,
            )

        # Add bonds and angles based on chain structure
        bond_pairs = []
        current_offset = 0
        for length in chain_lengths:
            for i in range(length - 1):
                idx1 = current_offset + i
                idx2 = idx1 + 1
                bonds_force.addBond(idx1, idx2, r0_bond_openmm, k_bond_openmm)
                topology.addBond(atom_list[idx1], atom_list[idx2])
                bond_pairs.append((idx1, idx2))

            for i in range(length - 2):
                idx1 = current_offset + i
                idx2 = idx1 + 1
                idx3 = idx1 + 2
                angles_force.addAngle(idx1, idx2, idx3, theta0_openmm, k_angle_openmm)

            current_offset += length

        # Exclude 1-2 and 1-3 pairs from nonbonded interactions
        nb_force.createExceptionsFromBonds(bond_pairs, 0.0, 0.0)

        system.addForce(bonds_force)
        system.addForce(angles_force)
        system.addForce(nb_force)

        # Positions
        positions = [
            openmm.Vec3(float(row[0]), float(row[1]), float(row[2]))
            for row in saw.particule_table
        ] * unit.angstrom

        return topology, system, positions

    def build_system_from_lammps_data(self, data_path):
        """Parse a LAMMPS .data file and reconstruct the OpenMM System, Topology, and Positions.

        Supports standard angle-style polymer data generated by SAW / LAMMPS.
        """
        with open(data_path, "r") as f:
            lines = [line.strip() for line in f.readlines()]

        num_atoms = 0
        num_bonds = 0
        num_angles = 0
        xlo, xhi = 0.0, 100.0
        ylo, yhi = 0.0, 100.0
        zlo, zhi = 0.0, 100.0
        mass = 44.5
        epsilon = 0.4988
        sigma = 4.628
        bond_k = 1463.4
        bond_l = 2.6
        angle_k = 2.0
        angle_theta = 180.0

        i = 0
        while i < len(lines):
            line = lines[i]
            if "atoms" in line and not line.startswith("Atoms"):
                num_atoms = int(line.split()[0])
            elif "bonds" in line and not line.startswith("Bonds"):
                num_bonds = int(line.split()[0])
            elif "angles" in line and not line.startswith("Angles"):
                num_angles = int(line.split()[0])
            elif "xlo xhi" in line:
                parts = line.split()
                xlo, xhi = float(parts[0]), float(parts[1])
            elif "ylo yhi" in line:
                parts = line.split()
                ylo, yhi = float(parts[0]), float(parts[1])
            elif "zlo zhi" in line:
                parts = line.split()
                zlo, zhi = float(parts[0]), float(parts[1])
            elif line.startswith("Masses"):
                i += 1
                while i < len(lines) and not lines[i]:
                    i += 1
                if i < len(lines) and lines[i]:
                    mass = float(lines[i].split()[1])
            elif line.startswith("Pair Coeffs"):
                i += 1
                while i < len(lines) and not lines[i]:
                    i += 1
                if i < len(lines) and lines[i]:
                    parts = lines[i].split()
                    epsilon = float(parts[1])
                    sigma = float(parts[2])
            elif line.startswith("Bond Coeffs"):
                i += 1
                while i < len(lines) and not lines[i]:
                    i += 1
                if i < len(lines) and lines[i]:
                    parts = lines[i].split()
                    bond_k = float(parts[1])
                    bond_l = float(parts[2])
            elif line.startswith("Angle Coeffs"):
                i += 1
                while i < len(lines) and not lines[i]:
                    i += 1
                if i < len(lines) and lines[i]:
                    parts = lines[i].split()
                    angle_k = float(parts[1])
                    angle_theta = float(parts[2])
            i += 1

        lx = (xhi - xlo) * unit.angstrom
        ly = (yhi - ylo) * unit.angstrom
        lz = (zhi - zlo) * unit.angstrom

        system = openmm.System()
        box_vectors = make_periodic_box_vectors(lx, ly, lz)
        system.setDefaultPeriodicBoxVectors(*box_vectors)

        topology = app.Topology()
        topology.setPeriodicBoxVectors(box_vectors)

        for _ in range(num_atoms):
            system.addParticle(mass * unit.dalton)

        # Parse Atoms, Bonds, Angles sections
        atoms_data = {}
        bonds_data = []
        angles_data = []

        mode = None
        for line in lines:
            if not line:
                continue
            if line.startswith("Atoms"):
                mode = "atoms"
                continue
            elif line.startswith("Bonds"):
                mode = "bonds"
                continue
            elif line.startswith("Angles"):
                mode = "angles"
                continue
            elif line.startswith("Velocities"):
                mode = "velocities"
                continue
            elif line.startswith("Masses") or line.startswith("Pair Coeffs") or line.startswith("Bond Coeffs") or line.startswith("Angle Coeffs"):
                mode = None
                continue

            parts = line.split()
            if mode == "atoms":
                if len(parts) >= 6:
                    atom_id = int(parts[0])
                    mol_id = int(parts[1])
                    x = float(parts[3]) - xlo
                    y = float(parts[4]) - ylo
                    z = float(parts[5]) - zlo
                    atoms_data[atom_id] = (mol_id, x, y, z)
            elif mode == "bonds":
                if len(parts) >= 4:
                    a1 = int(parts[2]) - 1
                    a2 = int(parts[3]) - 1
                    bonds_data.append((a1, a2))
            elif mode == "angles":
                if len(parts) >= 5:
                    a1 = int(parts[2]) - 1
                    a2 = int(parts[3]) - 1
                    a3 = int(parts[4]) - 1
                    angles_data.append((a1, a2, a3))

        # Build positions and topology
        positions = [None] * num_atoms
        chain_map = {}
        atom_objects = [None] * num_atoms

        sorted_ids = sorted(atoms_data.keys())
        for atom_id in sorted_ids:
            idx = atom_id - 1
            mol_id, x, y, z = atoms_data[atom_id]
            positions[idx] = openmm.Vec3(float(x), float(y), float(z))
            if mol_id not in chain_map:
                chain = topology.addChain(id=str(mol_id))
                chain_map[mol_id] = chain
            chain = chain_map[mol_id]
            res = topology.addResidue("POL", chain)
            atom_obj = topology.addAtom("C", app.Element.getBySymbol("C"), res)
            atom_objects[idx] = atom_obj

        bonds_force = openmm.HarmonicBondForce()
        bonds_force.setUsesPeriodicBoundaryConditions(True)
        k_bond_openmm = 2.0 * bond_k * unit.kilocalories_per_mole / (unit.angstrom**2)
        r0_bond_openmm = bond_l * unit.angstrom
        for a1, a2 in bonds_data:
            bonds_force.addBond(a1, a2, r0_bond_openmm, k_bond_openmm)
            topology.addBond(atom_objects[a1], atom_objects[a2])

        angles_force = openmm.HarmonicAngleForce()
        angles_force.setUsesPeriodicBoundaryConditions(True)
        k_angle_openmm = 2.0 * angle_k * unit.kilocalories_per_mole / (unit.radian**2)
        theta0_openmm = angle_theta * unit.degree if angle_theta > np.pi else angle_theta * unit.radian
        for a1, a2, a3 in angles_data:
            angles_force.addAngle(a1, a2, a3, theta0_openmm, k_angle_openmm)

        nb_force = openmm.NonbondedForce()
        nb_force.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
        nb_force.setCutoffDistance(20.0 * unit.angstrom)
        for _ in range(num_atoms):
            nb_force.addParticle(
                0.0 * unit.elementary_charge,
                sigma * unit.angstrom,
                epsilon * unit.kilocalories_per_mole,
            )
        nb_force.createExceptionsFromBonds(bonds_data, 0.0, 0.0)

        system.addForce(bonds_force)
        system.addForce(angles_force)
        system.addForce(nb_force)

        return topology, system, positions * unit.angstrom

    def simulate_relaxation(
        self,
        topology,
        system,
        positions,
        output_prefix,
        temp=500.0,
        pressure=1.0,
        steps=50000,
        timestep_fs=4.0,
        report_interval=5000,
    ):
        """Perform energy minimization and NPT equilibrium relaxation.

        Exports:
            - Initial configuration: {output_prefix}_init.pdb
            - Equilibration trajectory: {output_prefix}_eq.dcd
            - Equilibrated final configuration: {output_prefix}_eq.pdb
            - OpenMM System XML: {output_prefix}_sys.xml
        """
        os.makedirs(os.path.dirname(os.path.abspath(output_prefix)), exist_ok=True)

        # Write initial unminimized PDB
        init_pdb = f"{output_prefix}_init.pdb"
        with open(init_pdb, "w") as f:
            app.PDBFile.writeFile(topology, positions, f)
        print(f"[OpenMM] Saved initial structure to: {init_pdb}")

        # Add isotropic MonteCarloBarostat for NPT
        barostat = openmm.MonteCarloBarostat(
            pressure * unit.bar, temp * unit.kelvin, 25
        )
        system.addForce(barostat)

        # Integrator
        integrator = openmm.LangevinMiddleIntegrator(
            temp * unit.kelvin,
            1.0 / unit.picosecond,
            timestep_fs * unit.femtosecond,
        )

        simulation = app.Simulation(topology, system, integrator, self.platform)
        simulation.context.setPositions(positions)

        # Energy minimization
        print("[OpenMM] Minimizing energy...")
        simulation.minimizeEnergy(maxIterations=2000)

        # Velocity initialization
        simulation.context.setVelocitiesToTemperature(temp * unit.kelvin)

        # Reporters
        dcd_path = f"{output_prefix}_eq.dcd"
        simulation.reporters.append(app.DCDReporter(dcd_path, report_interval))
        simulation.reporters.append(
            app.StateDataReporter(
                sys.stdout,
                report_interval,
                step=True,
                potentialEnergy=True,
                temperature=True,
                volume=True,
                speed=True,
            )
        )

        print(f"[OpenMM] Running NPT equilibration for {steps} steps (T={temp}K, P={pressure}bar)...")
        simulation.step(steps)

        # Extract final state
        final_state = simulation.context.getState(
            getPositions=True, enforcePeriodicBox=True
        )
        final_positions = final_state.getPositions()
        final_box = final_state.getPeriodicBoxVectors()
        set_topology_box(topology, final_box)

        eq_pdb = f"{output_prefix}_eq.pdb"
        with open(eq_pdb, "w") as f:
            app.PDBFile.writeFile(topology, final_positions, f)
        print(f"[OpenMM] Equilibrated structure saved to: {eq_pdb}")
        print(f"[OpenMM] Trajectory saved to: {dcd_path}")

        # Save OpenMM System XML for subsequent simulations
        xml_path = f"{output_prefix}_sys.xml"
        with open(xml_path, "w") as f:
            f.write(openmm.XmlSerializer.serialize(system))
        print(f"[OpenMM] System definition saved to: {xml_path}")

        return eq_pdb, dcd_path, final_positions

    def simulate_tensile_test(
        self,
        topology,
        system,
        positions,
        output_prefix,
        temp=300.0,
        stress_mpa=30.0,
        steps=50000,
        timestep_fs=4.0,
        report_interval=2000,
    ):
        """Run uniaxial tensile test along the X direction.

        Applies tensile stress (negative pressure) along X while keeping Y and Z at atmospheric pressure:
        Px = -stress_mpa * 10 bar (e.g. -300 bar for 30 MPa tensile stress)
        Py = 1.0 bar, Pz = 1.0 bar

        Exports:
            - Tensile trajectory: {output_prefix}_tensile.dcd
            - Deformed structure: {output_prefix}_tensile.pdb
        """
        os.makedirs(os.path.dirname(os.path.abspath(output_prefix)), exist_ok=True)

        # Remove existing isotropic barostat if present
        for i in range(system.getNumForces() - 1, -1, -1):
            force = system.getForce(i)
            if isinstance(force, openmm.MonteCarloBarostat):
                system.removeForce(i)

        # Add Anisotropic Barostat for tensile test
        # 1 MPa = 10 bar => 30 MPa = 300 bar tensile stress -> Px = -300 bar
        px_bar = -float(stress_mpa) * 10.0
        barostat = openmm.MonteCarloAnisotropicBarostat(
            openmm.Vec3(px_bar, 1.0, 1.0) * unit.bar,
            temp * unit.kelvin,
            True,  # scaleX
            True,  # scaleY
            True,  # scaleZ
            25,    # frequency
            False  # scaleMoleculesAsRigid
        )
        system.addForce(barostat)

        integrator = openmm.LangevinMiddleIntegrator(
            temp * unit.kelvin,
            1.0 / unit.picosecond,
            timestep_fs * unit.femtosecond,
        )

        simulation = app.Simulation(topology, system, integrator, self.platform)
        simulation.context.setPositions(positions)

        # Minimize energy slightly at test temperature
        simulation.minimizeEnergy(maxIterations=500)
        simulation.context.setVelocitiesToTemperature(temp * unit.kelvin)

        dcd_path = f"{output_prefix}_tensile.dcd"
        simulation.reporters.append(app.DCDReporter(dcd_path, report_interval))

        # Track box lengths
        initial_state = simulation.context.getState(getPositions=True)
        initial_box = initial_state.getPeriodicBoxVectors()
        lx0 = initial_box[0][0].value_in_unit(unit.angstrom)

        print(f"[OpenMM Tensile] Starting test at T={temp}K, Tensile Stress={stress_mpa} MPa (Px={px_bar} bar)...")
        print(f"[OpenMM Tensile] Initial Lx = {lx0:.3f} A")

        steps_run = 0
        step_chunk = report_interval
        box_history = []

        while steps_run < steps:
            simulation.step(min(step_chunk, steps - steps_run))
            steps_run += min(step_chunk, steps - steps_run)
            state = simulation.context.getState(getPositions=True)
            box = state.getPeriodicBoxVectors()
            lx = box[0][0].value_in_unit(unit.angstrom)
            ly = box[1][1].value_in_unit(unit.angstrom)
            lz = box[2][2].value_in_unit(unit.angstrom)
            cur_strain = (lx - lx0) / lx0
            box_history.append((steps_run, lx, ly, lz, cur_strain))
            print(f"Step {steps_run:6d}/{steps} | Lx: {lx:.2f} A | Strain: {cur_strain:+.4f}")

        final_state = simulation.context.getState(getPositions=True, enforcePeriodicBox=True)
        final_box = final_state.getPeriodicBoxVectors()
        lxf = final_box[0][0].value_in_unit(unit.angstrom)

        strain = (lxf - lx0) / lx0
        if abs(strain) > 1e-6:
            young_modulus_mpa = stress_mpa / strain
        else:
            young_modulus_mpa = 0.0

        deformed_pdb = f"{output_prefix}_tensile.pdb"
        set_topology_box(topology, final_box)
        with open(deformed_pdb, "w") as f:
            app.PDBFile.writeFile(topology, final_state.getPositions(), f)

        print(f"[OpenMM Tensile] Finished. Final Lx = {lxf:.3f} A | Strain = {strain:.4f} | Young's Modulus = {young_modulus_mpa:.2f} MPa")
        print(f"[OpenMM Tensile] Trajectory: {dcd_path}, Structure: {deformed_pdb}")

        return {
            "lx_initial": lx0,
            "lx_final": lxf,
            "strain": strain,
            "young_modulus_mpa": young_modulus_mpa,
            "tensile_pdb": deformed_pdb,
            "tensile_dcd": dcd_path,
        }
