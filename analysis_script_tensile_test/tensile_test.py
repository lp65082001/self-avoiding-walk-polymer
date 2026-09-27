"""OpenMM Tensile Test & Young's Modulus Analysis Script.

Replaces the legacy LAMMPS tensile test with native OpenMM:
1. Loads polymer structures from PDB (or legacy LAMMPS data files).
2. Applies uniaxial tensile stress (-300 bar = 30 MPa in X) at 300 K using MonteCarloAnisotropicBarostat.
3. Records stress-strain behavior and saves:
   - Output trajectory in DCD format (*_tensile.dcd)
   - Final deformed structure in PDB format (*_tensile.pdb)
   - Calculated Young's Modulus and strain in properties.json.
"""

import os
import sys
import glob
import json
import numpy as np

# Ensure parent directory is in path to import src
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
if parent_dir not in sys.path:
    sys.path.insert(0, parent_dir)

from src.openmm_sim import OpenMMPolymerSimulation, set_topology_box
import openmm
from openmm import app, unit, XmlSerializer


def tensile_test(
    model_path,
    output_dir="./output_tensile",
    stress_mpa=30.0,
    temp=300.0,
    steps=10000,
    timestep_fs=4.0,
    report_interval=1000,
    device="cpu",
):
    """Run OpenMM uniaxial tensile test on a given model file (PDB or LAMMPS .data).

    Args:
        model_path (str): Path to input PDB or .data file.
        output_dir (str): Directory where DCD and deformed PDB will be stored.
        stress_mpa (float): Tensile stress applied along the X direction in MPa (default 30 MPa = 300 bar).
        temp (float): Simulation temperature in Kelvin (default 300 K).
        steps (int): Number of MD simulation steps (default 10000).
        timestep_fs (float): Integration timestep in femtoseconds (default 4.0 fs).
        report_interval (int): Reporting step interval (default 1000).
        device (str): OpenMM platform device ('cpu' or 'cuda').

    Returns:
        dict: Test results containing strain, young_modulus_mpa, deformed PDB path, DCD path.
    """
    os.makedirs(output_dir, exist_ok=True)
    sim = OpenMMPolymerSimulation(device=device)
    base_name = os.path.splitext(os.path.basename(model_path))[0]
    out_prefix = os.path.join(output_dir, base_name)

    print(f"\n==================================================")
    print(f"[Tensile Test] Starting test on: {model_path}")
    print(f"==================================================")

    # 1. Load System, Topology, and Positions
    if model_path.endswith(".pdb"):
        pdb = app.PDBFile(model_path)
        topology = pdb.topology
        positions = pdb.positions

        # Look for companion XML system definition
        xml_path = os.path.splitext(model_path)[0] + "_sys.xml"
        if not os.path.exists(xml_path):
            xml_path = os.path.join(os.path.dirname(model_path), base_name.replace("_eq", "") + "_sys.xml")

        if os.path.exists(xml_path):
            print(f"[Tensile Test] Loading system definition from: {xml_path}")
            with open(xml_path, "r") as f:
                system = XmlSerializer.deserialize(f.read())
        else:
            # Reconstruct system using parameters from config.json or defaults
            print(f"[Tensile Test] Companion XML not found. Building system from parameters...")
            config_path = os.path.join(parent_dir, "config.json")
            if os.path.exists(config_path):
                with open(config_path, "r") as f:
                    cfg = json.load(f)
            else:
                cfg = {}

            mass = cfg.get("monomer_M", 44.5)
            bond_k = cfg.get("bond_k", 1463.4)
            bond_l = cfg.get("bond_l", 2.6)
            angle_k = cfg.get("angle_k", 2.0)
            angle_l = cfg.get("angle_l", 2.0)
            epsilon = cfg.get("epsilon", 0.4988)
            sigma = cfg.get("sigma", 4.628)

            system = openmm.System()
            box_vecs = topology.getPeriodicBoxVectors()
            if box_vecs is not None:
                system.setDefaultPeriodicBoxVectors(*box_vecs)

            num_atoms = topology.getNumAtoms()
            for _ in range(num_atoms):
                system.addParticle(mass * unit.dalton)

            bonds_force = openmm.HarmonicBondForce()
            bonds_force.setUsesPeriodicBoundaryConditions(True)
            k_bond = 2.0 * bond_k * unit.kilocalories_per_mole / (unit.angstrom**2)
            r0_bond = bond_l * unit.angstrom

            angles_force = openmm.HarmonicAngleForce()
            angles_force.setUsesPeriodicBoundaryConditions(True)
            k_angle = 2.0 * angle_k * unit.kilocalories_per_mole / (unit.radian**2)
            theta0 = angle_l * unit.degree if angle_l > np.pi else angle_l * unit.radian

            bond_pairs = []
            for b in topology.bonds():
                i1, i2 = b.atom1.index, b.atom2.index
                bonds_force.addBond(i1, i2, r0_bond, k_bond)
                bond_pairs.append((i1, i2))

            # Detect consecutive angles along chains
            for chain in topology.chains():
                atoms_in_chain = list(chain.atoms())
                for i in range(len(atoms_in_chain) - 2):
                    angles_force.addAngle(
                        atoms_in_chain[i].index,
                        atoms_in_chain[i + 1].index,
                        atoms_in_chain[i + 2].index,
                        theta0,
                        k_angle,
                    )

            nb_force = openmm.NonbondedForce()
            nb_force.setNonbondedMethod(openmm.NonbondedForce.CutoffPeriodic)
            nb_force.setCutoffDistance(20.0 * unit.angstrom)
            for _ in range(num_atoms):
                nb_force.addParticle(0.0 * unit.elementary_charge, sigma * unit.angstrom, epsilon * unit.kilocalories_per_mole)
            nb_force.createExceptionsFromBonds(bond_pairs, 0.0, 0.0)

            system.addForce(bonds_force)
            system.addForce(angles_force)
            system.addForce(nb_force)

    elif model_path.endswith(".data"):
        print(f"[Tensile Test] Parsing legacy LAMMPS data file...")
        topology, system, positions = sim.build_system_from_lammps_data(model_path)
    else:
        raise ValueError(f"Unsupported model file format: {model_path}. Must be .pdb or .data")

    # 2. Run Tensile Test Simulation
    results = sim.simulate_tensile_test(
        topology=topology,
        system=system,
        positions=positions,
        output_prefix=out_prefix,
        temp=temp,
        stress_mpa=stress_mpa,
        steps=steps,
        timestep_fs=timestep_fs,
        report_interval=report_interval,
    )

    return results


if __name__ == "__main__":
    mypath = os.path.join(current_dir, "datafile")
    out_dir = os.path.join(current_dir, "output")
    os.makedirs(out_dir, exist_ok=True)

    # Search for PDB files first, fallback to .data if no PDBs
    pdb_files = [os.path.join(mypath, f) for f in os.listdir(mypath) if f.endswith(".pdb")]
    if not pdb_files:
        pdb_files = [os.path.join(mypath, f) for f in os.listdir(mypath) if f.endswith(".data")]

    print(f"Found {len(pdb_files)} model(s) to test in {mypath}: {pdb_files}")

    data_all = {}
    for fullpath in pdb_files:
        try:
            res = tensile_test(
                model_path=fullpath,
                output_dir=out_dir,
                stress_mpa=30.0,
                temp=300.0,
                steps=5000,
                timestep_fs=4.0,
                report_interval=1000,
                device="cpu",
            )
            data_all[fullpath] = {
                "strain": res["strain"],
                "young_modulus_mpa": res["young_modulus_mpa"],
                "lx_initial_angstrom": res["lx_initial"],
                "lx_final_angstrom": res["lx_final"],
                "deformed_pdb": res["tensile_pdb"],
                "tensile_dcd": res["tensile_dcd"],
            }
        except Exception as e:
            print(f"[Error] Failed to test {fullpath}: {e}")
            import traceback
            traceback.print_exc()

    prop_path = os.path.join(current_dir, "properties.json")
    with open(prop_path, "w") as f:
        json.dump(data_all, f, indent=4)
    print(f"\n[Done] All tensile test properties successfully saved to: {prop_path}")
