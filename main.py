"""Main execution script for Self-Avoiding Walk (SAW) polymer generation and OpenMM simulation.

Automates:
1. SAW polymer initial structure generation (Homopolymer or Heteropolymer).
2. Energy minimization & NPT equilibrium relaxation via OpenMM.
3. Structure and trajectory export in PDB and DCD formats.
4. Parameter distribution sampling for high-throughput batch generation.
"""

import json
import os
import random
import numpy as np

from src.SAW import SAW
from src.openmm_sim import OpenMMPolymerSimulation

# Load configuration
config_path = "config.json"
with open(config_path, "r") as f:
    data = json.load(f)

aMW = data["average_MW"]  # average molecular weight, unit: g/mol
mass = data["monomer_M"]  # monomer molecular weight, unit: g/mol
chains = data["chains"]  # number of chains
bond_k = data["bond_k"]  # bond spring constant, unit: kcal/mol/angstrom^2
bond_l = data["bond_l"]  # bond length, unit: angstrom
angle_k = data["angle_k"]  # angle spring constant, unit: kcal/mol/radian^2
angle_l = data["angle_l"]  # angle, unit: degree or radian
epsilon = data["epsilon"]  # LJ epsilon, unit: kcal/mol
sigma = data["sigma"]  # LJ sigma, unit: angstrom
mode = data["mode"]  # single_homo, single_heter, sample_homo, sample_heter
device = data.get("device", "cpu")  # cpu, cuda
num_sample = data.get("num_sample", 1)  # number of samples
model_init_save_path = data.get("model_init_save_path", "./model/")

# Simulation parameters
eq_steps = data.get("eq_steps", 20000)
eq_temperature = data.get("eq_temperature", 500.0)
eq_pressure = data.get("eq_pressure", 1.0)
timestep_fs = data.get("timestep_fs", 4.0)
report_interval = data.get("report_interval", 2000)


def create_number(num):
    """Generate a sampled number or range."""
    if isinstance(num, (int, float)):
        return num
    elif isinstance(num, list):
        sample_range = num[1] - num[0]
        sample = random.random()
        return num[0] + sample * sample_range
    else:
        raise TypeError(f"Invalid data type for sampling: {type(num)}")


def run_pipeline():
    os.makedirs(model_init_save_path, exist_ok=True)
    sim = OpenMMPolymerSimulation(device=device)

    print(f"[Pipeline] Initialized OpenMM simulation engine on platform: {sim.platform.getName()}")
    print(f"[Pipeline] Running mode: '{mode}'")

    if mode in ["single_homo", "single_heter"]:
        parameter_package = [bond_k, bond_l, angle_k, angle_l, epsilon, sigma]
        saw_universe = SAW(aMW, mass, chains, bond_l)

        if mode == "single_homo":
            print("[Pipeline] Generating 3D homopolymer model via Self-Avoiding Walk...")
            saw_universe.general_model()
            out_prefix = os.path.join(model_init_save_path, "homo_single")
        else:
            print("[Pipeline] Generating 3D heteropolymer model via Self-Avoiding Walk...")
            saw_universe.general_model_heter()
            out_prefix = os.path.join(model_init_save_path, "heter_single")

        # Build OpenMM System
        topology, system, positions = sim.build_system_from_saw(saw_universe, parameter_package)

        # Run NPT relaxation & export PDB / DCD
        sim.simulate_relaxation(
            topology=topology,
            system=system,
            positions=positions,
            output_prefix=out_prefix,
            temp=eq_temperature,
            pressure=eq_pressure,
            steps=eq_steps,
            timestep_fs=timestep_fs,
            report_interval=report_interval,
        )

    elif mode in ["sample_homo", "sample_heter"]:
        for sample in range(num_sample):
            print(f"\n--- Generating Sample [{sample + 1}/{num_sample}] ---")
            aMW_sample = create_number(aMW)
            mass_sample = create_number(mass)
            chains_sample = int(create_number(chains))
            bond_l_sample = create_number(bond_l)
            bond_k_sample = create_number(bond_k)
            angle_l_sample = create_number(angle_l)
            angle_k_sample = create_number(angle_k)
            epsilon_sample = create_number(epsilon)
            sigma_sample = create_number(sigma)

            parameter_package = [
                bond_k_sample,
                bond_l_sample,
                angle_k_sample,
                angle_l_sample,
                epsilon_sample,
                sigma_sample,
            ]

            saw_universe = SAW(aMW_sample, mass_sample, chains_sample, bond_l_sample)

            if mode == "sample_homo":
                saw_universe.general_model()
                out_prefix = os.path.join(model_init_save_path, f"homo_sample_{sample}")
            else:
                saw_universe.general_model_heter()
                out_prefix = os.path.join(model_init_save_path, f"heter_sample_{sample}")

            topology, system, positions = sim.build_system_from_saw(saw_universe, parameter_package)

            sim.simulate_relaxation(
                topology=topology,
                system=system,
                positions=positions,
                output_prefix=out_prefix,
                temp=eq_temperature,
                pressure=eq_pressure,
                steps=eq_steps,
                timestep_fs=timestep_fs,
                report_interval=report_interval,
            )

    else:
        raise ValueError(f"Unknown mode: {mode}. Must be one of ['single_homo', 'single_heter', 'sample_homo', 'sample_heter']")

    print("\n[Pipeline] All simulations finished successfully!")


if __name__ == "__main__":
    run_pipeline()