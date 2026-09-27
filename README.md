# self-avoiding-walk-polymer

A 3D polymer initial structure generator and molecular dynamics pipeline using **Self-Avoiding Walk (SAW)** and **OpenMM**.

---

## Overview

1. **Initial Model Generation**:
   - Constructs coarse-grained or atomistic polymer chains in 3D lattice space via Self-Avoiding Walk (SAW).
   - Eliminates steric clashes and initial atom overlap.
   - Supports both monodisperse homopolymers and polydisperse heteropolymers (Gaussian chain length distributions).
2. **Equilibration (OpenMM)**:
   - Energy minimization (`LocalEnergyMinimizer`).
   - NPT ensemble equilibrium relaxation using `LangevinMiddleIntegrator` and `MonteCarloBarostat`.
3. **Mechanical Testing (OpenMM)**:
   - Automated uniaxial tensile deformation using `MonteCarloAnisotropicBarostat` (constant tensile stress / ramping).
   - Real-time strain tracking and Young's modulus calculation.
4. **Standard File Formats**:
   - Outputs structures as standard **PDB** (`.pdb`) files.
   - Outputs trajectories as standard **DCD** (`.dcd`) files.

---

## Getting Started

### 1. Installation

#### Native Python
Install dependencies via pip:
```bash
pip install -r requirements.txt
```

#### Docker
Build the container image:
```bash
docker build -t saw_openmm:v1 .
```

Run container with GPU support:
```bash
docker run --gpus all -it -v $(pwd):/workspace saw_openmm:v1
```

Or run container on CPU:
```bash
docker run -it -v $(pwd):/workspace saw_openmm:v1
```

---

## Configuration (`config.json`)

All physical parameters, polymer dimensions, and simulation settings are configured in `config.json`:

```json
{
    "average_MW": 4450,
    "monomer_M": 44.5,
    "chains": 10,
    "bond_k": 1463.4,
    "bond_l": 2.6,
    "angle_k": 2.0,
    "angle_l": 2.0,
    "epsilon": 0.4988,
    "sigma": 4.628,
    "num_sample": 1,
    "mode": "single_homo",
    "device": "cpu",
    "model_init_save_path": "./model/",
    "eq_steps": 10000,
    "eq_temperature": 500.0,
    "eq_pressure": 1.0,
    "timestep_fs": 4.0,
    "report_interval": 1000
}
```

### Parameter Description
| Parameter | Unit | Description |
| :--- | :--- | :--- |
| `average_MW` | g/mol | Average molecular weight of polymer chains |
| `monomer_M` | g/mol | Molecular weight of a single monomer unit |
| `chains` | - | Number of polymer chains in the simulation box |
| `bond_k` | kcal/mol/Å² | Harmonic bond force constant |
| `bond_l` | Å | Equilibrium bond length |
| `angle_k` | kcal/mol/rad² | Harmonic angle force constant |
| `angle_l` | degree or rad | Equilibrium bond angle |
| `epsilon` | kcal/mol | Lennard-Jones potential depth |
| `sigma` | Å | Lennard-Jones finite distance at zero potential |
| `mode` | string | `single_homo`, `single_heter`, `sample_homo`, or `sample_heter` |
| `device` | string | `cpu` or `cuda` (GPU accelerated) |
| `eq_steps` | steps | Number of equilibration MD simulation steps |
| `eq_temperature` | K | Equilibration temperature |
| `eq_pressure` | bar | Equilibration pressure |
| `timestep_fs` | fs | Integration timestep |

---

## Usage

### 1. Generate Polymer & Run Equilibration

Execute `main.py`:
```bash
python main.py
```

Generated outputs will be saved in `./model/`:
- `<name>_init.pdb`: Initial polymer configuration from Self-Avoiding Walk.
- `<name>_eq.pdb`: Equilibrated polymer structure after NPT relaxation.
- `<name>_eq.dcd`: Full trajectory of the equilibration simulation.
- `<name>_sys.xml`: OpenMM System serialized definition.

### 2. Run Tensile Test & Calculate Young's Modulus

Execute `tensile_test.py`:
```bash
python analysis_script_tensile_test/tensile_test.py
```

The script will automatically detect models in `analysis_script_tensile_test/datafile/`, apply uniaxial tensile stress (30 MPa in X direction at 300 K), and record deformation:
- `output/<name>_tensile.pdb`: Deformed structure under tensile strain.
- `output/<name>_tensile.dcd`: Tensile deformation trajectory.
- `properties.json`: Computed strain and Young's modulus (MPa).

---

## Features
- **Self-Avoiding Walk (SAW)**: Eliminates initial steric overlap.
- **Native OpenMM Integration**: Pure Python MD pipeline with zero external compilation requirements.
- **Hardware Acceleration**: Automatic GPU (CUDA) and multi-core CPU platform support.
- **Standard Molecular Formats**: Direct PDB structure and DCD trajectory support.
- **High-Throughput Sampling**: Batch generation with randomized parameter distributions.
