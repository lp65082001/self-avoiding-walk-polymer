# self-avoiding-walk-polymer

## Introduction
This project implements a self-avoiding walk for a three-dimensional polymer. The code include two part: create model and equilibrium (using LAMMPS).

## Getting Started
To get started with this project, follow these few steps:

1. **Clone the repository:** Clone this repository to your local machine using the following command:
```
git clone https://github.com/lp65082001/self-avoiding-walk-polymer.git
```
2. **Build environment** Build the environment with docker:
```
docker build -f ./Dockerfile -t saw:v1 .
```
Run container (gpu version):
```
sudo docker run -i -t  --gpus all -v /in this folder/:/root/ --name=saw_cuda saw:v1
```
Add environment parameter to pass mpi under root:
```
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
```

3. **Install LAMMPS** Compile LAMMPS with gpu (option)
```
wget https://lammps.sandia.gov/tars/lammps-stable.tar.gz
```
unzip LAMMPS:
```
tar-xzvf lammps*.tar.gz
```
compile with gpu:
```
cd lammps                # change to the LAMMPS distribution directory
mkdir build; cd build    # create and use a build directory
cmake -D BUILD_MPI=yes -D LAMMPS_MACHINE=mpi -D PKG_GPU=on -D GPU_API=cuda -D GPU_ARCH=sm_75 -D PKG_OPT=yes -D PKG_PERI=yes -D PKG_KSPACE=yes -D PKG_MC=yes -D PKG_MOLECULE=yes -D PKG_RIGID=yes -D PKG_MPIIO=yes -D PKG_VORONOI=yes -D DOWNLOAD_PLUMED=yes -D PKG_USER-PLUMED=yes -D PKG_USER-OMP=yes ../cmake
```
cmake install:
```
cmake --build .          # compilation (or type "make")
```

## Features
- Create polymer model by using self avoiding walk
- Automate to using LAMMPS
- Automate create sample with random parameter
