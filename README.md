# self-avoiding-walk-polymer

## Introduction
This project implements a self-avoiding walk for a three-dimensional polymer. The code include two part: create model and equilibrium (using LAMMPS).

## Getting Started
To get started with this project, follow these few steps:

1. **Clone the repository:** Clone this repository to your local machine using the following command:
```
git clone https://github.com/lp65082001/self-avoiding-walk-polymer.git
```
2. **Build environment**
First, build the environment with docker:
```
docker build -f ./Dockerfile -t saw:v1 .
```
Second, run container (gpu version):
```
sudo docker run -i -t  --gpus all -v /in this folder/:/root/ --name=saw_cuda saw:v1
```
Third, add environment parameter to pass mpi under root:
```
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
```
Four, download LAMMPS:
```
wget https://lammps.sandia.gov/tars/lammps-stable.tar.gz
```
Fifth, compile LAMMPS with gpu:

## Features
- Create polymer model by using self avoiding walk
- Collaboration with other developers
