from src.SAW import SAW
import random
import numpy as np
import json
import os

# get information from config #
with open('config.json') as f:
    data = json.load(f)

aMW = data["average_MW"] # unit: g/mol
mass = data["monomer_M"] # unit: g/mol
chains = data["chains"]
bond_k = data["bond_k"] # unit: kcal/mol/angstrom
bond_l = data["bond_l"] # unit: angstrom
angle_k = data["angle_k"] # unit: kcal/mol/radian
angle_l = data["angle_l"] # unit: theta
epsilon = data["epsilon"] # unit: kcal/mol
sigma = data["sigma"] # unit angstrom
mode = data["mode"] # mode:str, (single_homo, single_heter, sample_homo, sample_heter).
device = data["device"] # device:str (cpu, cuda)
num_sample = data["num_sample"] # num_sample:int, how many sample do you want to create.

# data file name #
model_init_save_path = data["model_init_save_path"]

#=============================#

# create distribution or not #
def create_number(num):
    if (isinstance(num,int)):
        return num
    elif(isinstance(num,float)):
        return num
    elif(isinstance(num,list)):
        sample_range = num[1] - num[0]
        sample = random.random()
        return num + sample*sample_range

    else:
        raise Exception("Not correect data type.")
    pass

# run process #
def simulate(model_init_save_path,model_eq_save_path, device="cpu"):

    # change model name in input file #
    new_file = ""
    with open("lmp_system.in","r") as file:
        for line in file:
            if "read_data" in line:     
                line = line.replace(line.split(" ")[1],model_init_save_path)
                new_file += line
            elif "write_data" in line: 
                line = line.replace(line.split(" ")[1],model_eq_save_path)   
                new_file += line
            else:
                new_file += line
    with open("lmp_system.in","w") as file:
        file.write(new_file)

    if device=="cpu":
        # run in lammps (cpu) #
        os.system("mpirun -np 4 ./lmp_mpi -in lmp_system.in")
    elif(device=="cuda"):
        # run in lammps (gpu) #
        os.system("mpirun -np 4 ./lmp_mpi -sf gpu -pk gpu 1 -in lmp_system.in")
    else:
        raise Exception("Not correct device, (cpu or cuda).")

#=============================#

if __name__ == '__main__':

    if mode=="single_homo" or mode=="single_heter":

        # parameter #
        parameter_package = [bond_k,bond_l,angle_k,angle_l,epsilon,sigma]

        if mode=="single_homo":
            # generate model #
            SAW_universe = SAW(aMW, mass, chains, bond_l)
            SAW_universe.general_model()
            SAW_universe.build_lammps_datafile(parameter_package,model_init_save_path+"homo_single.data")
            
            # equilibrium #
            simulate(model_init_save_path+"homo_single.data", model_init_save_path+"homo_single_m.data", device=device)

        elif mode=="single_heter":
            # generate model #
            SAW_universe = SAW(aMW, mass, chains, bond_l)
            SAW_universe.general_model_heter()
            SAW_universe.build_lammps_datafile_heter(parameter_package,model_init_save_path+"heter_single.data")
            
            # equilibrium #
            simulate(model_init_save_path+"heter_single.data", model_init_save_path+"heter_single_m.data", device=device)
    
    elif mode=="sample_homo" or mode=="sample_heter":
        for sample in range(num_sample):
            aMW_sample = create_number(aMW)
            mass_sample = create_number(mass)
            chains_sample = int(create_number(chains))
            bond_l_sample = create_number(bond_l)
            bond_k_sample = create_number(bond_k)
            angle_l_sample = create_number(angle_l)
            angle_k_sample = create_number(angle_k)
            epsilon_sample = create_number(epsilon)
            sigma_sample = create_number(sigma)

            # parameter #
            parameter_package = [bond_k_sample,bond_l_sample,angle_k_sample,angle_l_sample,epsilon_sample,sigma_sample]
            
            # init SAW package #
            SAW_universe = SAW(aMW_sample, mass_sample, chains_sample, bond_l_sample)

            # generate model #
            if mode=="sample_homo":
                SAW_universe.general_model()
                SAW_universe.build_lammps_datafile(parameter_package,f"{model_init_save_path}homo_sample_{sample}.data")
            
                # equilibrium #
                simulate(f"{model_init_save_path}homo_sample_{sample}.data", f"{model_init_save_path}homo_sample_m_{sample}.data", device=device)
            
            elif mode=="sample_heter":
                SAW_universe.general_model_heter()
                SAW_universe.build_lammps_datafile_heter(parameter_package,f"{model_init_save_path}heter_sample_{sample}.data")
                
                # equilibrium #
                simulate(f"{model_init_save_path}heter_sample_{sample}.data", f"{model_init_save_path}homo_sample_m_{sample}.data", device=device)
        
    else:
        raise Exception("Not correct mode.")
    print("done")
    pass