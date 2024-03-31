from package.SAW import SAW
#from package.MD import simulator
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
mode = data["mode"] # mode:str, (single, distribution)

# data file name #
model_init_save_path = data["model_init_save_path"]
model_melt_save_path = data["model_melt_save_path"]

#=============================#

#    check config datatype    #

if mode=="single":
    if (not isinstance(aMW,int) and not isinstance(aMW,float)):
        raise Exception("The aMW arg is wrong in single mode\.")
    if (not isinstance(mass,int) and not isinstance(mass,float)):
        raise Exception("The mass arg is wrong in single mode.")
    if (not isinstance(chains,int)):
        raise Exception("The chains arg is wrong in single mode\.")
    if (not isinstance(bond_k,int) and not isinstance(bond_k,float)):
        raise Exception("The bond_k arg is wrong in single mode.")
    if (not isinstance(bond_l,int) and not isinstance(bond_l,float)):
        raise Exception("The bond_l arg is wrong in single mode.")
    if (not isinstance(angle_k,int) and not isinstance(angle_k,float)):
        raise Exception("The angle_k arg is wrong in single mode.")
    if (not isinstance(angle_l,int) and not isinstance(angle_l,float)):
        raise Exception("The angle_l arg is wrong in single mode.")
    if (not isinstance(epsilon,int) and not isinstance(epsilon,float)):
        raise Exception("The eplison arg is wrong in single mode.")
    if (not isinstance(sigma,int) and not isinstance(sigma,float)):
        raise Exception("The sigma arg is wrong in single mode.")
    parameter_package = [bond_k,bond_l,angle_k,angle_l,epsilon,sigma]
    print("Datetype: check")
elif mode=="distribution":
    pass
else:
    raise Exception("The mode arg is wrong.")


#=============================#

if __name__ == '__main__':

    
    if mode=="single":
        # generate model #
        #SAW_universe = SAW(aMW, mass, chains, bond_l)
        #SAW_universe.general_model()
        #SAW_universe.build_lammps_datafile(parameter_package,model_init_save_path)

        # change model name in input file #
        new_file = ""
        with open("lmp_system.in","r") as file:
            for line in file:
                if "read_data" in line:     
                    line = line.replace(line.split(" ")[1],model_init_save_path)
                    new_file += line
                elif "write_data" in line: 
                    line = line.replace(line.split(" ")[1],model_melt_save_path)   
                    new_file += line
                else:
                    new_file += line
        with open("lmp_system.in","w") as file:
            file.write(new_file)

        # run in lammps #
        os.system("mpirun -np 4 ./lmp_mpi -in lmp_system.in")
        
    elif mode=="distribution":
        pass
    else:
        pass

    
    print("done")
    pass