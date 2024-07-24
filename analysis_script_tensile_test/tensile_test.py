import os 
import time
import numpy as np
import json
from os import listdir
from os.path import isfile, isdir, join


def tensile_test(model_init_save_path):
    new_file = ""
    with open("lmp_tensile_test.in","r") as file:
            for line in file:
                if "read_data" in line:     
                    line = line.replace(line.split(" ")[1],model_init_save_path)
                    new_file += line
                else:
                    new_file += line
    with open("lmp_tensile_test.in","w") as file:
        file.write(new_file)

    os.system("mpirun -np 4 ./lmp_mpi -sf gpu -pk gpu 0 -in lmp_tensile_test.in")
    print("================ Done ================")
    time.sleep(1)
    with open("log.lammps","r") as f:
        data = f.read()

    data = data.split("\n")

    get_data = []
    sw = 0
    for i in range(len(data)):
        if "Loop time" in data[i]: 
            sw = 0
        if sw == 1:
            get_data.append(list(filter(None,data[i].split(" "))))

        if "Step" in data[i] and "CPU" in data[i]:
            sw = 1
    os.system("rm log.lammps")
    return np.array(get_data).astype(float)

# calculate young's modulus (unit MPa)#
def young_modulus(data,stress=30): # 30 MPa
    strain = (data[-1,2] - data[0,2])/data[0,2]
    return stress/strain


if __name__ == "__main__":
    mypath = "./datafile/"
    files = listdir(mypath)
    data_all = dict()
    # data_field: step temp lx ly lz vol cpu #
    for f in files:
        fullpath = join(mypath, f)
        data_field = tensile_test(fullpath)

        # calculate young's modulus under 30 MPa #
        ym = young_modulus(data_field)

        """ calculate other properties can add here"""
        data_all[fullpath] = {"ym":ym}

    with open('properties.json', 'w') as f:
        json.dump(data_all, f)
