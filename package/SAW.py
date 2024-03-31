import random
import numpy as np
from tqdm import trange
# implement self avoiding warking (SAW) in 3D polymer #

class SAW:
    def __init__(self, aMW, mass, chains, bond_l):
        self.aMW = aMW
        self.mass = mass
        self.chains = chains
        self.bond_l = bond_l

    # check and modify under PBC
    def PBC_check(self,x):
        if x >= self.grid_space.shape[0]:
            x = x - self.grid_space.shape[0]
        elif x < 0:
            x = x + self.grid_space.shape[0]
        return x

    # create random initial coordination #
    def random_init_coordination(self):
        non_use = np.where(self.grid_space==0)
        non_use_ = np.hstack((np.hstack((non_use[0].reshape(-1,1),non_use[1].reshape(-1,1))),non_use[2].reshape(-1,1)))
        rand_init = non_use_[int(random.random() * non_use_.shape[0]),:]
        return rand_init

    # unwrap the model #
    def unwrap_model(self):
        for i in range(self.chains*(self.res-1)):
            if self.particule_table[i+1,3]==self.particule_table[i,3]:
                check_point = 0
                vector = self.particule_table[i+1] - self.particule_table[i]
                if abs(vector[0])>=self.grid_space.shape[0]-1 or abs(vector[1])>=self.grid_space.shape[0]-1 or abs(vector[2])>=self.grid_space.shape[0]-1:
                    if vector[0] >= self.grid_space.shape[0]-1:
                        self.particule_table[i+1,0] = self.particule_table[i+1,0] - self.grid_space.shape[0]
                        check_point = 1
                    elif(vector[0] <= -1*self.grid_space.shape[0]+1):
                        self.particule_table[i+1,0] = self.particule_table[i+1,0] + self.grid_space.shape[0]
                        check_point = 1
                    if vector[1] >= self.grid_space.shape[0]-1:
                        self.particule_table[i+1,1] = self.particule_table[i+1,1] - self.grid_space.shape[0]
                        check_point = 1
                    elif(vector[1] <= -1*self.grid_space.shape[0]+1):
                        self.particule_table[i+1,1] = self.particule_table[i+1,1] + self.grid_space.shape[0]
                        check_point = 1
                    if vector[2] >= self.grid_space.shape[0]-1:
                        self.particule_table[i+1,2] = self.particule_table[i+1,2] - self.grid_space.shape[0]
                        check_point = 1
                    elif(vector[2] <= -1*self.grid_space.shape[0]+1):
                        self.particule_table[i+1,2] = self.particule_table[i+1,2] + self.grid_space.shape[0]
                        check_point = 1
            
            # unwrap bug #
            if np.linalg.norm(self.particule_table[i+1,0:3] - self.particule_table[i,0:3])>self.grid_space.shape[0]:
                vector_ = self.particule_table[i+1,0:3] - self.particule_table[i,0:3]
                vector_num = np.linalg.norm(vector)
                #self.particule_table[i+1,0] = self.particule_table[i+1,0]-(vector_[0]/vector_num)*self.grid_space.shape[0]
                #self.particule_table[i+1,1] = self.particule_table[i+1,1]-(vector_[1]/vector_num)*self.grid_space.shape[0]
                #self.particule_table[i+1,2] = self.particule_table[i+1,2]-(vector_[2]/vector_num)*self.grid_space.shape[0]
                print(vector_)
                print(vector_num)
                print(check_point)
                # have bug
        print("unwrap done")

    # random decide orientation #
    def random_decide_orientation(self,current_coordination):
        x_current = current_coordination[0]
        y_current = current_coordination[1]
        z_current = current_coordination[2]
        all_candidate = []
        for i in [x_current-1,x_current,x_current+1]:
            for j in [y_current-1,y_current,y_current+1]:
                for k in [z_current-1,z_current,z_current+1]:
                    all_candidate.append([self.PBC_check(i),self.PBC_check(j),self.PBC_check(k)])
        all_candidate  = np.array(all_candidate)
        
        non_use_cand = []
        for i in range(all_candidate.shape[0]):
            check_index_x = all_candidate[i,0]
            check_index_y = all_candidate[i,1]
            check_index_z = all_candidate[i,2]
            if self.grid_space[check_index_x,check_index_y,check_index_z] == 0:
                non_use_cand.append([check_index_x,check_index_y,check_index_z])
        rand_select = int(random.random() * len(non_use_cand))
        rand_orientation = non_use_cand[rand_select]
  
        return rand_orientation
        
    # model generation #
    def general_model(self):
        # calculate number unit #
        number_unit = int(abs(self.aMW/self.mass))
        self.res = number_unit

        # create grid space #
        eval_grid = int(pow(number_unit * self.chains, 1/3)+100)
        self.grid_space = np.zeros((eval_grid,eval_grid,eval_grid)).astype(int)

        # particule table #
        self.particule_table = np.zeros((number_unit*self.chains,4)).astype(int)

        # warking #
        num_index = 0 
        for num_chains in trange(self.chains):
            # create initial coordination #
            init_coordination = self.random_init_coordination()
            self.particule_table[num_index,0:3] = init_coordination
            self.particule_table[num_index,3] = num_chains
            self.grid_space[init_coordination[0],init_coordination[1],init_coordination[2]] = 1
            next_coordination = init_coordination
            num_index += 1
            for num_particule in range(number_unit-1):
                new_next_coordination = self.random_decide_orientation(next_coordination)
                self.particule_table[num_index,0:3] = new_next_coordination
                self.particule_table[num_index,3] = num_chains
                self.grid_space[new_next_coordination[0],new_next_coordination[1],new_next_coordination[2]] = 1
                next_coordination = new_next_coordination
                num_index += 1

        # unwrap model #
        self.unwrap_model()

        # rescaling grid space to real space #
        self.particule_table = self.particule_table.astype(float)
        self.particule_table[:,0:3] = self.particule_table[:,0:3] * self.bond_l
        self.init_box_size = eval_grid * self.bond_l

    # write the model into lammps data format #
    def build_lammps_datafile(self,parameter,path = None):
        if path == None:
            path = "./test.data"

        # bond and angle table #
        self.bond_num = self.chains * (self.res - 1)
        self.angle_num = self.chains * (self.res - 2)

        #header#
        f = open(f"{path}", 'w')
        f.write("Created by Amborse hui from M^5 lab\n")
        f.write("\n")
        f.write(f"{self.particule_table.shape[0]} atoms\n")
        f.write("1 atom types\n")
        f.write(f"{self.bond_num} bonds\n")
        f.write("1 bond types\n")
        f.write(f"{self.angle_num} angles\n")
        f.write("1 angle types\n")
        f.write("\n")
        f.write(f"0.0 {self.init_box_size} xlo xhi\n")
        f.write(f"0.0 {self.init_box_size} ylo yhi\n") 
        f.write(f"0.0 {self.init_box_size} zlo zhi\n")  
        f.write("\n")
        f.write("Masses\n")
        f.write("\n")
        f.write(f"1 {self.mass} \n")
        f.write("\n")
        f.write("Pair Coeffs # lj/cut\n")
        f.write("\n")
        f.write(f"1 {parameter[4]} {parameter[5]}\n")
        f.write("\n")
        f.write("Bond Coeffs # harmonic\n")
        f.write("\n")
        f.write(f"1 {parameter[0]} {parameter[1]}\n")
        f.write("\n")
        f.write("Angle Coeffs # harmonic\n")
        f.write("\n")
        f.write(f"1 {parameter[2]} {parameter[3]}\n")
        f.write("\n")

        # Atom list
        f.write("Atoms #angle \n")
        f.write("\n")
        for i in range(0,self.particule_table.shape[0]):
            f.write(str(i+1)+" 1 1 "+str(self.particule_table[i,0])+" "+str(self.particule_table[i,1])+" "+str(self.particule_table[i,2])+" 0 0 0\n")
        f.write("\n")

        # Bond list
        f.write("Bonds\n")
        f.write("\n")
        n = 1
        for j in range(0,self.chains):
            for i in range(0,self.res-1):
                f.write(str(n)+" 1 "+str(i+1+j*(self.res))+" "+str(i+2+j*(self.res))+"\n")
                n += 1
        f.write("\n")

        # angle list
        f.write("Angles\n")
        f.write("\n")
        n = 1
        for j in range(0,self.chains):
            for i in range(0,self.res-2):
                f.write(str(n)+" 1 "+str(i+1+j*(self.res))+" "+str(i+2+j*(self.res))+" "+str(i+3+j*(self.res))+"\n")
                n += 1
        f.write("\n")
        f.close()

