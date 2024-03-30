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

    # create random initial coordination #
    def random_init_coordination(self):
        non_use = np.where(self.grid_space==0)
        rand_init = non_use[int(random.random() * non_use.shape[0]),:]
        return rand_init

    # random decide orientation #
    def random_decide_orientation(self,current_coordination):
        xmin = current_coordination[0] - 1
        ymin = current_coordination[1] - 1
        zmin = current_coordination[2] - 1
        all_candidate = []
        for i in [xmin,xmin+1,xmin+2]:
            for j in [ymin,ymin+1,ymin+2]:
                for k in [zmin,zmin+1,zmin+2]:
                    all_candidate.append([i,j,k])
        all_candidate  = np.array(all_candidate)
        
        non_use_cand = []
        for check_index in all_candidate:
            if self.grid_space[check_index] != 0:
                non_use_cand.append(check_index)
        
        rand_orientation = non_use_cand[int(random.random() * len(non_use_cand)),:]

        return rand_orientation
    # model generation #
    def general_model(self):
        # calculate number unit #
        number_unit = abs(self.aMW/self.mass)
        
        # create grid space #
        eval_grid = int(pow(number_unit * self.chains, 1/2.5)+10)
        self.grid_space = np.zeros((eval_grid,eval_grid,eval_grid)).astype(int)

        # particule table #
        self.particule_table = np.zeros((number_unit*self.chains,4)).astype(int)

        # warking #
        num_index = 0 
        for num_chains in trange(self.chains):
            # create initial coordination #
            init_coordination = self.random_init_coordination()
            self.particule_table[num_index,:] = [init_coordination,num_chains]
            next_coordination = init_coordination
            num_index += 1
            for num_particule in range(number_unit-1):
                new_next_coordination = self.random_decide_orientation(next_coordination)
                self.particule_table[num_index,:] = [new_next_coordination,num_chains]
                next_coordination = new_next_coordination
                num_index += 1
                
        # rescaling grid space to real space #
        self.particule_table = self.particule_table.astype(float)
        self.particule_table[:,0:3] = self.particule_table[:,0:3] * self.bond_l
        self.init_box_size = eval_grid * self.bond_l

    # write the model into lammps data format #
    def build_lammps_datafile(self,parameter,path = None):
        if path == None:
            path = "./test.data"

        #header#
        f = open(path+".data", 'w')
        f.write("create by Amborse hui from M^5 lab\n")
        f.write("\n")
        f.write(str(self.res*self.sig)+" atoms\n")
        f.write("1 atom types\n")
        f.write(str(self.sig*(self.res-1))+" bonds\n")
        f.write("1 bond types\n")
        f.write(str(self.sig*(self.res-2))+" angles\n")
        f.write("1 angle types\n")
        f.write(str(self.sig*(self.res-3))+" dihedrals\n")
        f.write("1 dihedral types\n")
        f.write("0 impropers\n")
        f.write("0 improper types\n")
        f.write("\n")
        f.write("0.0 "+str(self.box[0])+" xlo xhi\n")
        f.write("0.0 "+str(self.box[1])+" ylo yhi\n") 
        f.write("0.0 "+str(self.box[2])+" zlo zhi\n")  
        f.write("\n")
        f.write("Masses\n")
        f.write("\n")
        f.write("1 "+str(self.bead_mass)+" \n")
        f.write("\n")
        f.write("Pair Coeffs # lj/cut\n")
        f.write("\n")
        f.write("1 "+str(self.epsilon_)+" "+str(self.sigma_)+"\n")
        f.write("\n")
        f.write("Bond Coeffs # harmonic\n")
        f.write("\n")
        f.write("1 "+str(self.bondk)+" "+str(self.bondi)+"\n")
        f.write("\n")
        f.write("Angle Coeffs # harmonic\n")
        f.write("\n")
        f.write("1 "+str(self.anglek[0])+" "+str(self.anglei)+"\n")
        f.write("\n")

        # Atom list
        f.write("Atoms #molecule \n")
        f.write("\n")
        for i in range(0,self.list_position.shape[1]):
            f.write(str(i+1)+" 1 1 "+str(self.get_pos()[-1][i][0])+" "+str(self.get_pos()[-1][i][1])+" "+str(self.get_pos()[-1][i][2])+" 0 0 0\n")
        f.write("\n")

        # Bond list
        f.write("Bonds\n")
        f.write("\n")
        n = 1
        for j in range(0,self.sig):
            for i in range(0,self.res-1):
                f.write(str(n)+" 1 "+str(i+1+j*(self.res))+" "+str(i+2+j*(self.res))+"\n")
                n += 1
        f.write("\n")

        # angle list
        f.write("Angles\n")
        f.write("\n")
        n = 1
        for j in range(0,self.sig):
            for i in range(0,self.res-2):
                f.write(str(n)+" 1 "+str(i+1+j*(self.res))+" "+str(i+2+j*(self.res))+" "+str(i+3+j*(self.res))+"\n")
                n += 1
        f.write("\n")

