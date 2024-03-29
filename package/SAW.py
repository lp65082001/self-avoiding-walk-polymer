import numpy as np
import random
# implement self avoiding warking (SAW) in 3D polymer #

class SAW:
    def __init__(self, aMW, mass, chains, bond_l):
        self.aMW = aMW
        self.mass = mass
        self.chains = chains
        self.bond_l = bond_l

    def decide_oriention(self):
        pass

    def general_model(self):
        # calculate number unit #
        number_unit = abs(self.aMW/self.mass)
        
        # create grid space #
        eval_grid = int(pow(number_unit * self.chains, 1/2.5)+10)
        grid_space = np.zero((eval_grid,eval_grid,eval_grid))

        # warking #

        pass

    def build_lammps_datafile(self):
        pass
