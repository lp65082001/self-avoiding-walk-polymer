from package.SAW import model_build
from package.MD import model_universe
import json

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
eplison = data["eplison"] # unit: kcal/mol
sigma = data["sigma"] # unit angstrom



if __name__ == '__main__':
    pass