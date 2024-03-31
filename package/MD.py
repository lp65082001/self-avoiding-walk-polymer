from lammps import lammps

class simulator:
    def __init__(self,parameter):
        self.bond_k = parameter[0]
        self.bond_l = parameter[1]
        self.angle_k = parameter[2]
        self.angle_l = parameter[3]
        self.epsilon = parameter[4]
        self.sigma = parameter[5]

    # run cooling # 
    def run(self,model_path,simulate_script="./lmp_system.in"):
        save_path = model_path.split(".data")[0]+"_m.data"
        lmp = lammps()
        # load system setting # 
        lmp.file(simulate_script)
        block = f"""
read_data	{model_path}

#minimize 
minimize 0.0 1.0e-8 1000 100000

# time integration settings #
timestep		20	# 20fs

# thermosystem set #
velocity 		all create 500.0 4928459 dist gaussian
  
# reset timestep #
reset_timestep 0

# run cooling #
thermo 10000
fix				1 all npt temp 500.0 500.0 200 iso 1.0 1.0 1000
run 50000 
write_data {save_path}
"""
        lmp.commands_string(block)
        lmp.close()
    