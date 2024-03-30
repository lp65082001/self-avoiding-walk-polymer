from lammps import lammps


class simulator:
    def __init__(self,parameter):
        self.mass = parameter[0]
        self.bond_k = parameter[1]
        self.bond_l = parameter[2]
        self.angle_k = parameter[3]
        self.angle_l = parameter[4]
        self.epsilon = parameter[5]
        self.sigma = parameter[6]

    # run cooling # 
    def run(self,model_path,simulate_script="./lmp_system.in"):
        save_path = model_path.split(".")[0]+"_m.data"
        lmp = lammps()
        # load system setting # 
        lmp.file(simulate_script)
        block = f"""
#read_data	{model_path}

#minimize 
minimize 0.0 1.0e-8 1000 100000

# time integration settings
timestep		20	# 20fs

#thermosystem set
velocity 		all create 500.0 4928459 dist gaussian
  
fix				2 all npt temp 500.0 300.0 200 iso 1.0 1.0 1000
"""
        lmp.commands_string(block)
        lmp.command("run 100000 pre no post no")
        lmp.command(f"write_data {save_path}")
        lmp.close()
