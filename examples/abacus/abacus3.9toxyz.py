import dpdata 
from ase import Atoms 
from ase.io import write 

sys = dpdata.LabeledSystem('./', fmt='abacus/scf')
traj: list[Atoms] = sys.to_ase_structure()
print(traj[0].get_stress(voigt=True))
write("train.xyz", traj)