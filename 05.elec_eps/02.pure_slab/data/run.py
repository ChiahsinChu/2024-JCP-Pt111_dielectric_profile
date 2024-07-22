from toolbox.utils import *

from intdielec.workflow.elec_eps import ElecEps


v_seq = np.arange(-1.0, 1.1, 1.0)

dnames = glob.glob("../calc/ref.*/")
dnames.sort()
atoms = io.read(os.path.join(dnames[0], "coord.xyz"))
wf = ElecEps(atoms, ".", v_seq=v_seq)
wf.v_tasks = dnames
wf.calculate(pos_vac=10., save_fname="eps_data_ref")

_dnames = glob.glob("../calc/task.*/")
_dnames.sort()
_dnames = np.reshape(_dnames, [3, -1])
_dnames = np.transpose(_dnames)
for ii, dnames in enumerate(_dnames):
    atoms = io.read(os.path.join(dnames[0], "coord.xyz"))
    wf = ElecEps(atoms, ".", v_seq=v_seq)
    wf.v_tasks = dnames
    wf.calculate(pos_vac=10., save_fname="eps_data.%03d" % ii)