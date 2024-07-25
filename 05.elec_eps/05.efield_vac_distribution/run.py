import glob
import numpy as np

from scipy import stats

from toolbox.utils.utils import load_dict

fnames = glob.glob("../01.eps_calc/calc/*/eps_data_*.pkl")
fnames.sort()
efield = []
for fname in fnames:
    data = load_dict(fname)[0.0]
    grid = data["v_grid"]
    mask = (grid > 50.0) & (grid < 60.0)
    for hartree in data["hartree"][-3:]:
        lsq_out = stats.linregress(grid[mask], hartree[mask])
        # print(lsq_out.rvalue)
        efield.append(lsq_out.slope)

efield = np.reshape(efield, [-1, 3])
ids = np.argmin(np.abs(efield), axis=-1)
# take efield at axis=1 with ids
efield_out = efield[np.arange(efield.shape[0]), ids]

print("E_vac=%.3f+-%.3f [V/Å]" % (np.mean(efield_out), np.std(efield_out)))
