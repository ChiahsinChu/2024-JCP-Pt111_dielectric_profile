import os
import glob

import numpy as np
from tqdm import trange
from deepmd.infer import DeepPot


dp = DeepPot("graph.pb")
type_map = dp.tmap

def make_atype(dname):
    _type_map = np.loadtxt(os.path.join(dname, "type_map.raw"), dtype=str)
    _atype = np.loadtxt(os.path.join(dname, "type.raw"), dtype=np.int32)
    symbols = _type_map[_atype]
    
    atype = np.ones_like(_atype, dtype=np.int32)
    for ii, _type in enumerate(type_map):
        atype[symbols == _type] = ii
    return atype

def dp_test(dname):
    print(dname)
    atype = make_atype(dname)
    # print(atype)
    
    fnames = glob.glob(os.path.join(dname, "*/coord.npy"))
    fnames.sort()
    for fname in fnames:
        coords = np.load(fname)
        cells = np.load(fname.replace("coord", "box"))
        nframes = len(coords)
        energies = []
        forces = []
        for ii in trange(nframes):
            out = dp.eval(coords[ii].reshape(1, -1), cells[ii].reshape(1, -1), atype)
            energies.append(out[0].reshape(-1))
            forces.append(out[1].reshape(-1))
        np.save(fname.replace("coord", "energy"), energies)
        np.save(fname.replace("coord", "force"), forces)
        
if __name__ == '__main__':
    fnames = glob.glob("ml_data/**/type.raw", recursive=True)
    fnames.sort()
    for fname in fnames:
        dname = os.path.dirname(fname)
        dp_test(dname)