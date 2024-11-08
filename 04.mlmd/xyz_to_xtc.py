"""
Ref:
- https://userguide.mdanalysis.org/stable/formats/index.html
- https://manual.gromacs.org/current/reference-manual/file-formats.html#xtc
"""
import MDAnalysis as mda

from toolbox.utils.utils import load_dict


u = mda.Universe("coord.xyz", 
                 "dpmd-pos-1.xyz", 
                 topology_format="XYZ", 
                 format="XYZ",
                 )
with mda.Writer("dump.xtc") as W:
    for ts in u.trajectory:
        W.write(u)
        
