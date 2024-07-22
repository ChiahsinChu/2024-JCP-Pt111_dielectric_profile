from toolbox.utils import *
from toolbox.utils.utils import load_dict, safe_makedirs
from toolbox.utils.math import handle_zero_division

from intdielec.watanalysis.elec_eps import WaterEDen


dnames_int = glob.glob("/public/home/jxzhu/workspace/2022_leiden/08.Pt_111/01.elec_eps/01.calc/calc/*")
dnames_int.sort()

fnames_metal_coord = glob.glob("/public/home/jxzhu/workspace/2022_leiden/08.Pt_111/01.elec_eps/02.pure_slab/calc/task.000.*/coord.xyz") 
fnames_metal_coord.sort()
fnames_metal_eps = glob.glob("/public/home/jxzhu/workspace/2022_leiden/08.Pt_111/01.elec_eps/02.pure_slab/data/eps_data.*.pkl")
fnames_metal_eps.sort()

for dname_int, fname_metal_coord, fname_metal_eps in zip(dnames_int, fnames_metal_coord, fnames_metal_eps):
    # pure metal
    metal_data = load_dict(fname_metal_eps)[0.0]
    metal_inveps = metal_data["inveps"].mean(axis=0)
    atoms = io.read(fname_metal_coord)
    z = atoms.get_positions()[:, 2]
    z_ave_metal = np.sort(z)[-36:].mean()
    metal_grid = metal_data["v_prime_grid"] - z_ave_metal
    metal_eps = handle_zero_division(1, metal_inveps, 1e-2)
    
    # interface
    work_dir = os.path.basename(dname_int)
    safe_makedirs(work_dir)
    data_dict = load_dict(os.path.join(dname_int, "task_info.json"))
    for prefix in ["lo", "hi"]:
        z_ref = data_dict[prefix]["z_ave"]
        atoms = io.read(os.path.join(dname_int, "ref_%s/coord.xyz" % prefix))
        data = load_dict(os.path.join(dname_int, "eps_data_%s.pkl" % prefix))[0.0]
        _x = data["v_prime_grid"] - z_ref
        mask = (_x >= 0)
        inveps = data["inveps"][-1][mask]
        grid = _x[mask]
        int_eps = handle_zero_division(1, inveps, 1e-2)
        interp_metal_eps = np.interp(grid, metal_grid, metal_eps)
        wat_chi = (int_eps - 1) - (interp_metal_eps - 1)

        task = WaterEDen(atoms)
        wat_eden = task.run(grid=grid, z_ref=z_ref)
        wat_norm_chi = handle_zero_division(wat_chi, wat_eden, 1e-2)

        np.save(os.path.join(work_dir, "%s_grid.npy" % prefix), grid)
        np.save(os.path.join(work_dir, "%s_int_eps.npy" % prefix), int_eps)
        np.save(os.path.join(work_dir, "%s_wat_eden.npy" % prefix), wat_eden)
        np.save(os.path.join(work_dir, "%s_wat_chi.npy" % prefix), wat_chi)
        np.save(os.path.join(work_dir, "%s_wat_norm_chi.npy" % prefix), wat_norm_chi)
