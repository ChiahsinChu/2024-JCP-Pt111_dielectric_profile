from toolbox.utils import *
from toolbox.calculator import CP2KDPDispatcher


task = CP2KDPDispatcher(work_dir=".")
machine_setup = {
    "remote_root": "/public/home/jxzhu/tmp_calc"
}
resources_setup = {
    "queue_name": "small_s",
    "number_node": 1,
    "cpu_per_node": 32,
    "group_size": 50,
    "module_list": [
        "gcc/9.3",
        "intel/2020.2",
        "cp2k/2022.1-intel-2020"
    ],
    "envs": {
        "OMP_NUM_THREADS": 1
    }
}
task_setup = {
    "command": "mpirun cp2k.psmp -i input.inp",
    "backward_files": [
        "output", "cp2k-RESTART.wfn",
        "cp2k-TOTAL_DENSITY-1_0.cube", "cp2k-v_hartree-1_0.cube",
        "cp2k-ELECTRON_DENSITY-1_0.cube"
    ]
}

work_dirs = glob.glob("calc/ref.*")
task.run(work_dirs, 
         machine_setup=machine_setup, 
         resources_setup=resources_setup, 
         task_setup=task_setup)

work_dirs = glob.glob("calc/task.*")
task.run(work_dirs, 
         machine_setup=machine_setup, 
         resources_setup=resources_setup, 
         task_setup=task_setup)
