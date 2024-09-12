# README

- `*aimd_prl*` for AIMD reference data while `*mlmd*` for MLMD data.
- `water_structure.*.pkl`: files for water structures, which can be read as python `Dict`. The dict contains the following keys:
  - `rho_water`: [grid, density of water] 2D array
  - `geo_dipole_water`: [grid, dipole of water] 2D array
  - `water_a_theta`: [grid, angle between chemisorbed water dipole and z-axis] 2D array
  - `water_b_theta`: [grid, angle between physisorbed water dipole and z-axis] 2D array
  - `water_a_coverage`: [chemisorbed water coverage at every timestep] 1D array
  - `water_b_coverage`: [physisorbed water coverage at every timestep] 1D array
- `water_sp*.pkl`: files for water survival probability, which can be read as python `Dict`. The dict contains the following keys:
  - `water_a`: [time interval, survival probability of chemisorbed water] 2D array
  - `water_a_and_b`: [time interval, survival probability of chemisorbed + physisorbed water] 2D array
- `hb.*.npy`: files for water hydrogen bonds number, which can be read as 2D `numpy.ndarray` [time, number of hydrogen bonds, time cummulation of hydrogen bonds]. `wata` and `watb` in the file names refer to chemisorbed and physisorbed water, respectively. `a` and `d` refer to acceptor and donor, respectively.
