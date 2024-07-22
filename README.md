# README

## to-do list

- [ ] 01.dpgen
- [ ] 03.mlp_validation
    - [ ] PES
- [ ] 05.elec_eps
- [ ] 06.orient_eps

- [ ] add `requirements.txt` for python codes (`pipreqs myproj/`)
 
[info of publication]

## `00.ffmd`

(Optional) Input for performing FFMD simulation of Pt(111)/water interfaces for initial training dataset.

``` bash
cd 00.ffmd
# initial config for Pt(111)/water interface
python make_config.py
# setup lsf / slurm script according to your cluster
bash run.sh
```

## `01.dpgen`

- version: `0.11.0`

```bash
dpgen run param.json machine.json
```


## `02.mlp_train`

## `03.mlp_validation`

placeholder

## `04.mlmd`

The ML-MD simulations were performed with modified cp2k based on v9.0 (development version), which can be installed with the following commands [[ref]](https://wiki.cheng-group.net/wiki/software_installation/deepmd-kit/deepmd-kit_installation_191/#dp-cp2k):

``` bash
git clone https://github.com/Cloudac7/cp2k.git -b deepmd_latest --recursive --depth=1
cd cp2k/tools/toolchain/
./install_cp2k_toolchain.sh --enable-cuda=no --with-deepmd=$deepmd_root --with-tfcc=$tensorflow_root --deepmd-mode=cuda --mpi-mode=no --with-libint=no --with-libxc=no --with-libxsmm=no
make -j 4 ARCH=local VERSION="ssmp sdbg"
```
The codes can also be found in `softwares/dp2.0-cp2k.tar.gz`.

The simulations can be run with the following commands:

``` bash
module add deepmd/2.0
# only serial version is available!
cp2k.sopt input.inp
```

The output trajectory file is in `xyz` format, which can be converted to `xtc` format with the `xyz_to_xtc.py` script for further analysis.


## `05.elec_eps`

## `06.orient_eps`

## `07.manuscript`

## `trajs`

- `aimd.xyz`: AIMD trajectory of Pt(111)/water interface

## `softwares`

- `dp2.0-cp2k.tar.gz`: Deep Potential training code

