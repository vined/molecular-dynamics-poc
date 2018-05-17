# molecular-dynamics-poc
Molecular dynamics simulations for master thesis poc


## Dependencies

- VMD - for visulisation


## Build

Process is usual to CMake based projects:

```bash
mkdir build
cd build
cmake ..
make
```


## Running

In the same directory run:
```bash
./molecular_dynamics_poc ../params.txt ../data/atoms.xyz
```
This will write simulation output to _out_ directory of parent folder.


## Viewing the results

Results can be viewed as plots or simulation video.

To view simulation result plots use:
```bash
./plot.sh
```

To view simulation video VMD is requierd. To show results in VMD use:
```bash
vmd -nt -e vmd/load.tcl 
```


## TODO
- Initialize random velocities
- Add other forces (electrostatic, angles and ...) and modify to allow use of different atoms
- Add neighbours if box is bigger than cut-off*4 radius
- Add velocity scaling (equilibration)
- Change to fluctuating charge force field
- Add MPI
- Water polarized model (hitch-13)
