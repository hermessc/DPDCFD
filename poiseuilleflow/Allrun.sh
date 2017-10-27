set -x 
blockMesh
decomposePar
export LD_LIBRARY_PATH=/home/hermes/Desktop/Programmi/Lammps/lammps-14May16/src:$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/hermes/Desktop/Programmi/MPI/install/bin:$LD_LIBRARY_PATH
$HOME/Desktop/Programmi/MPI/install/bin/mpirun -n 2 AnonNewtonianIcoFoam -parallel >& log
reconstructPar

