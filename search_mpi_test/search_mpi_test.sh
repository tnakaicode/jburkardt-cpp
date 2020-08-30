#! /bin/bash
#
mpirun -np 8 $HOME/bincpp/search_mpi > search_mpi_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Normal end of execution."

