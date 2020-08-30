#! /bin/bash
#
g++ -c -Wall -I/$HOME/include midpoint_test.cpp
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
g++ -o midpoint_test midpoint_test.o /$HOME/libcpp/midpoint.o -lm
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
rm midpoint_test.o
#
./midpoint_test > midpoint_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm midpoint_test
#
gnuplot < predator_prey_midpoint_commands.txt
gnuplot < stiff_midpoint_commands.txt
#
echo "Normal end of execution."
