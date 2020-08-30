#! /bin/bash
#
g++ -c -Wall -I/$HOME/include lagrange_interp_2d_test.cpp
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
g++ lagrange_interp_2d_test.o /$HOME/libcpp/lagrange_interp_2d.o \
                             /$HOME/libcpp/test_interp_2d.o \
                             /$HOME/libcpp/r8lib.o -lm
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm lagrange_interp_2d_test.o
#
mv a.out lagrange_interp_2d_test
./lagrange_interp_2d_test > lagrange_interp_2d_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm lagrange_interp_2d_test
#
echo "Normal end of execution."
