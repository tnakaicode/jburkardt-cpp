#! /bin/bash
#
g++ -c -Wall -I/$HOME/include brent_test.cpp
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
g++ brent_test.o /$HOME/libcpp/brent.o -lm
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm brent_test.o
#
mv a.out brent_test
./brent_test > brent_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
rm brent_test
#
echo "Normal end of execution."
