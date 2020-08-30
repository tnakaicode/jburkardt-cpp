#! /bin/bash
#
cp midpoint.hpp /$HOME/include
#
g++ -c -Wall -I/$HOME/include midpoint.cpp
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv midpoint.o ~/libcpp/midpoint.o
#
echo "Normal end of execution."
