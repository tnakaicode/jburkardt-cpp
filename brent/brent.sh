#! /bin/bash
#
cp brent.hpp /$HOME/include
#
g++ -c -Wall -I /$HOME/include brent.cpp
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
mv brent.o ~/libcpp/brent.o
#
echo "Normal end of execution."
