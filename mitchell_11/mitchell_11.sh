#! /bin/bash
#
## mitchell_11.sh compiles and runs the mitchell_11 dealii example.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    23 May 2020
#
#  Author:
#
#    John Burkardt
#
echo ""
echo "mitchell_11.sh"
echo "  Run the mitchell_11 dealii program." 
#
#  1: As a precaution, remove old files.
#
echo ""
echo "1) Remove old files created by make."
echo ""
rm CMakeCache.txt
rm -r CMakeFiles
rm Makefile
rm cmake_install.cmake
#
#  2: Run cmake.
#
echo ""
echo "2) Use cmake to create the Makefile."
echo ""
cmake -DDEAL_II_DIR=/usr/local/dealii .
#
#  3: Compile and run the program.
#
echo ""
echo "3) Use 'make run' to compile, link and run the mitchell_11 example."
echo ""
make
#
#  4: Run the program.
#
./mitchell_11 &> mitchell_11.txt
#
make clean
#
#  Terminate.
#
echo ""
echo "mitchell_11.sh:"
echo "  Normal end of execution."
 
