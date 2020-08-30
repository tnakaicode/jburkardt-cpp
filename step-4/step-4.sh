#! /bin/bash
#
## step-4.sh
#
#  Discussion:
#
#    Set up and run the step-4 example.
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
echo "step-4.sh"
echo "  A record of the basic commands necessary to set up the"
echo "  step-4 example, and then to run it, on my desktop running Linux"
echo "  and using a common file system." 
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
#  3: Compile the program.
#
echo ""
echo "3) Use 'make run' to compile, link and run the step-1 example."
echo ""
make
#
#  4: Run the program.
#
./step-4 &> step-4.txt
#
make clean
#
#  Terminate.
#
echo ""
echo "STEP-4.SH:"
echo "  Normal end of execution."
 
