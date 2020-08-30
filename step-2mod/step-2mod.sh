#! /bin/bash
#
## step-2mod.sh
#
#  Discussion:
#
#    Set up and run the step-2mod example.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    24 May 2020
#
#  Author:
#
#    John Burkardt
#
echo ""
echo "step-2mod.sh"
echo "  A record of the basic commands necessary to set up the"
echo "  step-2mod example, and then to run it, on my desktop running Linux"
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
#  3: Compile and run the program.
#
echo ""
echo "3) Use 'make run' to compile, link and run the step-1 example."
echo ""
make
#
#  4: Run the program.
#
./step-2mod &> step-2mod.txt
#
#  5: Use GNUPLOT script to turn sparsity patterns into plots.
#
echo ""
echo "4) Convert sparsity patterns to PNG files."
echo ""
./sparsity_pattern_gnuplot.sh sparsity_pattern.1
./sparsity_pattern_gnuplot.sh sparsity_pattern.2
#
make clean
#
#  Terminate.
#
echo ""
echo "STEP-2MOD.SH:"
echo "  Normal end of execution."
 
