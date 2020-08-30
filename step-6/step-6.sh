#! /bin/bash
#
## step-6.sh
#
#  Discussion:
#
#    Set up and run the step-6 example.
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
echo "step-6.sh"
echo "  Set up and run the step-6 example on my desktop running Linux"
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
./step-6 &> step-6.txt
#
#  5: Convert PS to PNG
#
ps2png grid-0.eps grid-0.png
ps2png grid-1.eps grid-1.png
ps2png grid-2.eps grid-2.png
ps2png grid-3.eps grid-3.png
ps2png grid-4.eps grid-4.png
ps2png grid-5.eps grid-5.png
ps2png grid-6.eps grid-6.png
ps2png grid-7.eps grid-7.png
rm *.eps
#
make clean
#
#  Terminate.
#
echo ""
echo "step-6.sh:"
echo "  Normal end of execution."
 
