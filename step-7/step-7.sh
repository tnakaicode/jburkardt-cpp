#! /bin/bash
#
## step-7.sh
#
#  Discussion:
#
#    Set up and run the step-7 example.
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
echo "step-7.sh"
echo "  Set up and run the step-7 example on my desktop running Linux"
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
./step-7 &> step-7.txt
#
#  5: Process the latex files.
#
latex convergence-global-q1.tex
dvipdf convergence-global-q1.dvi
latex convergence-global-q2.tex
dvipdf convergence-global-q2.dvi
latex error-adaptive-q1.tex
dvipdf error-adaptive-q1.dvi
latex error-adaptive-q2.tex
dvipdf error-adaptive-q2.dvi
latex error-global-q1.tex
dvipdf error-global-q1.dvi
latex error-global-q2.tex
dvipdf error-global-q2.dvi
#
rm *.aux
rm *.dvi
rm *.log
make clean
#
#  Terminate.
#
echo ""
echo "step-7.sh:"
echo "  Normal end of execution."
 
