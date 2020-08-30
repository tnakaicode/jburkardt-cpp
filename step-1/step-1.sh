#! /bin/bash
#
## step-1.sh compiles and runs the step-1 dealii example.
#
#  Discussion:
#
#    Before doing anything with deal.ii you have to:
#    1) ensure that "cmake" is in your path.
#       To do this, you could type "which cmake".  For me, the response is
#       /usr/bin/cmake
#       so I know cmake is installed.  
#    2) download and install deal.ii.  On my local system, I can't
#       install anything, but I do know where the system administrator
#       installed deal.ii:
#       /usr/common/deal.II
#
#    Then, you might want to run one of the deal.ii examples,
#    such as the first one, called step-1.  Rather than running
#    this example in the examples directory provided by deal.ii, namely:
#    /usr/common/deal.II/dealii-8.2.1/examples/step-1
#    I much prefer to make my own directory.  Then I must copy
#    from the deal.ii examples/step-1 directory the files
#      CMakeLists.txt
#    and
#      step-1.cc
#
#    If I can get step-1.cc to compile, link and run, then
#    I can probably understand how to do the same for a code
#    I write myself.
#
#    Assuming this directory already has the IMPORTANT files
#      CMakeLists.txt and step-1.cc
#    the following command will create Makefile (and much more!).
#    The files which are created are "UNIMPORTANT", in that they
#    can always be regenerated.  They are disgusting, because they
#    clutter up your directory, they are mysterious, and they are "fragile".
#    On the other hand, they do make an unimaginably complex program
#    accessible.  This is a price we pay these days.
#
#    Another note: if for any reason you MOVE or RENAME this directory,
#    then the files created by cmake will become "confused"...and worse,
#    may not work properly.  Moreover, this script itself may then fail,
#    because one of the imaginary virtues of makefiles is that they
#    skip any task that seems unnecessary, so even though the files
#    created by cmake have now become unusable, the make process
#    only checks that they already exist and proceeds to fall flat on
#    its face, instead of regenerating them.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    22 May 2020
#
#  Author:
#
#    John Burkardt
#
echo ""
echo "step-1.sh"
echo "  A record of the basic commands necessary to set up the"
echo "  step-1 example, and then to run it, on my desktop running Linux"
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
./step-1 &> step-1.txt
#
#  5: Convert PS to PNG
#
echo ""
echo "4) Convert useless PS files to PNG."
echo ""
ps2png grid-1.eps grid-1.png
ps2png grid-2.eps grid-2.png
#
make clean
#
#  Terminate.
#
echo ""
echo "step-1.sh:"
echo "  Normal end of execution."
 
