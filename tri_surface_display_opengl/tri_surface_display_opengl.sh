#! /bin/bash
#
g++ -c -Wall tri_surface_display_opengl.cpp
if [ $? -ne 0 ]; then
  echo "Compile error."
  exit
fi
#
#  Here is the load statement for Apple's OS X.
#
#g++ tri_surface_display_opengl.o -lm -framework GLUT -framework OpenGL
#
#  Here is the load statement for a normal UNIX system!
#
g++ tri_surface_display_opengl.o -lm -lGL -lGLU -lglut
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm tri_surface_display_opengl.o
mv a.out ~/bincpp/tri_surface_display_opengl
#
echo "Normal end of execution."
