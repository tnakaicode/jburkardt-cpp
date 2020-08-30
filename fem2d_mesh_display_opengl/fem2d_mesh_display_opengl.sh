#! /bin/bash
#
g++ -c -Wall fem2d_mesh_display_opengl.cpp
if [ $? -ne 0 ]; then
  echo "Compile error"
  exit
fi
#
g++ fem2d_mesh_display_opengl.o -lm -lGL -lGLU -lglut
if [ $? -ne 0 ]; then
  echo "Load error."
  exit
fi
#
rm fem2d_mesh_display_opengl.o
mv a.out ~/bincpp/fem2d_mesh_display_opengl
#
echo "Normal end of execution."
