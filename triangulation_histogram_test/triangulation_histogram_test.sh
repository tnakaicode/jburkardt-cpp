#! /bin/bash
#
$HOME/bincpp/triangulation_histogram house3 thousand_xy.txt > triangulation_histogram_test.txt
if [ $? -ne 0 ]; then
  echo "Run error."
  exit
fi
#
echo "Normal end of execution."
