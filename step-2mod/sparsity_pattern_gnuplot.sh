#! /bin/bash
#
cat <<EOF > commands.txt
set term png
set output "$1.png"
set title "$1"
set timestamp
plot "$1" with points
quit
EOF
#
gnuplot < commands.txt
echo "gnuplot created $1.png"
#
echo "Normal end of execution."

