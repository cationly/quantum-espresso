#!/bin/bash

for dist in 1.2 1.3 1.4 1.5 1.6 1.7; do

cat <<EOF | gnuplot

set term post
set da s lp
set out 'be_O2_${dist}_O1.ps'

set xtics ("150" sqrt(150.), "100" sqrt(100.), "70" sqrt(70.), "50" sqrt(50.), "35" sqrt(35.), "25" sqrt(25.), "18" sqrt(18.), "13" sqrt(13.))
set multiplot
set origin 0,0
set size 1./3,1
set title "2x O"
plot [][:-860]'O1/E_O1.dat' u (sqrt(\$1)):(\$2*2*13.606) t "NC1", 'O1/E_O1.dat' u (sqrt(\$1)):(\$3*2*13.606) t "NC2", 'O1/E_O1.dat' u (sqrt(\$1)):(\$4*2*13.606) t "US1", 'O1/E_O1.dat' u (sqrt(\$1)):(\$5*2*13.606) t "US2", 'O1/E_O1.dat' u (sqrt(\$1)):(\$6*2*13.606) t "PAW1", 'O1/E_O1.dat' u (sqrt(\$1)):(\$7*2*13.606) t "PAW2"

set origin 1./3,0
set size 1./3,1
set title "O2_${dist}"
plot [][:-870]'O2_${dist}/E_O2_${dist}.dat' u (sqrt(\$1)):(\$2*1*13.606) t "NC1", 'O2_${dist}/E_O2_${dist}.dat' u (sqrt(\$1)):(\$3*1*13.606) t "NC2", 'O2_${dist}/E_O2_${dist}.dat' u (sqrt(\$1)):(\$4*1*13.606) t "US1", 'O2_${dist}/E_O2_${dist}.dat' u (sqrt(\$1)):(\$5*1*13.606) t "US2", 'O2_${dist}/E_O2_${dist}.dat' u (sqrt(\$1)):(\$6*1*13.606) t "PAW1", 'O2_${dist}/E_O2_${dist}.dat' u (sqrt(\$1)):(\$7*1*13.606) t "PAW2"

! paste O2_${dist}/E_O2_${dist}.dat O1/E_O1.dat > /tmp/tmp.be.dat
set origin 2./3,0
set size 1./3,1
set title "O2_${dist} - 2x O"
plot [][-9:-4]'/tmp/tmp.be.dat' u (sqrt(\$1)):((\$2-2*\$9)*13.606) t "NC1", '/tmp/tmp.be.dat' u (sqrt(\$1)):((\$3-2*\$10)*13.606) t "NC2", '/tmp/tmp.be.dat' u (sqrt(\$1)):((\$4-2*\$11)*13.606) t "US1", '/tmp/tmp.be.dat' u (sqrt(\$1)):((\$5-2*\$12)*13.606) t "US2", '/tmp/tmp.be.dat' u (sqrt(\$1)):((\$6-2*\$13)*13.606) t "PAW1", '/tmp/tmp.be.dat' u (sqrt(\$1)):((\$7-2*\$14)*13.606) t "PAW2"

set nomultiplot

! gv be_O2_${dist}_O1.ps  &
EOF
done
