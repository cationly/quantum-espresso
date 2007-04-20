#!/bin/bash

NCONF=11
CALCLIST="NC1 NC2 US1 US2 PAW1 PAW2 PAW3 PAW4"

first="NC1"
last="PAW4"

echo "# occ AE "$CALCLIST > s.dat
echo "# occ AE "$CALCLIST > p.dat
echo "# occ AE "$CALCLIST > e.dat

for ((i=0;((i<NCONF));i++)); do
occ=$(dc -e "2k 2 $((NCONF-1))/ $i* _1 * 2+ p")
#occ=$(dc -e "2k 1 $((NCONF-1))/ $i* _1 * 2+ p")
echo $occ

for calc in $CALCLIST; do
awk -v calc=$calc -v occ=$occ -v first=$first -v last=$last '
/     1 0     2S   1\(/ {sae=$6; sps=$7}
/     2 1     2P   1\(/ {pae=$6; pps=$7}
/dEtot_ae/ {eae=$3}
/dEtot_ps/ {eps=$3}
END {
if (calc==first) printf "%f ", occ >>"s.dat";
if (calc==first) printf "%f ", occ >>"p.dat";
if (calc==first) printf "%f ", occ >>"e.dat";
if (calc==first) printf "%f ", sae >>"s.dat";
if (calc==first) printf "%f ", pae >>"p.dat";
if (calc==first) printf "%f ", eae >>"e.dat";
printf "%f ",sps-sae >>"s.dat";
printf "%f ",pps-pae >>"p.dat";
printf "%f ",eps-eae >>"e.dat";
if (calc==last) printf "\n" >>"s.dat";
if (calc==last) printf "\n" >>"p.dat";
if (calc==last) printf "\n" >>"e.dat";
}' $calc/$calc.tst_$occ.out
done
done

cat <<EOF | gnuplot
set term post eps dl 2
set out 'all.eps'
set da s l
set xlabel "2s occupation"
set multiplot
set origin 0,0
set size 0.5,1
set ylabel "E(2s): PS - AE (Ry)"
set key bottom right
plot 's.dat' u 1:3 t "NC1", 's.dat' u 1:4 t "NC2", 's.dat' u 1:5 t "US1", 's.dat' u 1:6 t "US2", 's.dat' u 1:7 t "PAW1" w p, 's.dat' u 1:8 t "PAW2" w p, 's.dat' u 1:9 t "PAW3" w p, 's.dat' u 1:10 t "PAW4" w p
set nokey
set origin 0.55,0.5
set size 0.45,0.5
set ylabel "E(2p): PS - AE (Ry)"
plot 'p.dat' u 1:3 t "NC1", 'p.dat' u 1:4 t "NC2", 'p.dat' u 1:5 t "US1", 'p.dat' u 1:6 t "US2", 'p.dat' u 1:7 t "PAW1" w p, 'p.dat' u 1:8 t "PAW2" w p, 'p.dat' u 1:9 t "PAW3" w p, 'p.dat' u 1:10 t "PAW4" w p
set origin 0.55,0.0
set size 0.45,0.5
set ylabel "deltaE: PS - AE (Ry)"
plot 'e.dat' u 1:3 t "NC1", 'e.dat' u 1:4 t "NC2", 'e.dat' u 1:5 t "US1", 'e.dat' u 1:6 t "US2", 'e.dat' u 1:7 t "PAW1" w p, 'e.dat' u 1:8 t "PAW2" w p, 'e.dat' u 1:9 t "PAW3" w p, 'e.dat' u 1:10 t "PAW4" w p
set nomultiplot
! ggv -antialias -scale 3 all.eps &
EOF

exit

#############boh da guardare################
cat <<EOF | gnuplot
set term post eps dl 2
set out 'all2.eps'
set da s l
set xlabel "2s occupation"
set multiplot
set origin 0,0
set size 0.5,1
set ylabel "E(2s): PAW - reference (Ry)"
set key bottom right
plot 's.dat' u 1:7 t "PAW-AE-QAE", 's.dat' u 1:(\$8-\$4) t "PAW-NChard-Qhard" w p, 's.dat' u 1:(\$9-\$4) t "PAW-NChard-Qsmooth"
set nokey
set origin 0.55,0.5
set size 0.45,0.5
set ylabel "E(2p):  PAW - reference (Ry)"
plot 'p.dat' u 1:7 t "PAW-AE-QAE", 'p.dat' u 1:(\$8-\$4) t "PAW-NChard-Qhard" w p, 'p.dat' u 1:(\$9-\$4) t "PAW-NChard-Qsmooth"
set origin 0.55,0.0
set size 0.45,0.5
set ylabel "deltaE:  PAW - reference (Ry)"
plot 'e.dat' u 1:7 t "PAW-AE-QAE", 'e.dat' u 1:(\$8-\$4) t "PAW-NChard-Qhard" w p, 'e.dat' u 1:(\$9-\$4) t "PAW-NChard-Qsmooth"
set nomultiplot
! ggv -antialias -scale 3 all2.eps &
EOF
