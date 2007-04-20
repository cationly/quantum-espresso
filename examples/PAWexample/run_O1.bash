#!/bin/bash

PW=~/Codes/espresso-develop_PAW/bin/pw.x

fileout="E_O1.dat"

PSEUDOLIST="NC1.UPF NC2.UPF US1.UPF US2.UPF PAW1.PAW PAW2.PAW PAW3.PAW PAW4.PAW"
ECUTWFCLIST="100 70 50 35 25 18 13"

echo -n "#Ecutwfc" > $fileout
for pseudo in $PSEUDOLIST; do
echo -n " $pseudo " >> $fileout
done
echo >> $fileout

for ecutwfc in $ECUTWFCLIST; do
echo -n $ecutwfc >> $fileout

for pseudo in $PSEUDOLIST; do
echo Now: $ecutwfc $pseudo

cat <<EOF | $PW > O1_${ecutwfc}_${pseudo}.out
O1
Oxygen
 &control
    calculation = 'scf'
    restart_mode='from_scratch',
    prefix='o1',
    tstress = .false.,
    tprnfor = .false.,
    pseudo_dir = './PSEUDO/'
    outdir='./tmp/'
    disk_io='low'
 /
 &system
    ibrav=  1, celldm(1) =10.0, nat=  1, ntyp= 1,
    ecutwfc= $ecutwfc
    !!ecutrho=...
    occupations = 'smearing'
    smearing = 'methfessel-paxton'
    degauss = 0.03
 /
 &electrons
    mixing_mode = 'plain'
    mixing_beta = 0.7
    mixing_ndim = 1
    conv_thr =  1.0d-6
 /
ATOMIC_SPECIES
 O   1.000  $pseudo
ATOMIC_POSITIONS {alat}
O   0.5000   0.0000   0.0000
K_POINTS automatic
 1 1 1 0 0 0
EOF

energy=$(awk '/\!/ {print $5}' O1_${ecutwfc}_${pseudo}.out)
echo -n " $energy " >> $fileout

done
echo >> $fileout
done

cat <<EOF | gnuplot
set term post
set out 'E_O1.ps'
set title "O1"
set xlabel "ecutwfc (Ryd)  (no dual)"
set ylabel "Total energy (Ryd)"
set da s lp
plot 'E_O1.dat' u 1:2 t "NC1", 'E_O1.dat' u 1:3 t "NC2", 'E_O1.dat' u 1:4 t "US1", 'E_O1.dat' u 1:5 t "US2", 'E_O1.dat' u 1:6 t "PAW1", 'E_O1.dat' u 1:7 t "PAW2", 'E_O1.dat' u 1:8 t "PAW3
EOF
