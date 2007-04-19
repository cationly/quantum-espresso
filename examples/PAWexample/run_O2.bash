#!/bin/bash

PW=~/develop_PAW/bin/pw.x

for dist in 1.2 1.3 1.4 1.5 1.6 1.7; do

name=O2_$dist

fileout=E_$name.dat

PSEUDOLIST="NC1.UPF NC2.UPF US1.UPF US2.UPF PAW1.PAW PAW2.PAW PAW4.PAW" 
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

pwout=${name}_${ecutwfc}_${pseudo}.out
if [ $1 ]; then
if [ ! -f $pwout ]; then
cat <<EOF | $PW > $pwout
${name}
Oxygen-hydrogen
 &control
    calculation = 'scf'
    restart_mode='from_scratch',
    prefix='o2',
    tstress = .false.,
    tprnfor = .false.,
    pseudo_dir = './PSEUDO/'
    outdir='./tmp/'
    disk_io='low'
 /
 &system
    ibrav=  1, celldm(1) =10.0, nat=  2, ntyp= 1,
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
 H   1.000  H.pbe-rrkjus.UPF
ATOMIC_POSITIONS {bohr}
O  -$dist   0.0000   0.0000
O   $dist   0.0000   0.0000
K_POINTS automatic
 1 1 1 0 0 0
EOF
fi
fi

energy=$(awk '/\!/ {print $5}' $pwout)
echo -n " $energy " >> $fileout

done
echo >> $fileout
done


cat <<EOF | gnuplot
set term post
set out 'E_${name}.ps'
set title "${name}"
set xlabel "ecutwfc (Ryd)  (no dual)"
set ylabel "Total energy (Ryd)"
set da s lp
plot 'E_${name}.dat' u 1:2 t "NC1", 'E_${name}.dat' u 1:3 t "NC2", 'E_${name}.dat' u 1:4 t "US1", 'E_${name}.dat' u 1:5 t "US2", 'E_${name}.dat' u 1:6 t "PAW1", 'E_${name}.dat' u 1:7 t "PAW2"
#, 'E_${name}.dat' u 1:8 t "PAW3
EOF
done
