#!/bin/bash
PWROOT=$HOME/QE/PAW/espresso/

PW=$PWROOT/bin/pw.x
SCRATCH=/scratch/degironc/o1/
echo $PW, $SCRATCH

PSEUDOLIST="NC1.UPF NC2.UPF US1.UPF US2.UPF PAW1.PAW PAW2.PAW PAW3.PAW PAW4.PAW"
PSEUDOLIST="NC1.UPF NC2.UPF US1.UPF US2.UPF PAW1.PAW PAW3.PAW PAW4.PAW"
ECUTWFCLIST="100 84 70 59 50 42 35 29 25 21 18 15 13"
ECUTWFCLIST="100 70 50 35 25 18 13"


for nelec in 6.0 5.6 5.2 4.8 4.4 4.0; do

fileout="OAtom_${nelec}.dat"

echo -n "#Ecutwfc" > $fileout
for pseudo in $PSEUDOLIST; do
echo -n " $pseudo " >> $fileout
done
echo >> $fileout


for ecutwfc in $ECUTWFCLIST; do
echo -n $ecutwfc >> $fileout

for pseudo in $PSEUDOLIST; do
echo Now: $nelec $ecutwfc $pseudo

cat <<EOF | $PW > OAtom_${nelec}_${pseudo}_${ecutwfc}.out
O1
Oxygen
 &control
    calculation = 'scf'
    prefix='o1',
    pseudo_dir = './PSEUDO/'
    outdir='$SCRATCH/'
    disk_io='low'
 /
 &system
    ibrav=  2, celldm(1) =14.0, nat=  1, ntyp= 1, nelec=$nelec,
    ecutwfc= $ecutwfc
    ecutrho=400,
    occupations = 'smearing'
    smearing = 'gaussian'
    degauss = 0.01
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
O   0.0   0.0   0.0
K_POINTS automatic
 1 1 1 0 0 0
EOF

energy=$(awk '/\!/ {print $5}' OAtom_${nelec}_${pseudo}_${ecutwfc}.out)

echo -n " $energy " >> $fileout

done
echo >> $fileout
done

cat <<EOF | gnuplot
set term post
set out 'OAtom_${nelec}.ps'
set title "O1"
set xlabel "ecutwfc (Ryd)  (ecutrho=400)"
set ylabel "Total energy (Ryd)"
set da s lp
plot '$fileout' u 1:2 t "NC1", '$fileout' u 1:3 t "NC2", '$fileout' u 1:4 t "US1", '$fileout' u 1:5 t "US2", '$fileout' u 1:6 t "PAW1", '$fileout' u 1:7 t "PAW3", '$fileout' u 1:8 t "PAW4"
EOF

done

