#!/bin/bash
PWROOT=$HOME/QE/PAW/espresso/

PW=$PWROOT/bin/pw.x
SCRATCH=/scratch/degironc/o1/
echo $PW, $SCRATCH

PSEUDOLIST="NC1.UPF NC2.UPF US1.UPF US2.UPF PAW1.PAW PAW2.PAW PAW3.PAW PAW4.PAW"
PSEUDOLIST="NC1.UPF NC2.UPF US1.UPF US2.UPF PAW1.PAW PAW3.PAW PAW4.PAW"
ECUTWFCLIST="100 84 70 59 50 42 35 29 25 21 18 15 13"
ECUTWFCLIST="100"

CALC[1]=NC1.UPF
CALC[2]=NC2.UPF
CALC[3]=US1.UPF
CALC[4]=US2.UPF
CALC[5]=PAW1.PAW
CALC[6]=PAW3.PAW
CALC[7]=PAW4.PAW
ECUTWFC[2]=100
ECUTWFC[3]=100
ECUTWFC[4]=100
ECUTWFC[5]=100
ECUTWFC[6]=100
ECUTWFC[7]=100

fileout="OAtom_scann.dat"

echo -n "#nelec" > $fileout
for pseudo in $PSEUDOLIST; do
echo -n " $pseudo " >> $fileout
done
echo >> $fileout

ecutwfc=70

for nelec in 6.0 5.6 5.2 4.8 4.4 4.0; do
echo -n $nelec >> $fileout

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
set out 'OAtom_scann.ps'
set title "O1"
set xlabel "ecutwfc (Ryd)  (ecutrho=400)"
set ylabel "Total energy (Ryd)"
set da s lp
plot '$fileout' u 1:2 t "NC1", '$fileout' u 1:3 t "NC2", '$fileout' u 1:4 t "US1", '$fileout' u 1:5 t "US2", '$fileout' u 1:6 t "PAW1", '$fileout' u 1:7 t "PAW3", '$fileout' u 1:8 t "PAW4"
EOF


