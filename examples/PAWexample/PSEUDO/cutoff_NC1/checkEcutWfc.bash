#!/bin/bash
PW=~/develop_PAW/bin/pw.x
pseudo=NC1.UPF

for ecutwfc in 10 15 20 25 30 35 40 45 50 60 70 80 90 100 150 200; do
if [ ! -f $pseudo_$ecutwfc.out ]; then
cat <<EOF | $PW > $pseudo_$ecutwfc.out
O1
Oxygen
 &control
    calculation = 'scf'
    restart_mode='from_scratch',
    prefix='o1',
    tstress = .false.,
    tprnfor = .false.,
    pseudo_dir = '../PSEUDO/'
    outdir='../tmp/'
    disk_io='low'
 /
 &system
    ibrav=  1, celldm(1) =8.0, nat=  1, ntyp= 1,
    ecutwfc= $ecutwfc
    !! ecutrho=150
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
O   0.0000   0.0000   0.0000
K_POINTS automatic
 1 1 1 0 0 0
EOF
fi
done
