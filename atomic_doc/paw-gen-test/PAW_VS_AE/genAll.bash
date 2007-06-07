#!/bin/bash 
#
# NB: edit file config.bash in order to set the XC functional used
#
. config.bash
#
echo "running the Pseudopotential generation script for DFT =" $FUNC

LD1=../../../bin/ld1.x

RNC=1.40
RNChard=1.20
RUS=1.70
Rd=1.40
Rmatch_augfun=1.2

NLCC="true"

### NC ###

name=NC
mkdir -p $name
cat <<EOF > $name/$name.gen.in
 &input
        title='O',
        zed=8.0,
        rel=0,
        beta=0.5,
        iswitch=3,
        dft='$FUNC',
        prefix='./$name/$name'
        nld=0
 /
4
1S  1  0  2.0  1
2S  2  0  2.0  1
2P  2  1  4.0  1
3D  3  2 -1.0  1
 &inputp
   pseudotype=3,
   nlcc=.$NLCC., rcore=0.50
   lloc=2,
   file_pseudopw='./$name/$name.UPF'
   zval=6.d0,
   lpaw=.false.
 /
5
2S  2  0  2.00  0.00  $RNC  $RNC
2S  2  0  0.00  0.05  $RNC  $RNC
2P  2  1  4.00  0.00  $RNC  $RNC
2P  2  1  0.00  0.05  $RNC  $RNC
3D  3  2 -2.00  0.15  $Rd  $Rd
 &test
  nconf=1,
 / 
2
2S  1  0  2.00  0.00  $RNC  $RUS  1
2P  2  1  4.00  0.00  $RNC  $RUS  1
EOF
$LD1 < $name/$name.gen.in > $name/$name.gen.out

### NC hard ###

name=NChard
mkdir -p $name
cat <<EOF > $name/$name.gen.in
 &input
        title='O',
        zed=8.0,
        rel=0,
        beta=0.5,
        iswitch=3,
        dft='$FUNC',
        prefix='./$name/$name'
        nld=0
 /
4
1S  1  0  2.0  1
2S  2  0  2.0  1
2P  2  1  4.0  1
3D  3  2 -1.0  1
 &inputp
   pseudotype=3,
   nlcc=.$NLCC., rcore=0.50
   lloc=2,
   file_pseudopw='./$name/$name.UPF'
   zval=6.d0,
   lpaw=.false.
 /
5
2S  2  0  2.00  0.00  $RNChard  $RNChard
2S  2  0  0.00  0.05  $RNChard  $RNChard
2P  2  1  4.00  0.00  $RNChard  $RNChard
2P  2  1  0.00  0.05  $RNChard  $RNChard
3D  3  2 -2.00  0.15  $Rd  $Rd
 &test
  nconf=1,
 / 
2
2S  1  0  2.00  0.00  $RNC  $RUS  1
2P  2  1  4.00  0.00  $RNC  $RUS  1
EOF
$LD1 < $name/$name.gen.in > $name/$name.gen.out

### US ###

name=US
mkdir -p $name
cat <<EOF > $name/$name.gen.in
 &input
        title='O',
        zed=8.0,
        rel=0,
        beta=0.5,
        iswitch=3,
        dft='$FUNC',
        prefix='./$name/$name'
        nld=0
 /
4
1S  1  0  2.0  1
2S  2  0  2.0  1
2P  2  1  4.0  1
3D  3  2 -1.0  1
 &inputp
   pseudotype=3,
   nlcc=.$NLCC., rcore=0.50
   lloc=2,
   file_pseudopw='./$name/$name.UPF'
   zval=6.d0,
   lpaw=.false.
   file_qvan = './$name/$name.qvan'
 /
5
2S  2  0  2.00  0.00  $RNC  $RUS
2S  2  0  0.00  0.05  $RNC  $RUS
2P  2  1  4.00  0.00  $RNC  $RUS
2P  2  1  0.00  0.05  $RNC  $RUS
3D  3  2 -2.00  0.15  $Rd  $Rd
 &test
  nconf=1,
 / 
2
2S  1  0  2.00  0.00  $RNC  $RUS  1
2P  2  1  4.00  0.00  $RNC  $RUS  1
EOF
$LD1 < $name/$name.gen.in > $name/$name.gen.out

### US with hard augmentation charges ###

name=UShard
mkdir -p $name
cat <<EOF > $name/$name.gen.in
 &input
        title='O',
        zed=8.0,
        rel=0,
        beta=0.5,
        iswitch=3,
        dft='$FUNC',
        prefix='./$name/$name'
        nld=0
 /
4
1S  1  0  2.0  1
2S  2  0  2.0  1
2P  2  1  4.0  1
3D  3  2 -1.0  1
 &inputp
   pseudotype=3,
   nlcc=.$NLCC., rcore=0.50
   lloc=2,
   file_pseudopw='./$name/$name.UPF'
   zval=6.d0,
   lpaw=.false.
   file_qvan = './$name/$name.qvan'
 /
5
2S  2  0  2.00  0.00  $RNChard  $RUS
2S  2  0  0.00  0.05  $RNChard  $RUS
2P  2  1  4.00  0.00  $RNChard  $RUS
2P  2  1  0.00  0.05  $RNChard  $RUS
3D  3  2 -2.00  0.15  $Rd  $Rd
 &test
  nconf=1,
 / 
2
2S  1  0  2.00  0.00  $RNC  $RUS  1
2P  2  1  4.00  0.00  $RNC  $RUS  1
EOF
$LD1 < $name/$name.gen.in > $name/$name.gen.out

### PAW from AE with AE augm func ###

name=PAW-AE-QAE
mkdir -p $name
cat <<EOF > $name/$name.gen.in
 &input
        title='O',
        zed=8.0,
        rel=0,
        beta=0.5,
        iswitch=3,
        dft='$FUNC',
        prefix='./$name/$name'
        nld=0
 /
4
1S  1  0  2.0  1
2S  2  0  2.0  1
2P  2  1  4.0  1
3D  3  2 -1.0  1
 &inputp
   pseudotype=3,
   nlcc=.$NLCC., rcore=0.50
   lloc=2,
   file_pseudopw='./$name/$name.PAW'
   zval=6.d0,
   lpaw=.true.
   file_qvan = './$name/$name.qvan'
 /
5
2S  2  0  2.00  0.00  $RNC  $RUS
2S  2  0  0.00  0.05  $RNC  $RUS
2P  2  1  4.00  0.00  $RNC  $RUS
2P  2  1  0.00  0.05  $RNC  $RUS
3D  3  2 -2.00  0.15  $Rd  $Rd
 &test
  nconf=1,
 / 
2
2S  1  0  2.00  0.00  $RNC  $RUS  1
2P  2  1  4.00  0.00  $RNC  $RUS  1
EOF
$LD1 < $name/$name.gen.in > $name/$name.gen.out

### PAW from NC hard with NChard augmentation functions ###

name=PAW-AE-BESS
mkdir -p $name
cat <<EOF > $name/$name.gen.in
 &input
        title='O',
        zed=8.0,
        rel=0,
        beta=0.5,
        iswitch=3,
        dft='$FUNC',
        prefix='./$name/$name'
        nld=0
 /
4
1S  1  0  2.0  1
2S  2  0  2.0  1
2P  2  1  4.0  1
3D  3  2 -1.0  1
 &inputp
   pseudotype=3,
   nlcc=.$NLCC., rcore=0.50
   lloc=2,
   file_pseudopw='./$name/$name.PAW'
   zval=6.d0,
   lpaw=.true.
   which_paw_augfun ='BESSEL'
   paw_rmatch_augfun=$Rmatch_augfun
   file_qvan = './$name/$name.qvan'
 /
5
2S  2  0  2.00  0.00  $RNC  $RUS
2S  2  0  0.00  0.05  $RNC  $RUS
2P  2  1  4.00  0.00  $RNC  $RUS
2P  2  1  0.00  0.05  $RNC  $RUS
3D  3  2 -2.00  0.15  $Rd  $Rd
 &test
  nconf=1,
 / 
2
2S  1  0  2.00  0.00  $RNC  $RUS  1
2P  2  1  4.00  0.00  $RNC  $RUS  1
EOF
$LD1 < $name/$name.gen.in > $name/$name.gen.out

### PAW from NC hard with smoother NC augmentation functions (same radius as US calculation) ###

name=PAW-AE-GAUSS
mkdir -p $name
cat <<EOF > $name/$name.gen.in
 &input
        title='O',
        zed=8.0,
        rel=0,
        beta=0.5,
        iswitch=3,
        dft='$FUNC',
        prefix='./$name/$name'
        nld=0
 /
4
1S  1  0  2.0  1
2S  2  0  2.0  1
2P  2  1  4.0  1
3D  3  2 -1.0  1
 &inputp
   pseudotype=3,
   nlcc=.$NLCC., rcore=0.50
   lloc=2,
   file_pseudopw='./$name/$name.PAW'
   zval=6.d0,
   lpaw=.true.
   which_paw_augfun ='GAUSS'
   paw_rmatch_augfun=$Rmatch_augfun
   file_qvan = './$name/$name.qvan'
 /
5
2S  2  0  2.00  0.00  $RNC  $RUS
2S  2  0  0.00  0.05  $RNC  $RUS
2P  2  1  4.00  0.00  $RNC  $RUS
2P  2  1  0.00  0.05  $RNC  $RUS
3D  3  2 -2.00  0.15  $Rd  $Rd
 &test
  nconf=1,
 / 
2
2S  1  0  2.00  0.00  $RNC  $RUS  1
2P  2  1  4.00  0.00  $RNC  $RUS  1
EOF
$LD1 < $name/$name.gen.in > $name/$name.gen.out
