#!/bin/bash

LD1=~/develop_PAW/bin/ld1.x

RNC=1.40
RNChard=1.20
RUS=1.70
Rd=1.40

FUNC="PBE"
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

name=PAW-NChard-Qhard
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
   lnc2paw=.true.
   rcutnc2paw(1)=$RNChard,
   rcutnc2paw(2)=$RNChard,
   rcutnc2paw(3)=$RNChard,
   rcutnc2paw(4)=$RNChard,
   rcutnc2paw(5)=$RNChard,
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

name=PAW-NChard-Qsmooth
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
   lnc2paw=.true.
   rcutnc2paw(1)=$RNChard,
   rcutnc2paw(2)=$RNChard,
   rcutnc2paw(3)=$RNChard,
   rcutnc2paw(4)=$RNChard,
   rcutnc2paw(5)=$RNChard,
   which_paw_augfun ='QVAN'
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

