#!/bin/bash
#
# NB: edit file config.bash in order to set the XC functional used and 
#     number of configurations
#
. config.bash
echo "running script $0 (DFT = $FUNC) on $NCONF configurations"

LD1=../../../bin/ld1.x

CALC[1]=NC
PROG[1]=$LD1
SUFF[1]=UPF

CALC[2]=NChard
PROG[2]=$LD1
SUFF[2]=UPF

CALC[3]=US
PROG[3]=$LD1
SUFF[3]=UPF

CALC[4]=UShard
PROG[4]=$LD1
SUFF[4]=UPF

CALC[5]=PAW-AE-QAE
PROG[5]=$LD1
SUFF[5]=PAW

CALC[6]=PAW-AE-BESS
PROG[6]=$LD1
SUFF[6]=PAW

CALC[7]=PAW-AE-GAUSS
PROG[7]=$LD1
SUFF[7]=PAW

for ((i=0;((i<NCONF));i++)); do
occ=$(dc -e "2k 2 $((NCONF-1))/ $i* _1 * 2+ p")
occ2=$(dc -e "2k 4 2 $((NCONF-1))/ $i* -  p")
echo "2S occupation 2.00 ; 2P occupation $occ2"
cat <<EOF > testconfigs
2
2S  1  0  2.00  0.00  1.40  1.60  1
2P  2  1  4.00  0.00  1.40  1.60  1
2
2S  1  0  2.00  0.00  1.40  1.60  1
2P  2  1  $occ2  0.00  1.40  1.60  1
EOF

for ((c=1;((c<=7));c++)); do

cat <<EOF > ${CALC[$c]}/${CALC[$c]}.tst.in
 &input
        title='O',
        zed=8.0,
        rel=0,
        beta=0.5,
        iswitch=2,
        dft='$FUNC',
        prefix='./${CALC[$c]}/${CALC[$c]}'
        nld=0
 /
4
1S  1  0  2.0  1
2S  2  0  2.0  1
2P  2  1  4.0  1
3D  3  2 -1.0  1
 &test
  file_pseudo='./${CALC[$c]}/${CALC[$c]}.${SUFF[$c]}'
  nconf=2,
 / 
EOF
cat testconfigs >> ${CALC[$c]}/${CALC[$c]}.tst.in
${PROG[$c]} < ${CALC[$c]}/${CALC[$c]}.tst.in > ${CALC[$c]}/${CALC[$c]}.tst_$occ.out

done

done

rm -f testconfigs
