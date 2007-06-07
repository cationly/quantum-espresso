#!/bin/bash
sh genAll.bash

for test in s p sp; do
#
# run the test
sh testAll-$test.bash
# collect results 
sh viewAll.bash $test
# save them with proper suffix
mv e.dat  e.dat-$test
mv s.dat  s.dat-$test
mv p.dat  p.dat-$test
mv all.eps all-$test.eps
#
done
