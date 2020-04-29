#!/bin/bash
./make.sh

for S in 400
do 

for SB in 50 100 150 200
do 

for ITERATION in $(seq 1)
do
 
SC=$(($S-$SB))

SEED=$(($SB * $SC * $ITERATION))

export SEED
export SB
export SC
export ITERATION

./atn $SEED $SB $SC $ITERATION

done
done 
done 
