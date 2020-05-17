#!/bin/bash
./make.sh

for WEB in 1 2 3 4 5
do

for S in 400
do 

for SB in 50
do 

SC=$(($S-$SB))

SEED=$(($SB * $SC *$WEB))

export SEED
export SB
export SC
export WEB

./atn $WEB $SEED $SB $SC

done 
done 
done
