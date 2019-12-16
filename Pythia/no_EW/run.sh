#!/bin/bash
source env.sh

echo 'channel'      	$1
echo 'DM mass'      	$2 
echo 'Number of bins'	$3
echo 'seed'         	$4
echo 'Event Number' 	$5
echo 'location' 	$6

python set_up.py --channel $1 --mass $2 --Nevent $5 --location $6

./main $1 $2 $3 $4 $6 > ./output/$1_$2_$3_$4_$5_$6.out & 
