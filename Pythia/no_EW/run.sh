#!/bin/bash

echo 'channel'      	$1
echo 'DM mass'      	$2 
echo 'Number of bins'	$3
echo 'seed'         	$4
echo 'Event Number' 	$5
echo 'location' 	$6

python set_up.py --channel $1 --mass $2 --Nevent $5

./main $1 $2 $3 $4 $6
