#!/bin/bash
# Author  :  Q.R. Liu
python set_up.py --channel $1 --mass $2 --Nevent $3
echo 'channel'      $1
echo 'DM mass'      $2 
echo 'Event Number' $3
echo 'seed'         $4
./main $1 $2 $3 $4 > $1_$2_$3_$4.output
