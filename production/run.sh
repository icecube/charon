#!/bin/bash

echo 'channel'      	$1
echo 'DM mass'      	$2 
echo 'Number of bins'	$3
echo 'Event Number' 	$4
echo 'location' 	$5
echo 'process'  	$6
echo 'type'  	        $7

if [ $7 == 'secluded' ]
then	
    python set_up.py --channel $1 --mass $9 --bins $3  --Nevent $4  --location $5 --process $6 --type $7
  echo 'density'  	$8
  echo 'mediator mass'  $9
  if [ -z "${10}" ]
  then 
  ./main $1 $9 $3 $5 $6 $8 1   > ./output/$1_$9_$3_$4_$8_$7.out & 
  else
  echo 'seed' ${10}	  
  ./main $1 $9 $3 $5 $6 $8 ${10} > ./output/$1_$9_$3_$4_$8_$7.out & 
  fi
elif [ $7 != 'secluded' ]
  then
  python set_up.py --channel $1 --mass $2 --bins $3  --Nevent $4  --location $5 --process $6 --type $
  if [ $1 != 'nuenue' ] && [ $1 != 'numunumu' ] && [ $1 != 'nutaunutau' ]
  then	
    if [ $5 == 'Sun' ]
    then
      if [ -z "$8" ]
      then 
      ./main $1 $2 $3 $5 $6 148.9 1> ./output/$1_$2_$3_$4_$6.out & 
      else
      echo 'seed' $8	  
      ./main $1 $2 $3 $5 $6 148.9 $8> ./output/$1_$2_$3_$4_$6.out &
      fi 
    elif [ $5 == 'Earth' ]
    then
      if [ -z "$8" ]
      then 
      ./main $1 $2 $3 $5 $6 13.08849999999999802 1> ./output/$1_$2_$3_$4_$6.out & 
      else
      echo 'seed' $8	  
      ./main $1 $2 $3 $5 $6 13.08849999999999802 $8 > ./output/$1_$2_$3_$4_$6.out & 
      fi 
    elif [ $5 == 'Galactic' ]
    then
      if [ -z "$8" ]
      then 
      ./main $1 $2 $3 $5 $6 0. 1> ./output/$1_$2_$3_$4_$6.out & 
      else
      echo 'seed' $8	  
      ./main $1 $2 $3 $5 $6 0. $8> ./output/$1_$2_$3_$4_$6.out &
      fi 
    fi
  fi  
fi

echo "Start running ..."

