#!/bin/bash

echo 'channel'      		$1
echo 'DM Mass'      		$2 
echo 'Location' 			$3
echo 'Process'  			$4
echo 'Type'  	        	$5
echo 'Event Number' 		$6
echo 'Number of Bins'		$7
echo 'Lower Energy Bound'   $8
echo 'Binning Scale'    	$9


if [ $5 == 'secluded' ]
then	
  echo 'Production Spot Density'  	${10}
  echo 'Mediator Mass'              ${11}
  
  python set_up.py --channel $1 --mass $2 --bins $7  --Nevent $6  --location $3 --process $4 --type $5 --mass_phi ${11}
  if [ -z "${12}" ]
  then 
    ./mainv2 $1 $2 $7 $3 $4 ${10} 1 $9 $8 ${11} $5 > ./output/$1_$2_${11}_$6_${10}_$4_$5_$3.out & 
  else
    echo 'Seed' ${12}	  
    ./mainv2 $1 $2 $7 $3 $4 ${10} ${12} $9 $8 ${11} $5 > ./output/$1_$2_${11}_$6_${10}_$4_$5_$3.out & 
  fi
    
elif [ $5 != 'secluded' ]
  then
  python set_up.py --channel $1 --mass $2 --bins $7  --Nevent $6  --location $3 --process $4 --type $5
  if [ $1 != 'nuenue' ] && [ $1 != 'numunumu' ] && [ $1 != 'nutaunutau' ]
  then	
    if [ $3 == 'Sun' ]
    then
      if [ -z "${10}" ]
      then 
      ./mainv2 $1 $2 $7 $3 $4 148.9 1 $9 $8 $2 $5  > ./output/$1_$2_$6_${10}_$4_$5_$3.out & 
      else
      echo 'seed' ${10}	  
      ./mainv2 $1 $2 $7 $3 $4 148.9 ${10} $9 $8 $2 $5  > ./output/$1_$2_$6_${10}_$4_$5_$3.out & 
      fi 
    elif [ $3 == 'Earth' ]
    then
      if [ -z "${10}" ]
      then 
      ./mainv2 $1 $2 $7 $3 $4 13.08849999999999802 1 $9 $8 $2 $5  > ./output/$1_$2_$6_${10}_$4_$5_$3.out & 
      else
      echo 'seed' ${10}	  
      ./mainv2 $1 $2 $7 $3 $4 13.08849999999999802 ${10} $9 $8 $2 $5 > ./output/$1_$2_$6_${10}_$4_$5_$3.out & 
      fi 
    elif [ $3 == 'Halo' ]
    then
      if [ -z "${10}" ]
      then 
      ./mainv2 $1 $2 $7 $3 $4 0. 1 $9 $8 $2 $5  > ./output/$1_$2_$6_${10}_$4_$5_$3.out & 
      else
      echo 'seed' ${10}	  
      ./mainv2 $1 $2 $7 $3 $4 0. ${10} $9 $8 $2 $5 > ./output/$1_$2_$6_${10}_$4_$5_$3.out & 
      fi 
    fi
  fi  
fi

echo "Start running ..."
