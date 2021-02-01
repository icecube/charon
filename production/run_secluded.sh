#!/bin/bash
if [ -z "${8}" ]
then
	seed=1		
else
	seed=${8}
fi		

for i in 0.0 3.0 6.0 9.0 12.0 15.0 25.0 35.0 45.0 55.0 65.0 75.0 85.0 95.0 105.0 115.0 125.0 135.0 145.0 155.0
do
  ./run.sh $1 $2 $3 $4 secluded  $5 $6 0.0 - $i $7 $seed
done
