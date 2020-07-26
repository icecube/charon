#!/bin/bash
if [ -z "$8" ]
then
	seed=1		
else
	seed=$8
fi		

for i in 0. 3. 6. 9. 12. 15. 25. 35. 45. 55. 65. 75. 85. 95. 105. 115. 125. 135. 145. 155.
do
  ./run.sh 	$1 $2 Sun decay secluded $3 $4 $5 $6 $i $7 $seed
done
