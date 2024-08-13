#!/bin/bash


IN=${1}
awk 'BEGIN{max=0} {if( $3+0 > max+0) max=$3} END {print "MAX=", max}' ${IN} > maxrgstep.txt

filename=$(cat maxrgstep.txt)
i=0
for step in ${filename}
do
   maxstep[i]=${step}
done

#echo ${maxstep[0]}
