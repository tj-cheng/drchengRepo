#!/bin/bash

#IN=${1}
var[1]=128
var[2]=16
var[3]=32
var[4]=64
var[5]=128
prob[1]=1.0
prob[2]=0.9
prob[3]=0.8
prob[4]=0.7
prob[5]=0.6
prob[6]=0.5
prob[7]=0.4
prob[8]=0.3
prob[9]=0.2
prob[10]=0.1
prob[11]=0
rm p10*
rm p20*
rm p30*
rm p40*
rm p50*
rm p60*
rm p70*
rm p80*
rm p90*
for (( i=1; i<=1; i=i+1 ))
do
   j=$(( ${i}*8 ))
   for (( k=1; k<=11; k=k+1 ))
   do
      echo "{if(\$1==${var[i]} && \$8==${prob[k]}) print \$0}" > gather.awk
      pp=$(echo "(${prob[k]}*100) / 1" | bc)
      awk -f gather.awk res*.dat > p${pp}${var[i]}res.dat
      awk -f unsus.awk p${pp}${var[i]}res.dat > p${pp}${var[i]}unsus.dat
      awk -f locsus.awk p${pp}${var[i]}res.dat > p${pp}${var[i]}locsus.dat
   done
done
#python3 unsusdiagram.py
#python3 locsusdiagram.py
exit 0 
