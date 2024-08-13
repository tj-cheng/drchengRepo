#!/bin/bash

#IN=${1}
dd[1]=1.0
dd[2]=1.5
dd[3]=2.0
num_pelements=21  # 由於包含 0 和 1，所以有 11 個元素
for ((i=0; i<num_pelements; i++)); 
do
    value=$(echo "scale=1; $i * 0.05" | bc)  # 使用 bc 進行浮點運算
    prob+=($value)
done
#rm p10*
#rm p20*
#rm p30*
#rm p40*
#rm p50*
#rm p60*
#rm p70*
#rm p80*
#rm p90*
for (( i=1; i<=3; i=i+1 ))
do
   #j=$(( ${i}*8 ))
   for (( k=0; k<num_pelements; k++ ))
   do
      echo "{if(\$5==${dd[i]} && \$4==${prob[k]}) print \$0}" > gather.awk
      pp=$(echo "(${prob[k]}*100) / 1" | bc)
      dd=$(echo "(${dd[i]}*100) / 1" | bc)
      awk -f gather.awk ladderevol*.dat > p${pp}dd${dd}res.dat
      #awk -f unsus.awk p${pp}${var[i]}res.dat > p${pp}${var[i]}unsus.dat
      #awk -f locsus.awk p${pp}${var[i]}res.dat > p${pp}${var[i]}locsus.dat
   done
done
#python3 unsusdiagram.py
#python3 locsusdiagram.py
exit 0 
