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
for (( i=3; i<=3; i=i+1 ))
do
   #j=$(( ${i}*8 ))
   #for (( k=0; k<num_pelements; k++ ))
   for (( k=0; k<num_pelements; k++ ))
   do
      pp=$(echo "(${prob[k]}*100) / 1" | bc)
      dd=$(echo "(${dd[i]}*100) / 1" | bc)
      echo "p${pp}dd${dd}res.dat"
      awk -f max.awk p${pp}dd${dd}res.dat > max.in 
      ##################
      initial=$(cat max.in)
      j=0
      for condition in ${initial}
      do
         con[j]=${condition}
         j=$(( ${j}+1 ))
      done
      mstep=${con[0]}
      echo "${mstep}"
      ##################
      echo "{print int(\$3/$((${mstep}))*360),\$2,\$4,\$5,\$1}" > draw.awk
      echo "${i}, p${pp}dd${dd}res0.dat"
      awk -f draw.awk p${pp}dd${dd}res.dat > p${pp}dd${dd}res0.dat
      awk -f av.awk p${pp}dd${dd}res0.dat > p${pp}dd${dd}evo.dat
   done
done
exit 0
