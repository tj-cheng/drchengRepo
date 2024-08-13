#!/bin/bash

#sh mpicr.sh abcdefg.f90  
#1.mpif90 -o abcdefg.out abcdefg
#2.read the initial setting from initial.in
#  read defferent size lx and threads from (run.in)
#3.write the new information(lx) in data (read.in)
#4.mpirun -np threads ./abcdefg.out 
#5.wait previous task
# loop #3 #4 #5
#run.in 需寫入所要計算的尺寸。和要平行的數量
#initial.in 需寫入狀態 如samples mm0 等。

#initial setting
filename='read.in'
initial=$(cat initial.in)
#lx,ly,mm0,dd,init,samples,msteps,isteps
i=0
s=0
for condition in ${initial}
do
   con[i]=${condition}
   i=$(( ${i}+1 ))
done
#echo ${i}
lx=${con[0]}
ly=${con[1]}
mm0=${con[2]}
dd=${con[3]}
init=${con[4]}
samples=${con[5]}
mstps=${con[6]}
istps=${con[7]}

#echo "${lx} ${ly} ${mm0} ${dd}" > ${filename}
#echo "${init} ${samples} ${mstps} ${istps}" >> ${filename}  
#unset value
#value=$(<read.in)
#echo "${value}" 
#echo '----------'
runfile=$(cat run.in)
#echo ${runfile}
di=0
#ds=0
for variable in ${runfile}
do                      #save the variable in array var[]
   var[di]=${variable}
   di=$(( ${di}+1 ))
   #ds=$(( ${ds}+1 ))
done
ng=$(( ${di} / 2 ))   #number of lines in run.in
#nv=$(( ${di} / ${ng} ))   #number of columns in run.in
#echo $(( ${ds} / 3 ))
#echo ${nv}
#echo "lx=${var[3]}" #lx_1
#echo ${var[3]} #lx_2
#unset value 
#-------------------------------------

#exit 0
# >> append ; > new to write

#compile mpif90 -o cat.out cat.f90
IN=${1}

MPI=mpiifort

RUN=mpiexec

THREAD=-np
OUT=$(basename ${IN} .f90).out

${MPI} ${IN} -o ${OUT}
echo "${MPI} ${IN} -o ${OUT}"
##${RUN} ${THREAD} ${numnp} ${OUT}
#exit 0
#---------------------------------

#---run the code---
for (( i=0; i<=${ng}-1; i=i+1 ))
do
   lx=${var[i*2]}
   echo "read.in="
   echo "${lx} ${ly} ${mm0} ${dd}" > ${filename}
   echo "${init} ${samples} ${mstps} ${istps}" >> ${filename}  
   value=$(<read.in)
   echo "${value}" 
   unset value
   numth=${var[i*2+1]}
   time ${RUN} ./${OUT}
   echo "${RUN} ./${OUT}"
   echo " ${i} "
   PID=$!
   wait ${PID}             #wait the previous task
done
exit 0
#------------------

