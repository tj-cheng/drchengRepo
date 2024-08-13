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
#lx,ly,beta,dd,samples,msteps,isteps
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
beta=${con[2]}
dd=${con[3]}
rgprob=${con[4]}
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

${MPI} ${IN} -o ${OUT} -O3
echo "${MPI} ${IN} -o ${OUT} -O3"
##${RUN} ${THREAD} ${numnp} ${OUT}
#exit 0
#---------------------------------

#---run the code---
#for (( i=0; i<=${ng}-1; i=i+1 ))
#do
#   beta=${var[i*2]}
#   echo "read.in="
#   #echo "${lx} ${ly} ${beta} ${dd}" > ${filename}
#   echo "8 ${ly} ${beta} ${dd} ${rgprob}" > ${filename}
#   echo "${samples} ${mstps} ${istps}" >> ${filename}  
#   value=$(<read.in)
#   echo "${value}" 
#   unset value
#   numth=${var[i*2+1]}
#   time ${RUN} ./${OUT}
#   echo "${RUN} ./${OUT}"
#   PID=$!
#   wait ${PID}             #wait the previous task
#done
#for (( i=0; i<=${ng}-1; i=i+1 ))
#do
#   beta=${var[i*2]}
#   echo "read.in="
#   echo "16 ${ly} ${beta} ${dd} ${rgprob}" > ${filename}
#   echo "${samples} ${mstps} ${istps}" >> ${filename}  
#   value=$(<read.in)
#   echo "${value}" 
#   unset value
#   numth=${var[i*2+1]}
#   time ${RUN} ./${OUT}
#   echo "${RUN} ./${OUT}"
#   echo " ${i} "
#   PID=$!
#   wait ${PID}             #wait the previous task
#done
#for (( i=0; i<=${ng}-1; i=i+1 ))
#do
#   beta=${var[i*2]}
#   echo "read.in="
#   echo "32 ${ly} ${beta} ${dd} ${rgprob}" > ${filename}
#   echo "${samples} ${mstps} ${istps}" >> ${filename}  
#   value=$(<read.in)
#   echo "${value}" 
#   unset value
#   numth=${var[i*2+1]}
#   time ${RUN} ./${OUT}
#   echo "${RUN} ./${OUT}"
#   echo " ${i} "
#   PID=$!
#   wait ${PID}             #wait the previous task
#done
#for (( i=0; i<=${ng}-1; i=i+1 ))
#do
#   beta=${var[i*2]}
#   echo "read.in="
#   echo "64 ${ly} ${beta} ${dd} ${rgprob}" > ${filename}
#   echo "${samples} ${mstps} ${istps}" >> ${filename}  
#   value=$(<read.in)
#   echo "${value}" 
#   unset value
#   numth=${var[i*2+1]}
#   time ${RUN} ./${OUT}
#   echo "${RUN} ./${OUT}"
#   echo " ${i} "
#   PID=$!
#   wait ${PID}             #wait the previous task
#done
#for (( i=0; i<=${ng}-1; i=i+1 ))
#do
#   beta=${var[i*2]}
#   echo "read.in="
#   echo "128 ${ly} ${beta} ${dd} 1.0" > ${filename}
#   echo "${samples} ${mstps} ${istps}" >> ${filename}  
#   value=$(<read.in)
#   echo "${value}" 
#   unset value
#   numth=${var[i*2+1]}
#   time ${RUN} ./${OUT}
#   echo "${RUN} ./${OUT}"
#   echo " ${i} "
#   PID=$!
#   wait ${PID}             #wait the previous task
#done
#for (( i=0; i<=${ng}-1; i=i+1 ))
#do
#   beta=${var[i*2]}
#   echo "read.in="
#   echo "128 ${ly} ${beta} ${dd} 0.9" > ${filename}
#   echo "${samples} ${mstps} ${istps}" >> ${filename}  
#   value=$(<read.in)
#   echo "${value}" 
#   unset value
#   numth=${var[i*2+1]}
#   time ${RUN} ./${OUT}
#   echo "${RUN} ./${OUT}"
#   echo " ${i} "
#   PID=$!
#   wait ${PID}             #wait the previous task
#done
#for (( i=0; i<=${ng}-1; i=i+1 ))
#do
#   beta=${var[i*2]}
#   echo "read.in="
#   echo "128 ${ly} ${beta} ${dd} 0.8" > ${filename}
#   echo "${samples} ${mstps} ${istps}" >> ${filename}  
#   value=$(<read.in)
#   echo "${value}" 
#   unset value
#   numth=${var[i*2+1]}
#   time ${RUN} ./${OUT}
#   echo "${RUN} ./${OUT}"
#   echo " ${i} "
#   PID=$!
#   wait ${PID}             #wait the previous task
#done
#for (( i=0; i<=${ng}-1; i=i+1 ))
#do
#   beta=${var[i*2]}
#   echo "read.in="
#   echo "128 ${ly} ${beta} ${dd} 0.7" > ${filename}
#   echo "${samples} ${mstps} ${istps}" >> ${filename}  
#   value=$(<read.in)
#   echo "${value}" 
#   unset value
#   numth=${var[i*2+1]}
#   time ${RUN} ./${OUT}
#   echo "${RUN} ./${OUT}"
#   echo " ${i} "
#   PID=$!
#   wait ${PID}             #wait the previous task
#done
#for (( i=0; i<=${ng}-1; i=i+1 ))
#do
#   beta=${var[i*2]}
#   echo "read.in="
#   echo "128 ${ly} ${beta} ${dd} 0.6" > ${filename}
#   echo "${samples} ${mstps} ${istps}" >> ${filename}  
#   value=$(<read.in)
#   echo "${value}" 
#   unset value
#   numth=${var[i*2+1]}
#   time ${RUN} ./${OUT}
#   echo "${RUN} ./${OUT}"
#   echo " ${i} "
#   PID=$!
#   wait ${PID}             #wait the previous task
#done
#for (( i=0; i<=${ng}-1; i=i+1 ))
#do
#   beta=${var[i*2]}
#   echo "read.in="
#   echo "128 ${ly} ${beta} ${dd} 0.5" > ${filename}
#   echo "${samples} ${mstps} ${istps}" >> ${filename}  
#   value=$(<read.in)
#   echo "${value}" 
#   unset value
#   numth=${var[i*2+1]}
#   time ${RUN} ./${OUT}
#   echo "${RUN} ./${OUT}"
#   echo " ${i} "
#   PID=$!
#   wait ${PID}             #wait the previous task
#done
#for (( i=0; i<=${ng}-1; i=i+1 ))
#do
#   beta=${var[i*2]}
#   echo "read.in="
#   echo "128 ${ly} ${beta} ${dd} 0.4" > ${filename}
#   echo "${samples} ${mstps} ${istps}" >> ${filename}  
#   value=$(<read.in)
#   echo "${value}" 
#   unset value
#   numth=${var[i*2+1]}
#   time ${RUN} ./${OUT}
#   echo "${RUN} ./${OUT}"
#   echo " ${i} "
#   PID=$!
#   wait ${PID}             #wait the previous task
#done
#for (( i=0; i<=${ng}-1; i=i+1 ))
#do
#   beta=${var[i*2]}
#   echo "read.in="
#   echo "128 ${ly} ${beta} ${dd} 0.3" > ${filename}
#   echo "${samples} ${mstps} ${istps}" >> ${filename}  
#   value=$(<read.in)
#   echo "${value}" 
#   unset value
#   numth=${var[i*2+1]}
#   time ${RUN} ./${OUT}
#   echo "${RUN} ./${OUT}"
#   echo " ${i} "
#   PID=$!
#   wait ${PID}             #wait the previous task
#done
#for (( i=0; i<=${ng}-1; i=i+1 ))
#do
#   beta=${var[i*2]}
#   echo "read.in="
#   echo "128 ${ly} ${beta} ${dd} 0.2" > ${filename}
#   echo "${samples} ${mstps} ${istps}" >> ${filename}  
#   value=$(<read.in)
#   echo "${value}" 
#   unset value
#   numth=${var[i*2+1]}
#   time ${RUN} ./${OUT}
#   echo "${RUN} ./${OUT}"
#   echo " ${i} "
#   PID=$!
#   wait ${PID}             #wait the previous task
#done
#for (( i=0; i<=${ng}-1; i=i+1 ))
#do
#   beta=${var[i*2]}
#   echo "read.in="
#   echo "128 ${ly} ${beta} ${dd} 0.1" > ${filename}
#   echo "${samples} ${mstps} ${istps}" >> ${filename}  
#   value=$(<read.in)
#   echo "${value}" 
#   unset value
#   numth=${var[i*2+1]}
#   time ${RUN} ./${OUT}
#   echo "${RUN} ./${OUT}"
#   echo " ${i} "
#   PID=$!
#   wait ${PID}             #wait the previous task
#done
for (( i=0; i<=${ng}-1; i=i+1 ))
do
   beta=${var[i*2]}
   echo "read.in="
   echo "128 ${ly} ${beta} ${dd} 0.0" > ${filename}
   echo "${samples} ${mstps} ${istps}" >> ${filename}  
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

