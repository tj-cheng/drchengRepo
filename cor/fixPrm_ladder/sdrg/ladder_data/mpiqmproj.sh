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
trialcase=${con[2]}
init=${con[3]}
samples=${con[4]}


runfile2=$(cat lseed.in)
di=0
for variable2 in ${runfile2}
do 
   var2[di]=${variable2}
   di=$(( ${di}+1 ))
done
nsize=${di}
#-------------------------------------

#exit 0
# >> append ; > new to write

#compile mpif90 -o cat.out cat.f90
IN=${1}

#MPI=mpiifort

MPI=mpif90

#RUN=mpiexec
RUN=mpirun

THREAD=-np

OUT=$(basename ${IN} .f90).out

${MPI} ${IN} -o ${OUT} -O3
echo "${MPI} ${IN} -o ${OUT}"
##${RUN} ${THREAD} ${numnp} ${OUT}
#exit 0
#---------------------------------

#---run the code---
#for (( i=0; i<=${ng}-1; i=i+1 ))
#do
#   rgprob=${var[i*2]}
#for (( j=0; j<=${nsize}/2-1; j=j+1 ))
#do
#   lx=${var2[2*j]}
#   echo "read.in="
#   echo "${lx} ${ly} 1.0" > ${filename}
#   echo "${init} ${samples}" >> ${filename}
#   value=$(<read.in)
#   echo "${value}" 
#   unset value
#   numth=${var2[2*j+1]}
#   #${RUN} ${THREAD} ${s} ${OUT} &
#   time ${RUN} ${THREAD} ${numth} ./${OUT} &
#   #time ${RUN} ./${OUT}
#   #echo "${RUN} ./${OUT}"
#   echo "${RUN} ${THREAD} ${numth} ./${OUT} "
#   echo " ${i} "
#   PID=$!
#   wait ${PID}             #wait the previous task
#done
##done
##awk -f av3.awk chainlogcor*.dat > rcor.dat
#for (( j=0; j<=${nsize}/2-1; j=j+1 ))
#do
#   lx=${var2[2*j]}
#   echo "read.in="
#   echo "${lx} ${ly} 0.98" > ${filename}
#   echo "${init} ${samples}" >> ${filename}
#   value=$(<read.in)
#   echo "${value}" 
#   unset value
#   numth=${var2[2*j+1]}
#   #${RUN} ${THREAD} ${s} ${OUT} &
#   time ${RUN} ${THREAD} ${numth} ./${OUT} &
#   #time ${RUN} ./${OUT}
#   #echo "${RUN} ./${OUT}"
#   echo "${RUN} ${THREAD} ${numth} ./${OUT} "
#   echo " ${i} "
#   PID=$!
#   wait ${PID}             #wait the previous task
#done
#for (( j=0; j<=${nsize}/2-1; j=j+1 ))
#do
#   lx=${var2[2*j]}
#   echo "read.in="
#   echo "${lx} ${ly} 0.96" > ${filename}
#   echo "${init} ${samples}" >> ${filename}
#   value=$(<read.in)
#   echo "${value}" 
#   unset value
#   numth=${var2[2*j+1]}
#   #${RUN} ${THREAD} ${s} ${OUT} &
#   time ${RUN} ${THREAD} ${numth} ./${OUT} &
#   #time ${RUN} ./${OUT}
#   #echo "${RUN} ./${OUT}"
#   echo "${RUN} ${THREAD} ${numth} ./${OUT} "
#   echo " ${i} "
#   PID=$!
#   wait ${PID}             #wait the previous task
#done
#for (( j=0; j<=${nsize}/2-1; j=j+1 ))
#do
#   lx=${var2[2*j]}
#   echo "read.in="
#   echo "${lx} ${ly} 0.94" > ${filename}
#   echo "${init} ${samples}" >> ${filename}
#   value=$(<read.in)
#   echo "${value}" 
#   unset value
#   numth=${var2[2*j+1]}
#   #${RUN} ${THREAD} ${s} ${OUT} &
#   time ${RUN} ${THREAD} ${numth} ./${OUT} &
#   #time ${RUN} ./${OUT}
#   #echo "${RUN} ./${OUT}"
#   echo "${RUN} ${THREAD} ${numth} ./${OUT} "
#   echo " ${i} "
#   PID=$!
#   wait ${PID}             #wait the previous task
#done
#done
#awk -f av3.awk chainlogcor*.dat > rcor.dat
for (( j=0; j<=${nsize}/2-1; j=j+1 ))
do
   lx=${var2[2*j]}
   echo "read.in="
   echo "${lx} ${ly} 1.0" > ${filename}
   echo "${init} ${samples}" >> ${filename}
   value=$(<read.in)
   echo "${value}" 
   unset value
   numth=${var2[2*j+1]}
   #${RUN} ${THREAD} ${s} ${OUT} &
   time ${RUN} ${THREAD} ${numth} ./${OUT} &
   #time ${RUN} ./${OUT}
   #echo "${RUN} ./${OUT}"
   echo "${RUN} ${THREAD} ${numth} ./${OUT} "
   echo " ${i} "
   PID=$!
   wait ${PID}             #wait the previous task
done
for (( j=0; j<=${nsize}/2-1; j=j+1 ))
do
   lx=${var2[2*j]}
   echo "read.in="
   echo "${lx} ${ly} 0.99" > ${filename}
   echo "${init} ${samples}" >> ${filename}
   value=$(<read.in)
   echo "${value}" 
   unset value
   numth=${var2[2*j+1]}
   #${RUN} ${THREAD} ${s} ${OUT} &
   time ${RUN} ${THREAD} ${numth} ./${OUT} &
   #time ${RUN} ./${OUT}
   #echo "${RUN} ./${OUT}"
   echo "${RUN} ${THREAD} ${numth} ./${OUT} "
   echo " ${i} "
   PID=$!
   wait ${PID}             #wait the previous task
done
for (( j=0; j<=${nsize}/2-1; j=j+1 ))
do
   lx=${var2[2*j]}
   echo "read.in="
   echo "${lx} ${ly} 0.98" > ${filename}
   echo "${init} ${samples}" >> ${filename}
   value=$(<read.in)
   echo "${value}" 
   unset value
   numth=${var2[2*j+1]}
   #${RUN} ${THREAD} ${s} ${OUT} &
   time ${RUN} ${THREAD} ${numth} ./${OUT} &
   #time ${RUN} ./${OUT}
   #echo "${RUN} ./${OUT}"
   echo "${RUN} ${THREAD} ${numth} ./${OUT} "
   echo " ${i} "
   PID=$!
   wait ${PID}             #wait the previous task
done
for (( j=0; j<=${nsize}/2-1; j=j+1 ))
do
   lx=${var2[2*j]}
   echo "read.in="
   echo "${lx} ${ly} 0.94" > ${filename}
   echo "${init} ${samples}" >> ${filename}
   value=$(<read.in)
   echo "${value}" 
   unset value
   numth=${var2[2*j+1]}
   #${RUN} ${THREAD} ${s} ${OUT} &
   time ${RUN} ${THREAD} ${numth} ./${OUT} &
   #time ${RUN} ./${OUT}
   #echo "${RUN} ./${OUT}"
   echo "${RUN} ${THREAD} ${numth} ./${OUT} "
   echo " ${i} "
   PID=$!
   wait ${PID}             #wait the previous task
done
for (( j=0; j<=${nsize}/2-1; j=j+1 ))
do
   lx=${var2[2*j]}
   echo "read.in="
   echo "${lx} ${ly} 0.9" > ${filename}
   echo "${init} ${samples}" >> ${filename}
   value=$(<read.in)
   echo "${value}" 
   unset value
   numth=${var2[2*j+1]}
   #${RUN} ${THREAD} ${s} ${OUT} &
   time ${RUN} ${THREAD} ${numth} ./${OUT} &
   #time ${RUN} ./${OUT}
   #echo "${RUN} ./${OUT}"
   echo "${RUN} ${THREAD} ${numth} ./${OUT} "
   echo " ${i} "
   PID=$!
   wait ${PID}             #wait the previous task
done
for (( j=0; j<=${nsize}/2-1; j=j+1 ))
do
   lx=${var2[2*j]}
   echo "read.in="
   echo "${lx} ${ly} 0.8" > ${filename}
   echo "${init} ${samples}" >> ${filename}
   value=$(<read.in)
   echo "${value}" 
   unset value
   numth=${var2[2*j+1]}
   #${RUN} ${THREAD} ${s} ${OUT} &
   time ${RUN} ${THREAD} ${numth} ./${OUT} &
   #time ${RUN} ./${OUT}
   #echo "${RUN} ./${OUT}"
   echo "${RUN} ${THREAD} ${numth} ./${OUT} "
   echo " ${i} "
   PID=$!
   wait ${PID}             #wait the previous task
done
exit 0
#------------------

