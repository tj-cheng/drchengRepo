#!/bin/bash

#initial setting
filename='read.in'
lx=8
ly=2
dd=1
mm0=12
init=0
samples=2
mstps=100
istps=2
echo "${lx} ${ly} ${dd} ${mm0}" > ${filename}
echo "${init} ${samples} ${mstps} ${istps}" >> ${filename}  
unset value
value=$(<read.in)
echo "${value}" 
echo '----------'
runfile=$(cat run.in)
#echo ${runfile}
di=0
ds=0
for variable in ${runfile}
do                      #save the variable in array var[]
   var[di]=${variable}
   di=$(( ${di}+1 ))
   ds=$(( ${ds}+1 ))
done
ng=$(( ${ds} / 2 ))   #number of lines in run.in
nv=$(( ${ds} / ${ng} ))   #number of columns in run.in
#echo $(( ${ds} / 3 ))
#echo ${nv}
#echo "lx=${var[3]}" #lx_1
#echo ${var[3]} #lx_2
unset value 
#-------------------------------------

#exit 0
# >> append ; > new to write

#compile mpif90 -o cat.out cat.f90
IN=${1}

#read -p " please keyin the number of threads = " numnp

MPI=mpif90

RUN=mpirun

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
   echo "${lx} ${ly} ${dd} ${mm0}" > ${filename}
   echo "${init} ${samples} ${mstps} ${istps}" >> ${filename}  
   value=$(<read.in)
   echo "${value}" 
   unset value
   numth=${var[i*2+1]}
   #${RUN} ${THREAD} ${s} ${OUT} &
   ${RUN} ${THREAD} ${numth} ${OUT} &
   echo "${RUN} ${THREAD} ${numth} ${OUT} "
   echo " ${i} "
   PID=$!
   wait ${PID}             #wait the previous task
done
exit 0
#------------------

