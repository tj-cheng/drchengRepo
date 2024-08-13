{a[$3]+=$5; a2[$3]+=($5*$5); b[$3]+=$6; b2[$3]+=(($6)*($6)); n[$3]++; size=$1}#; prob=$8}
END{
for(i=1; i<=2048; i=i+1)
 if(n[i]) 
   {
     aa=a[i]/n[i]; aa2=a2[i]/n[i]; bb=b[i]/n[i]; bb2=b2[i]/n[i]; 
     e1=sqrt((aa2-aa*aa)/n[i]); e2=sqrt((bb2-bb*bb)/n[i]); 
     #printf("%d\t%.30f\t%.30f\t%d\t%d\t%.3f\n", i, bb, e2 , n[i],size,prob) #locsus
     printf("%d\t%.30f\t%.30f\t%d\t%d\n", i, bb, e2 , n[i],size) #locsus
   } 
}
