{a[$1]+=sqrt($2*$2); a2[$1]+=($2*$2); n[$1]++}
END{
for(i=8; i<=300; i=i+2)
 if(n[i]) 
   {
     aa=a[i]/n[i]; aa2=a2[i]/n[i]; 
     e1=sqrt((aa2-aa*aa)/n[i]); 
     #printf("%d\t%.16f\t%.16f\t%d\n", i, 0.5*(aa+bb), sqrt(0.5*(e1*e1+e2*e2)), 2*n[i])
     printf("%d\t%.16f\t%.16f\t%d\n", i, aa, sqrt(e1*e1), n[i])
   } 
}
