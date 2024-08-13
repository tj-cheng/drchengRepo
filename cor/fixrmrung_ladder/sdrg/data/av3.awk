#{a[$1]+=exp($4); a2[$1]+=(exp($4)*exp($4)); n[$1]++}
{a[$1]+=$3; a2[$1]+=($3*$3); n[$1]++; prob[$1]=$4}
END{
for(i=8; i<=10000; i=i+2)
 if(n[i]) 
 {
    if(prob[i] == 8) 
    {
       aa=-a[i]/n[i]; aa2=a2[i]/n[i]; 
       e1=sqrt((aa2-aa*aa)/n[i]);; 
       #printf("%d\t%.16f\t%.16f\t%d\n", i, 0.5*(aa+bb), sqrt(0.5*(e1*e1+e2*e2)), 2*n[i])
       printf("%d\t%.16f\t%.16f\t%d\t%.16f\n", i, sqrt(aa*aa), sqrt(e1*e1), n[i], prob[i])
    } 
 }
 }
