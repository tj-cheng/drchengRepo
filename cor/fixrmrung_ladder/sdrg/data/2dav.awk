{
   a[$1","int($4)]+=$3;
   a2[$1","int($4)]+=($3*$3);
   #x1[$2","$3]=$4; 
   #y1[$2","$3]=$5; 
   #x2[$2","$3]=$6; 
   #y2[$2","$3]=$7; 
   n[$1","int($4)]++;
   prob[$1","int($4)]=$4;
}
END{
   #printf("%d,%.16f",n[1","2],a[1","65]);
   #printf("%d,%.16f",n[32","8],a[32","8]);
   for(i=1; i<=2048; i=i+1)
   {
      for(j=1; j<=100; j=j+1)
      {
         #dx=x1[i","j]-x2[i","j]; dy=y1[i","j]-y2[i","j];dd= dx*dx+dy*dy;
         #printf("%d\t%d\n",i,j)
         #if(int(dd) == 1 || int(dd) == 5)
         if(prob[i","j] == 20)
         {
            #printf(n[i","j])
            aa=a[i","j]/n[i","j]; aa2=a2[i","j]/n[i","j]; 
            e1=sqrt((aa2-aa*aa)/n[i","j]); 
            #printf("%d\t%.16f\t%.16f\t%.16f\t%d\n", i, aa, sqrt(e1*e1),prob[i], n[i])
            printf("%d\t%.16f\t%.16f\t%.16f\t%d\n", i, aa, sqrt(e1*e1),prob[i","j], n[i","j])
            #printf("%d\t%d\t%d\t%d\t%.16f\t%d\n", x1[i","j], y1[i","j], x2[i","j],y2[i","j], e1 ,n[i","j])
         }
      }
   }
}
