{
   a[$1","$4]+=$5;
   a2[$1","$4]+=($5*$5);
   #x1[$2","$3]=$4; 
   #y1[$2","$3]=$5; 
   #x2[$2","$3]=$6; 
   #y2[$2","$3]=$7; 
   n[$1","$4]++;
}
END{
   #printf("%d,%.16f",n[1","2],a[1","65]);
   for(i=1; i<=300; i=i+1)
   {
      k=-1
      for(j=1; j<=300; j=j+1)
      {

         if(n[i","j] != 0)
         {
            k=k+1
            aa=a[i","j]/n[i","j]; aa2=a2[i","j]/n[i","j]; 
            e1=sqrt((aa2-aa*aa)/n[i","j]); 
            printf("%d\t%.16f\t%.16f\t%d\t%d\t%d\n", i, aa, sqrt(e1*e1), n[i","j], j, k)
         }
         #printf("%d\t%.16f\t%.16f\t%.16f\t%d\n", i, aa, sqrt(e1*e1),prob[i], n[i])
         #dx=x1[i","j]-x2[i","j]; dy=y1[i","j]-y2[i","j];dd= dx*dx+dy*dy;
         #printf("%d\t%d\n",i,j)
         #if(int(dd) == 1 || int(dd) == 5)
         #{
         #   aa=a[i","j]/n[i","j]; aa2=a2[i","j]/n[i","j]; 
         #   e1=sqrt((aa2-aa*aa)/n[i","j]); 
         #   #printf("%d\t%.16f\t%.16f\t%.16f\t%d\n", i, aa, sqrt(e1*e1),prob[i], n[i])
         #   printf("%d\t%d\t%d\t%d\t%.16f\t%d\n", x1[i","j], y1[i","j], x2[i","j],y2[i","j], e1 ,n[i","j])
         #}
      }
   }
}
