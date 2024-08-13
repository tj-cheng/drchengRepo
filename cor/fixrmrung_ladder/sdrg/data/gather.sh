awk '{if($4 == 1) print $0}' 2legladder*.dat > dp100.dat
awk '{if($4 == 0.98) print $0}' 2legladder*.dat > dp098.dat
awk '{if($4 == 0.96) print $0}' 2legladder*.dat > dp096.dat
awk '{if($4 == 0.94) print $0}' 2legladder*.dat > dp094.dat
awk '{if($4 == 0.92) print $0}' 2legladder*.dat > dp092.dat
awk '{if($4 == 0.9) print $0}' 2legladder*.dat > dp090.dat
awk '{if($4 == 0.8) print $0}' 2legladder*.dat > dp080.dat
awk '{if($4 == 0.88) print $0}' 2legladder*.dat > dp088.dat
awk '{if($4 == 0.7) print $0}' 2legladder*.dat > dp070.dat
awk '{if($4 == 0.5) print $0}' 2legladder*.dat > dp050.dat
awk '{if($4 == 0) print $0}' 2legladder*.dat > dp000.dat
#awk -f av3.awk dp100.dat 
#awk -f av3.awk dp098.dat 
#awk -f av3.awk dp096.dat 
#awk -f av3.awk dp094.dat 
#awk -f av3.awk dp092.dat 
#awk -f av3.awk dp090.dat 
#awk -f av3.awk dp080.dat 
#awk -f av3.awk dp088.dat 
#awk -f av3.awk dp070.dat 
#awk -f av3.awk dp050.dat 
#awk -f av3.awk dp000.dat 
awk -f av3.awk dp100.dat > rcorp100.dat
awk -f av3.awk dp098.dat > rcorp098.dat
awk -f av3.awk dp096.dat > rcorp096.dat
awk -f av3.awk dp094.dat > rcorp094.dat
awk -f av3.awk dp092.dat > rcorp092.dat
awk -f av3.awk dp090.dat > rcorp090.dat
awk -f av3.awk dp088.dat > rcorp088.dat
awk -f av3.awk dp080.dat > rcorp080.dat
awk -f av3.awk dp070.dat > rcorp070.dat
awk -f av3.awk dp050.dat > rcorp050.dat
awk -f av3.awk dp000.dat > rcorp000.dat
