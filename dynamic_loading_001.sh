for ii in $(seq 1 11)
do
 trial=$(echo "scale=0; (100000+$ii);" | bc)
 modtwo=$(echo "scale=0; (($trial+1)%2);" | bc)
 modthree=$(echo "scale=0; (($trial+2)%3);" | bc)
 modfour=$(echo "scale=0; (($trial+3)%4);" | bc)
 modfive=$(echo "scale=0; (($trial+4)%5);" | bc)
 modsix=$(echo "scale=0; (($trial+5)%6+1);" | bc)
 modseven=$(echo "scale=0; (($trial+6)%7);" | bc)
 modnine=$(echo "scale=0; (($trial+8)%9);" | bc)
 modten=$(echo "scale=0; (($trial+9)%10);" | bc)
 modeleven=$(echo "scale=0; (($trial+10)%11);" | bc)
 modthirteen=$(echo "scale=0; (($trial+12)%13);" | bc)
 modfourteen=$(echo "scale=0; (($trial+13)%14);" | bc)
 modfifteen=$(echo "scale=0; (($trial+14)%15);" | bc)
 modnineteen=$(echo "scale=0; (($trial+18)%19);" | bc)
 spring=$(echo "scale=8; (200.);" | bc)
 pexcvol=$spring
 force=$(echo "scale=8; (0.0);" | bc)
 lcontrol=$(echo "scale=0; (0);" | bc)
 pstiffness=$(echo "scale=2; (100.);" | bc)
 pvelocity=$(echo "scale=8; (0.0001);" | bc)
 radius=$(echo "scale=8; (7.07106781);" | bc)
 seed=$(echo "scale=0; (2+$modeleven);" | bc)
 lowbound=$(echo "scale=0; (4);" | bc)
 upbound=$(echo "scale=0; (8);" | bc)
 therm=$(echo "scale=0; (1);" | bc)
 ldep=$(echo "scale=0; (0);" | bc)
 springsonly=$(echo "scale=0; (0);" | bc)
 variablelength=$(echo "scale=0; (1);" | bc)
 boxlength=$(echo "scale=4; (4.1*$radius+0.35*$modthirteen*$radius);" | bc)
 polylength=$(echo "scale=0; (552);" | bc)
 steps=$(echo "scale=0; (20000001);" | bc)
 nshell=$(echo "scale=0; (1000);" | bc)
 springc=$(echo "scale=1; (100.);" | bc)
 shellbond=$(echo "scale=1; (100.);" | bc)
 loadfrac=$(echo "scale=6; (0.001);" | bc)
 kbt=$(echo "scale=8; (0.001);" | bc)
 dt=$(echo "scale=8; (0.0005);" | bc)
 pmonorad=$(echo "scale=2; (0.4);" | bc)
 crosslinkdensity=$(echo "scale=4; (0.2);" | bc)
 shelllinks=$(echo "scale=0; (40);" | bc)
 restartbool=$(echo "scale=0; (0);" | bc)
./runnucleus -t $trial -f $force -b $spring -v $pexcvol -r $radius -z $seed -l $lowbound -T $therm -L $ldep -u $upbound -P $springsonly -V $variablelength -x $boxlength -N $polylength -n $steps -M $nshell -k $springc -a $shellbond -F $loadfrac -y $kbt -d $dt -q $lcontrol -g $pstiffness -h $pvelocity -w $pmonorad -c $crosslinkdensity -C $shelllinks -R $restartbool  
done