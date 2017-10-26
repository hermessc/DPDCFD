set -x 

rm -r 0.*
rm -r 1*
rm -r 2*
rm -r 3*
rm -r 5*
rm -r 4*
rm -r 6*
rm -r processor*
cd 0/
mv U.bkp U; cp U U.bkp
mv p.bkp p; cp p p.bkp
mv shearRate.bkp sharRate; cp shearRate shearRate.bkp
mv viscoDPD.bkp viscoDPD; cp viscoDPD viscoDPD.bkp
cd ..



