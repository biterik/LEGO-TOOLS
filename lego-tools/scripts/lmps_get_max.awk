#!/usr/bin/awk -f
# adjusted for lammps format: nr typ  x y z vx vy vz
#
# Author:  Erik Bitzek <e.bitzek@mpi-susmat.de>
#          Max-Planck-Institut fuer Nachhaltige Materialien,
#          Duesseldorf, Germany
# Funding: NFDI-MatWerk

BEGIN {
OFMT = "%.16g";
CONVFMT ="%.16g";
xmax=0;ymax=0;zmax=0;epotmax=0;epotmin=-10;spropmax=0;spropmin=100;
vxmax=0;vymax=0;vzmax=0;vxmin=0;vymin=0;vzmin=0;
xmin=1000;ymin=1000;zmin=1000;
Tmax=0;
i=0;
}

{}{
#if($1 !~ /#/)
if(NF==8)
{
 i++;
if($3>xmax)
	xmax=$3;
if($4>ymax)
	ymax=$4;
if($5>zmax)
	zmax=$5;
if($3<xmin)
	xmin=$3;
if($4<ymin)
	ymin=$4;
if($5<zmin)
	zmin=$5;
}
}
END{
printf(" %d Atoms\n",i);
printf("Maximum Values:\n x:%.16f y:%.16f z:%.16f \n",xmax,ymax,zmax);
printf("Minimum Values:\n x:%.16f y:%.16f z:%.16f \n",xmin,ymin,zmin);
}

#absolute value and checking of periodicity
function abs(u){
v=sqrt(u*u);
if(v<20)
return(v);
else
return(0);
}
