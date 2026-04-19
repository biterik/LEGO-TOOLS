#!/usr/bin/awk -f
#
# Author:  Erik Bitzek <e.bitzek@mpi-susmat.de>
#          Max-Planck-Institut fuer Nachhaltige Materialien,
#          Duesseldorf, Germany
# Funding: NFDI-MatWerk
BEGIN {
OFMT = "%.16f";
CONVFMT ="%.16f";
xmin=ARGV[1];
ARGV[1]="";
xmax=ARGV[2];
ARGV[2]="";
ymin=ARGV[3];
ARGV[3]="";
ymax=ARGV[4];
ARGV[4]="";
zmin=ARGV[5];
ARGV[5]="";
zmax=ARGV[6];
ARGV[6]="";
}

{}{
if($0 !~ /#/)
{
  if(($4>=xmin)&&($4<=xmax)&&($5>=ymin)&&($5<=ymax)&&($6>=zmin)&&($6<=zmax))
    print $0;

}
else
{print $0;}
}
END{}



