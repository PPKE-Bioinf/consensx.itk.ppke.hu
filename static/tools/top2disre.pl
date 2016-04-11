#!/usr/bin/perl

print STDERR"
##########################################################
# disre2tbl.pl < gromacs_topology_file > xplor_tbl_file  #
#                                                        #
# An attempt to generate an x-plor format distance       #
# restraint file from the remarks in the distance        #
# restraint section of a gromacs topology file.          #
#                                                        #
##########################################################

";

use Getopt::Std;
getopts ('h') || die;
exit if $opt_h;

$readdisre=0;
while(<>){
    chomp;
    next if (($_ =~/^ *\;/) || ($_ !~/[a-z]/i));
    if ($_=~/\[ *distance_restraints *\]/){$readdisre=1}
    elsif ($_=~/\[ atoms \]/){$readatoms=1}
    elsif ($_=~/^ *\[/){$readdisre=0;$readatoms=0}
    if (($readatoms) && ($_ =~ /^ *[0-9]/)){
	$atomline=$_;
	$atomline=~s/^ +//;
	$rn=(split(/ +/,$atomline))[2];
	$RNAME{$rn}=(split(/ +/,$atomline))[3];
    }
    if (($readdisre) && ($_ =~ /^ *[0-9]/)){
	next if ($_=~/^\;/);
        #$_=~s/^ +//;
        #@allcol=split(/[ \;\:]+/,$_);
	$dline=$_;
	$dline=~s/\;.*$//;$dline=~s/^ +//;
        ($ai,$aj,$type,$index,$typea,$low,$up1,$up2,$fac)=split(/ +/,$dline);
        $remark=$_;
        $remark=~s/^.*\; *//;
        ($ra1,$ra2)=split(/\: +/,$remark);
	$ra1=~s/^ +//;
	$ra2=~s/ +$//;
	$atom1=substr($ra1,-4);
	$atom2=substr($ra2,-4);
	$res1=$ra1;$res1=~s/$atom1$//;
	$res2=$ra2;$res2=~s/$atom2$//;

	printf("%3d %4d %3s $atom1  %4d %3s $atom2 %5.2f\n",$index,$res1,$RNAME{$res1},$res2,$RNAME{$res2},$up1*10);
    }
}

