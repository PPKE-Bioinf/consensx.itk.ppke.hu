#!/usr/bin/perl

print STDERR"
######################################################################
# disre_renum.pl -o<offset>  < disre_in > renumbered_disre_out       #
#                                                                    #
# Adds the offset to the residue numbers in the input disre file     #
#                                                                    #
# Version 0.8                                            17.01.2016. #
######################################################################

";

use Getopt::Std;
$|=1;

$opt_o=0;
getopts('o:hv')  || die ("Try again please\n");
exit if $opt_h;

if ($opt_o == 0){
    die "No or zero offset given: nothing to do, exiting (won't act as 'cp)\n";
}

# Reading disre file
print "# Renumbered by $opt_o\n";
while(<>){
    chomp;
    if ($_ =~ /^ +[0-9]+ /){
       $_=~s/^ +//; 
       ($dn,$rnum1,$rname1,$atom1,$rnum2,$rname2,$atom2,$dist)=split(/ +/,$_);
       printf (" %4d   %3d %3s %4s   %3d %3s %4s   %5.2f\n",$dn,$rnum1+$opt_o,$rname1,$atom1,$rnum2+$opt_o,$rname2,$atom2,$dist);
    } 
    else{
	print "$_\n";
    }
}
