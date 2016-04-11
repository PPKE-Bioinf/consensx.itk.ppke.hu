#!/usr/bin/perl

print STDERR"
######################################################################
# disre_sort.pl -b<bins>  < disre_in > sorted_disre_out              #
#                                                                    #
# Replaces distances with that of the closest bin defined.           #
#                                                                    #
# Version 0.8                                            17.01.2016. #
######################################################################

";

use Getopt::Std;
$|=1;

$opt_b="2.5:3.5:5:6";
getopts('b:hv')  || die ("Try again please\n");
exit if $opt_h;

if ($opt_b ne ""){
    @BINS=(sort {$a <=> $b} (split(/\:/,$opt_b)));
}
print "# Restraint distances binned as:\n";
$prevd=0;
foreach $b (@BINS){
    printf ("# %5.2f - %5.2f\n",$prevd,$b);
    $prevd=$b;
}

# Reading disre file
$prevdn=-1;

while(<>){
    chomp;
    if ($_ =~ /^ +[0-9]+ /){
       $_=~s/^ +//; 
       ($dn,$rnum1,$rname1,$atom1,$rnum2,$rname2,$atom2,$dist)=split(/ +/,$_);
        foreach $b (@BINS){
           if ($dist < $b){$dist=$b;last}
        }
	printf (" %4d   %3d %3s %4s   %3d %3s %4s   %5.2f\n",$dn,$rnum1,$rname1,$atom1,$rnum2,$rname2,$atom2,$dist);
    } 
    else{
	print "$_\n";
    }
}
