#!/usr/bin/perl

print STDERR"
######################################################################
# disre_sort.pl [-u] < disre_in > sorted_disre_out                   #
#                                                                    #
# The program sorts distance restraints according to residue numbers #
# and atom names.                                                    #
# -u causes to omit duplicate restraints.                            #
#                                                                    #
# Version 0.8                                            17.01.2016. #
######################################################################

";

use Getopt::Std;
$|=1;

getopts('uhv')  || die ("Try again please\n");
exit if $opt_h;


# Reading disre file
$prevdn=-1;

while(<>){
    chomp;
    next if ($_=~/^\#/);
    $_=~s/\#.*$//;
    $_=~s/^ +//;
    ($dn,$rnum1,$rname1,$atom1,$rnum2,$rname2,$atom2,$dist)=split(/ +/,$_);
    # Will store everything ordered: first by residues, if equal, by atoms 
    if (($rnum2 < $rnum1) || (($rnum1 == $rnum2) && ($atom2 lt $atom1))){
	$atomid1=$rnum2."_".$atom2;
	$atomid2=$rnum1."_".$atom1;
    }
    else{
	$atomid1=$rnum1."_".$atom1;
	$atomid2=$rnum2."_".$atom2;
    }

    $pair = $atomid1."|".$atomid2;

    # Storing pairs as a comma-separated list for DR
    if ($dn != $prevdn){
	$PAIRSFORDR[$dn]=$pair;
    }
    else{
	$PAIRSFORDR[$dn].=",".$pair;
    }


    # Checking & bookkeeping
    if (($RESNAME{$rnum1} =~ /[A-Z]/) && ($RESNAME{$rnum1} ne $rname1)){
	die "In DR $dr residue $rnum1 is a $rname1 but was prevously identified a $RESNAME{$rnum1}\n";
    }
    else{
	$RESNAME{$rnum1}=$rname1;
    }
    if (($RESNAME{$rnum2} =~ /[A-Z]/) && ($RESNAME{$rnum2} ne $rname2)){
	die "In DR $dr residue $rnum2 is a $rname2 but was prevously identified a $RESNAME{$rnum2}\n";
    }
    else{
	$RESNAME{$rnum2}=$rname2;
    }

    $DIST[$dn]=$dist;

    $DRSREAD{$dn}=1;

    $prevdn=$dn;
}

$dnout=0;
foreach $dn (sort { disrepriority($a,$b) } keys %DRSREAD){
#                     if (disrepriority($a,$b) == $a){return -1}
#                     else {return 1}
#                   } keys %DRSREAD){


    if ((!$WRITTEN{$PAIRSFORDR[$dn]}) || (!$opt_u)){ 

     print "# was $dn";
     if ($WRITTEN{$PAIRSFORDR[$dn]}){
        print " # Duplicate of $WRITTEN{$PAIRSFORDR[$dn]}";
     }
     print "\n";
     @pfd=split(/\,/,$PAIRSFORDR[$dn]);
     foreach $pair (sort @pfd){
        ($atomid1,$atomid2)=split(/\|/,$pair);
        ($rnum1,$atom1)=split(/\_/,$atomid1);
        ($rnum2,$atom2)=split(/\_/,$atomid2);
	printf (" %4d   %3d %3s %4s   %3d %3s %4s   %5.2f\n",$dnout,$rnum1,$RESNAME{$rnum1},$atom1,$rnum2,$RESNAME{$rnum2},$atom2,$DIST[$dn]);
     }
     $WRITTEN{$PAIRSFORDR[$dn]}=$dnout unless ($WRITTEN{$PAIRSFORDR[$dn]});
     $dnout++;
    }
    else{
	print "# Restraint $dn is a duplicate of $WRITTEN{$PAIRSFORDR[$dn]}\n";
    }
}

# This function returns -1 if restraint a has to come before restraint b, 1 otherwise
sub disrepriority{

    my $a=shift;
    my $b=shift;
    @Pa=split(/\,/,$PAIRSFORDR[$a]);
    @Pb=split(/\,/,$PAIRSFORDR[$b]);
    $minrn_a1=(split(/[\|\_]/,$Pa[0]))[0];
    $minan_a1=(split(/[\|\_]/,$Pa[0]))[1];
    $minrn_a2=(split(/[\|\_]/,$Pa[0]))[2];
    $minan_a2=(split(/[\|\_]/,$Pa[0]))[3];
    $minrn_b1=(split(/[\|\_]/,$Pb[0]))[0];
    $minan_b1=(split(/[\|\_]/,$Pb[0]))[1];
    $minrn_b2=(split(/[\|\_]/,$Pb[0]))[2];
    $minan_b2=(split(/[\|\_]/,$Pb[0]))[3];
    #print "$a $minrn_a1 $minrn_a2 | $b $minrn_b1 $minrn_b2\n";
    if ($minrn_a1 < $minrn_b1){return -1}
    elsif(($minrn_a1 == $minrn_b1) && ($minrn_a2 < $minrn_b2)){return -1}
    elsif(($minrn_a1 == $minrn_b1) && ($minrn_a2 == $minrn_b2) && ($minan_a1 lt $minan_b1)){return -1}
    elsif(($minrn_a1 == $minrn_b1) && ($minrn_a2 == $minrn_b2) && ($minan_a1 eq $minan_b1) && ($minan_a2 lt $minan_b2)){return -1}
    else{
	return 1;
    }
}



