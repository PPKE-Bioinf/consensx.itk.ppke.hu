#!/usr/bin/perl

print STDERR"
######################################################################
# disre_removeambig.pl [-r][-t][-a] < disre_in > disre_out           #
#                                                                    #
# Removes restraints whith ambiguities.                              #
# -r causes to remove those with residue-level ambiguities.          #
# -t causes to remove ambiguities of hydrogens attached to carbon    #
#    at a different topological distance from CA (like HB/HG/HD/)    #
# -a causes to remove all ambiguities, i.e. keeps only those with a  #
#    single well-defined atom pair.                                  #
#                                                                    #
# Version 0.2                                            30.03.2016. #
######################################################################

";

use Getopt::Std;
$|=1;

getopts('rtahv')  || die ("Try again please\n");
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
    $atomstem1=$atom1;$atomstem1=~s/[0-9]//g;
    $atomstem2=$atom2;$atomstem2=~s/[0-9]//g;

#print STDERR "$atom1 -> $atomstem1\n";
    # Even if it is a single pair, there should be no amibuities
    # Pseudoatoms are not allowed in disre files, but just in case. 
    if (($atom1=~/[MQ\#\%\+]/) || ($atom2=~/[MQ\#\%\+]/)){
	    $ANYAMBIG{$dn}=1 if $opt_a;
    }

    # Identifying restraints with ambiguity of given type
    if ($prevdn == $dn){
        # Making use of ordered storing
	if (($prevrnum1 != $rnum1) || ($prevrnum2 != $rnum2)){
	    $RESAMBIG{$dn}=1 if $opt_r;
	}
	if (($prevatomstem1 ne $atomstem1) || ($prevatomstem2 ne $atomstem2)){
	    $ATSAMBIG{$dn}=1 if $opt_t;
	}
	if ($opt_a){
	    $ANYAMBIG{$dn}=1;
	}

    }

    $prevdn=$dn;
    $prevrnum1=$rnum1;
    $prevrnum2=$rnum2;
    $prevatomstem1=$atomstem1;
    $prevatomstem2=$atomstem2;
}

$dnout=0;
$maxdn=$dn;
#foreach $dn (sort {$a <=> $b} keys %DRSREAD){
for ($dn = 0; $dn<=$maxdn; $dn++){

    if ((!$RESAMBIG{$dn}) && (!$ATSAMBIG{$dn}) && (!$ANYAMBIG{$dn})){ 

     print "# was $dn\n";
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
    elsif ($RESAMBIG{$dn}){
	print "# Restraint $dn contains ambiguous residue assignment\n";
    }
    elsif ($ATSAMBIG{$dn}){
	print "# Restraint $dn contains ambiguous atom name assignment\n";
    }
    else{
	print "# Restraint $dn contains other type of ambiguitiy\n";
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



