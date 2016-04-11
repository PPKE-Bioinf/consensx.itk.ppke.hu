#!/usr/bin/perl

print STDERR"
######################################################################
# disre2stereoambig_xplor.pl < disre_in > disre_out                  #
#                                                                    #
# Takes a disre file with X-PLOR atom nomenclature and symmetrizes   #
# stereospecifically given restraints as:                            #
# - HB1 -> HB1 or HB2 etc.                                           #
# - if distance for HB1 and HB2 is different, assigns the larger to  #
#   both                                                             #
# The program works for the atom nomenclature of X-PLOR, make sure   #
# to convert your input if necessary!                                #
#                                                                    #
# Version 3.2                                            24.03.2016. #
######################################################################

";

use Getopt::Std;
$|=1;

$opt_b='';
getopts('b:ohv')  || die ("Try again please\n");
exit if $opt_h;

$ORONLY=$opt_o;

if ($opt_b ne ""){
    @BINS=(sort {$a <=> $b} (split(/\:/,$opt_b)));
}
#exit;


# The following table contains atoms and atom groups that can
# be affected by non-stereospecific assignment or other ambiguity usually not evident
# at the stage of resonance assignment / NOE identification  
# Each atom group is assigned to all relevant amino acids

%AMBIG=("|GLY|",                                                        "HA1|HA2",
        "|ARG|ASP|ASN|CYS|GLU|GLN|HIS|LEU|LYS|MET|PHE|PRO|SER|TRP|TYR|","HB1|HB2",
        "|ARG|GLU|GLN|LYS|MET|PRO|",                                    "HG1|HG2",
        "|ARG|LYS|PHE|PRO|TYR|",                                        "HD1|HD2",
        "|LYS|PHE|TYR|",                                                "HE1|HE2",
        "|ILE|",                                                        "HG11|HG12", # added: Mar 24 2016
        "|VAL|",                                                        "HG11|HG12|HG13|HG21|HG22|HG23",
        "|LEU|",                                                        "HD11|HD12|HD13|HD21|HD22|HD23",
        "|ARG|",                                                        "HH11|HH12|HH21|HH22",
        "|ASN|",                                                        "HD21|HD22",
        "|GLN|",                                                        "HE21|HE22"
    );


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
    # Per-pair, per-atom and per-residue info
    if ($DRFORPAIR{$pair} !~ /[0-9A-Z]/){
	$DRFORPAIR{$pair}=$dn;
    }
    else{
	$DRFORPAIR{$pair}.=",".$dn;
	print STDERR "$pair already present in restraints: $DRFORPAIR{$pair}\n";
	$NONUNIQUEPAIRS++;
    }


    $DIST[$dn]=$dist;
    # recording dn last seen
    $prevdn=$dn;
}
$maxdn=$dn;
# For now I assume that NONUNIQUEPAIRS == 0 
$dnout=0;
#foreach $dn (sort {$DRPR{$a} <=> $DRPR{$b}} keys %DRPR){
for ($dn = 0; $dn<=$maxdn; $dn++){

    print STDERR "\n>> $dn << DIST $DIST[$dn] <<\n$PAIRSFORDR[$dn]\n" if $opt_v;
    $prev_NOP=0;
    @pfd=split(/\,/,$PAIRSFORDR[$dn]);
    $NOP=@pfd;
    # Repeat check & add until no pairs can be added
    while ($NOP > $prev_NOP){
     foreach $pair (sort @pfd){
        ($atomid1,$atomid2)=split(/\|/,$pair);
        ($rnum1,$atom1)=split(/\_/,$atomid1);
        ($rnum2,$atom2)=split(/\_/,$atomid2);
	print STDERR " pair $pair -> $rnum1 $RESNAME{$rnum1} $atom1 | $rnum2 $RESNAME{$rnum2} $atom2 \n" if $opt_v;

	if (&listedambig($RESNAME{$rnum1},$atom1)){

	    foreach $aa (sort @AMBIGATOMS){
		print STDERR "-1-> $aa\n" if $opt_v;
		$ambigatomid1=$rnum1."_".$aa;
		if ($ambigatomid1 eq $atomid2){last;}
		$ambigpair=$ambigatomid1."|".$atomid2;
		if ($ambigpair ne $pair){
		    if ($DRFORPAIR{$ambigpair} eq $dn){
			print STDERR "$ambigpair in self\n" if $opt_v;
		    }
		    elsif ($DRFORPAIR{$ambigpair} =~ /[0-9]/){
			print STDERR "$ambigpair in $DRFORPAIR{$ambigpair} DIST $DIST[$DRFORPAIR{$ambigpair}]\n" if $opt_v;
			&makedistancesequal($dn, $DRFORPAIR{$ambigpair});
                        &addpair($dn, $ambigpair,0); 
		    }
		    else{
			print STDERR "$ambigpair not found\n" if $opt_v;
                        &addpair($dn, $ambigpair,1); 
		    }
		}
	    }#_foreach $aa

	}#_if atom1 ambig

	if (&listedambig($RESNAME{$rnum2},$atom2)){
	    foreach $aa (sort @AMBIGATOMS){
		print STDERR "-2-> $aa\n" if $opt_v;
		$ambigatomid2=$rnum2."_".$aa;
		if ($ambigatomid2 eq $atomid1){last;}
		$ambigpair=$atomid1."|".$ambigatomid2;
		if ($ambigpair ne $pair){
		    if ($DRFORPAIR{$ambigpair} eq $dn){
			print STDERR "$ambigpair in self\n" if $opt_v;
		    }
		    elsif ($DRFORPAIR{$ambigpair} =~ /[0-9]/){
			print STDERR "$ambigpair in $DRFORPAIR{$ambigpair} DIST $DIST[$DRFORPAIR{$ambigpair}]\n" if $opt_v;
			&makedistancesequal($dn, $DRFORPAIR{$ambigpair});
                        &addpair($dn, $ambigpair,0); 
		    }
		    else{
			print STDERR "$ambigpair not found\n" if $opt_v;
                        &addpair($dn, $ambigpair,1); 
		    }
		}
		else{$AmbigPairInDR{$dn}++;}
	    }#_foreach $aa

	}# atom2 ambig

     
     }#_foreach pair
     $prev_NOP=$NOP;
     @pfd=split(/\,/,$PAIRSFORDR[$dn]);
     $NOP=@pfd;
    }# while NOP > prev_NOP (pairs have been added) 


# writing updated restraint
    print "$AREMARK[$dn]$DREMARK[$dn]\n" if (($AREMARK[$dn] =~ /[a-z]/i) || ($DREMARK[$dn] =~ /[a-z]/i)) ;

    @pfd=split(/\,/,$PAIRSFORDR[$dn]);
    foreach $pair (sort @pfd){
        ($atomid1,$atomid2)=split(/\|/,$pair);
        ($rnum1,$atom1)=split(/\_/,$atomid1);
        ($rnum2,$atom2)=split(/\_/,$atomid2);
	printf (" %4d   %3d %3s %4s   %3d %3s %4s   %5.2f\n",$dnout,$rnum1,$RESNAME{$rnum1},$atom1,$rnum2,$RESNAME{$rnum2},$atom2,$DIST[$dn]);
    }
    $dnout++;
}

# test so far
exit;


sub listedambig{
    my $listed=0;
    $resname=shift;
    $atomname=shift;
    @AMBIGATOMS=();
    foreach $ambigres (keys %AMBIG){
       if (($ambigres =~ /\|$resname\|/) && ($AMBIG{$ambigres} =~ /$atomname/)){
	   $listed=1;
	   @AMBIGATOMS=split(/\|/,$AMBIG{$ambigres});
       }  
    }
    return $listed;
}

# Setting the lower distance equal to the larger one for tho restraints
sub makedistancesequal{

    $dn1=shift;
    $dn2=shift;

    if ($DIST[$dn1] < $DIST[$dn2]){
	$DIST[$dn1]=$DIST[$dn2];
	$DREMARK[$dn1].="# dist set as in $dn2 ";
    }

}

# Add a resstrained atom pair to a restraint
sub addpair{

    my $temp_dn=shift;
    my $temp_pair=shift;
    my $addfp=shift;
    $pairfound=0;
    foreach $pf (split(/\,/,$PAIRSFORDR[$temp_dn])){
	if ($pf eq $temp_pair){$pairfound=1}
    }
    if (!$pairfound){
	$PAIRSFORDR[$temp_dn].=",".$temp_pair;
       if ($addfp){
         if ($DRFORPAIR{$temp_pair} !~ /[0-9]/){
	  $DRFORPAIR{$temp_pair}.=$temp_dn;
         }
         else{
	  $DRFORPAIR{$temp_pair}.=",".$temp_dn;
         }
       }
	$AREMARK[$temp_dn].="# $temp_pair added ";
    }#_if pair not found


}


