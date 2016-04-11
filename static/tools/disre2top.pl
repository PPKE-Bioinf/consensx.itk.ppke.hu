#!/usr/bin/perl -I/home/pinson/szpari/programs/Perl/

print STDERR "
                    :-] G R O M A C S [-:
                        helpful tools

                    :-] disre2top.pl  [-:

Option     Filename                Type          Description
-------------------------------------------------------------
";


use Getopt::Std;

$opt_d="";
$opt_p="topol.top";
$opt_o="topol_res.top";
$opt_U=99.9;
getopts('d:p:o:hv') || die "ERROR in options, exiting\n";

$disrefile=$opt_d;
if ($opt_u){$uout="yes"}
else{$uout="no"}

if ($opt_h){$HO="yes"} else{$HO="no"}
if($opt_v){$V=1;$VO="yes"}else{$VO="no"}

printf STDERR ("  -p  %20s         Input   Topology file\n",$opt_p);
printf STDERR ("  -o  %20s         Output  Topology file with NOE restraints\n",$opt_o);
printf STDERR ("  -d  %20s         Input   disre file (NOE restraints)\n",$disrefile);

print STDERR"
      Option   Type  Value  Description
------------------------------------------------------
";

printf STDERR ("          -h   bool %9s  Print help info and quit\n",$HO);
printf STDERR ("          -v   bool %9s  Be loud and noisy\n",$VO);
printf STDERR ("          -U   real %9.2f  Second upper limit (r2) for NOE restraints (nm)\n",$opt_U);
#printf STDERR ("          -u   bool %9s  Use only unambiguous restraints in NOE list\n",$uout);


if ($opt_h) {
print STDERR "
Puts NOE-based distance restraints read from a disre format file into
a GROMACS topology file.

";
exit;
}

# Read disre file


# Currently only a 'clean', one-restraint-per-line file with clearly seperable columns is read!!
open (DF,"$disrefile") || die "Cannot open disre file with NOE restraints: $disrefile\n";
$drnum=0;
#$skipped_ambig=0;
$prevdn=-1;
while(<DF>){
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

    $prevdn=$dn;
}
$drnum=$dn;

$readatoms=0;
$restraintsadded=0;
# Reading and processing topology file(s)
# Restraints will be added before the position restraint include statement
open (TOPIN,"$opt_p") || die "FATAL: Toplogy file $opt_p not found!\n";

if (-e $opt_o){
   print STDERR "Back off! Backing up $opt_o to \#$opt_o\#\n";
   system "mv $opt_o \#$opt_o\#";
}
open (TOPOUT,">$opt_o") || die "FATAL: Unable to create output file $opt_o!\n";
while(<TOPIN>){
    $line=$_;
    if ($_=~/^\[ +atoms +\]/){$readatoms=1}
    elsif ($_=~/^\[/){$readatoms=0}
    if ($readatoms && ($_=~/^ *[0-9]+ +/)){
	$_=~s/^ +//;
       ($anum,$atype,$resnum,$resname,$aname,$cgnr,$charge,$mass,$typeB,$chargeB,$massB)=split(/ +/,$_);
        $ATOMNUM[$resnum]{$aname}=$anum;
	if ($aname =~ /[0-9]{2}$/){
	    $altname=$aname;
	    $altname=~s/^([A-Z0-9]+)([0-9])$/$2$1/;
	    $ATOMNUM[$resnum]{$altname}=$anum;
	}
	$foundres{$resnum}=1;
	if (($RESNAME{$resnum} =~ /[A-Z]/) && ($resname ne $RESNAME{$resnum})){
            die "FATAL: Residue $resnum in topology file is a $resname but it was $RESNAME{$resnum} in the disre file.\n"; 
	}
    }

    # Adding restraints HERE
    if ($_ =~ /; Include Position restraint file/){

        # NOE distance restraints
        if ($drnum > 0){
           print TOPOUT "[ distance_restraints ]\n";
           print TOPOUT "; ai   aj  type index  type'  low   up1  up2   fac\n";
           # from file
           for ($dn=0; $dn <= $drnum; $dn++){
	       foreach $pair (split/\,/,$PAIRSFORDR[$dn]){
		   #print STDERR "$pair\n";
                   ($rnum1,$atom1,$rnum2,$atom2)=split(/[\_\|]/,$pair);
		   if ($ATOMNUM[$rnum1]{$atom1} !~ /[0-9]/){
		       die "Atom $atom1 in residue $rnum1 not found in topology file (restraint $dn)\n";
		   }
		   elsif ($ATOMNUM[$rnum2]{$atom2} !~ /[0-9]/){
		       die "Atom $atom2 in residue $rnum2 not found in topology file (restraint $dn)\n";
		   }
		   else{
		       # Note that the division by 10 is necessary to convert Angstroms (X-PLOR) to nm (GROMACS)!
		       printf TOPOUT (" %4d %4d    1 %4d     1  %5.2f %5.2f %5.2f  1.0 ; %4d %3s %4s %4d %3s %4s\n",$ATOMNUM[$rnum1]{$atom1},$ATOMNUM[$rnum2]{$atom2},$dn,0.18,$DIST[$dn]/10,$opt_U, $rnum1,$RESNAME{$rnum1},$atom1,$rnum2,$RESNAME{$rnum2},$atom2);
		       }
		   }
	       }
	}#if drnum

    }#_at inserting position in top file

    print TOPOUT "$line";

}#_while TOPIN

close(TOPIN);
close(TOPOUT);

