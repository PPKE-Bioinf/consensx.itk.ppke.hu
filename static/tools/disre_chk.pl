#!/usr/bin/perl -I/home/szpari/programs/Perl

$|=1;

print STDERR "
                    :-] G R O M A C S [-:
                        legacy tools

                     :-] disrechk.pl  [-:

                   version as of 11 Feb 2016

Option     Filename  Type          Description
------------------------------------------------------------
";


use Getopt::Std;
use Vector3D;

$opt_f="eiwit.pdb";
$opt_d="noe.disre";
$opt_o="disre_chk.out";
$opt_U="unviolated.disre";
$opt_V="violated.disre";
$opt_a=-6;
$opt_t=0;
$opt_e=-6;

getopts('f:d:o:U:V:a:e:t:rvh') || die "Please check your options\n";

if ($opt_h){$HO="yes"} else{$HO="no"}
if ($opt_r){$RO="yes"} else{$RO="no"}
if($opt_v){$V=1;$VO="yes"}else{$VO="no"}

printf STDERR ("  -f  %12s   Input         generic structure file (multi-model PDB at the moment\n",$opt_f);
printf STDERR ("  -d  %12s   Input         disre file\n",$opt_d);
printf STDERR ("  -o  %12s   Output        NOE check report file\n",$opt_o);
printf STDERR ("  -U  %12s   Output        disre file with unviolated restraints\n",$opt_U);
printf STDERR ("  -V  %12s   Output        disre file with violated restraints\n",$opt_V);

print STDERR"
      Option   Type  Value  Description
------------------------------------------------------
";

printf STDERR ("          -h   bool %9s  Print help info and quit\n",$HO);
printf STDERR ("          -v   bool %9s  Be loud and noisy\n",$VO);
printf STDERR ("          -a   int  %9s  Exponent for ambiguous restraint averaging\n",$opt_a);
printf STDERR ("          -e   int  %9s  Exponent for ensemble averaging\n",$opt_e);
printf STDERR ("          -t   float  %5.2f  Tolerance for violation\n",$opt_t);
printf STDERR ("          -r   bool %9s  Rename hydrogens like 2HB -> HB2 in PDB file\n\n",$RO);


if ($opt_h) {
die "

Trying to do meaningful NOE analysis on a bunch of structures. 
No warranty.

";
}

$tolerance=$opt_t;
open (TF,"$opt_d") || die "Cannot open disre file $opt_d for reading\n";
$readdisre=0;
$disrenum=-1;
$pairnum=0;
$previndex=-1;
$pdbok=1;
$violated=0;
while (<TF>){
    chomp;
    next if ($_=~/^ *\#/);
    if ($_=~/^ *[0-9]+/){

       ($disre,$comment)=split(/\#/,$_);
       #if ($comment !~ /[A-Z]/){$pdbok=0}
       $disre=~s/^ +//;
       ($index,$rnum1,$rname1,$atomname1,$rnum2,$rname2,$atomname2,$up1)=split(/ +/,$disre);
       #($ai,$aj,$type,$index,$typea,$low,$up1,$up2,$fac)=split(/ +/,$disre);
       #($atom1,$atom2)=split(/\:/,$comment);
       #$atom1=~s/^ *([0-9]+) *([A-Z][A-Z0-9]*) *$/$1\_$2/;
       #$atom2=~s/^ *([0-9]+) *([A-Z][A-Z0-9]*) *$/$1\_$2/;
       #print "$_ :: $comment :: $atom1 <> $atom2\n" if $V;
       $atom1=$rnum1."_".$atomname1;
       $atom2=$rnum2."_".$atomname2;

       # Storing atoms to read
       $atomtoread{$atom1}=1;
       $atomtoread{$atom2}=1;

       $PAIR[$pairnum]=$atom1."|".$atom2;

       # Keeping intrenal disre numbering with a reference to the original one
       if ($index != $previndex){$disrenum++;$INDEX[$disrenum]=$index} 
       $PAIRSFORDR[$disrenum].=":".$pairnum;

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

       $DIST[$disrenum]=$up1;

       print "$_ | $index DR $disrenum PR $pairnum PFDR $PAIRSFORDR[$disrenum]\n" if $V;
       $pairnum++;
       $previndex=$index;
    }

}#_while TF
   
close(TF);
# numbering starts with 0...
$dno=$disrenum+1;
print STDERR "Read $pairnum atom pairs in $dno distance restraints\n";
if (!$pdbok){
    print STDERR "Warning: one or more restraint lines contain invalid atom names, treating PDB file as input may fail\n";
}

# Guessing input structure file format from extension
if ($opt_f =~ /\.pdb$|\.ent$/){$SFORMAT="pdb"}
else {
    die "Unrecognized struture file format, currently only PDB files (.pdb, .ent) are accepted\n";
}


open (PDB,"$opt_f") || die "Cannot open PDB file $opt_f\n";
$MODEL=0; 
while (<PDB>){
  chop;
  if ($_ =~ /^MODEL/){$MODEL++}
  elsif ($_=~/^ATOM/) {

   $resnum=substr($_,22,5);$resnum=~s/ //g;
   $name=substr($_,12,4);$name=~s/ //g;

   #if ($name eq "H"){$name="HN"}
   $name=~s/([0-9])H([A-Z0-9]+)/H$2$1/ if $opt_r;

   $atomexp=$resnum."_".$name;
   #print "$_ :: $atomexp\n";
 
   if ($atomtoread{$atomexp}) {
       $coordx=substr($_,30,8);$coordx=~s/ //g;
       $coordy=substr($_,38,8);$coordy=~s/ //g;
       $coordz=substr($_,46,8);$coordz=~s/ //g;
       $atomv{$MODEL}{$atomexp}=new Vector3D($coordx,$coordy,$coordz);
       $foundatom{$atomexp}=1;
       #print "  \#MODEL $MODEL atomexp $atomexp\#";$atomv{$MODEL}{$atomexp}->PrintInfo();
   }#_if
 }#_if ATOM
}#_while PDB
close(PDB);
# And now, more music

open (OF,">$opt_o") || die "Cannot create output file $opt_o\n";
open (UF,">$opt_U") || die "Cannot create output file $opt_U\n";
open (VF,">$opt_V") || die "Cannot create output file $opt_V\n";
$unvioln=0;
for ($drn=0; $drn <= $disrenum; $drn++){
    $pairs=$PAIRSFORDR[$drn];
    $pairs=~s/^://;
    @pr=split(/\:/,$pairs);
    $pn=@pr;

    $ensembleavedist=0;
    @AAVE=();

    for ($M=1; $M <= $MODEL; $M++){
	$ambigavedist=0;
	printf OF ("# DR %4d MODEL %3d",$drn,$M) if ($pn > 1);
	foreach $pair (@pr){
	($atomid1,$atomid2)=split(/\|/,$PAIR[$pair]);
	printf OF (" | PAIR $pair $atomid1 $atomid2") if ($pn > 1);
            if ($M == 1){
	     #print STDERR "  \#MODEL $M atom1 $atomid1\#";$atomv{$M}{$atomid1}->PrintInfo() if $opt_v;
	     #print STDERR "  \#MODEL $M atom2 $atomid2\#";$atomv{$M}{$atomid2}->PrintInfo() if $opt_v;
            }
            if ($foundatom{$atomid1}*$foundatom{$atomid2}){
	     $dist=($atomv{$M}{$atomid1}->Diff($atomv{$M}{$atomid2}))->Length();
	     $ambigavedist+=exp($opt_a*log($dist));
             printf OF (" %5.3f",$dist) if ($pn > 1); 
           }
           if (!$foundatom{$atomid1}){print STDERR "ATOM $atomid1 not found\n"}
           if (!$foundatom{$atomid2}){print STDERR "ATOM $atomid2 not found\n"}
	}
	$ambigavedist=exp((1/$opt_a)*log($ambigavedist/$pn));
	$AAVE[$M]=$ambigavedist;
        printf OF (" aaved=%5.3f\n",$ambigavedist) if ($pn > 1); 
	$ensembleavedist+=exp($opt_e*log($ambigavedist));
    }
    $ensembleavedist=exp((1/$opt_e)*log($ensembleavedist/$MODEL));
    printf OF ("%5d %5.2f  <> ",$drn,$DIST[$drn]);
    for ($M=1; $M <= $MODEL; $M++){
        printf OF (" %5.3f",$AAVE[$M]);
    }
    #printf OF ("   %5.3f\n",$ensembleavedist);

    $violation=0;
    if ($ensembleavedist > ($DIST[$drn]+$tolerance)){
	$violated++;
	#$violation=$ensembleavedist-($DIST[$drn]+$tolerance);
	$violation=$ensembleavedist-($DIST[$drn]);
        #print STDERR "violation: $violation\n";
	#$aveviol+=$ensembleavedist-$DIST[$drn];
	$aveviol+=$violation;
	#if (($ensembleavedist-$DIST[$drn]) > $maxviol){
	#    $maxviol=$ensembleavedist-$DIST[$drn];
	#    $maxvioldr=$drn;
	#}
	if ($violation > $maxviol){
	    $maxviol=$violation;
	    $maxvioldr=$drn;
	}
        # writing violated list
	printf VF ("#Was $INDEX[$drn] in original list | ACTUAL %7.3f \n",$ensembleavedist);
	foreach $pair (@pr){
	  ($atomid1,$atomid2)=split(/\|/,$PAIR[$pair]);
          ($rnum1,$atom1)=split(/\_/,$atomid1);
          ($rnum2,$atom2)=split(/\_/,$atomid2);
          # ensuring for output disre list to start with 1 (was violated-1 in the original code)
	  printf VF (" %4d   %3d %3s %4s   %3d %3s %4s   %5.2f\n",$violated,$rnum1,$RESNAME{$rnum1},$atom1,$rnum2,$RESNAME{$rnum2},$atom2,$DIST[$drn]);
	}
    }
    # writing unviolated list
    else{
	printf UF ("#Was $INDEX[$drn] in original list | ACTUAL %7.3f \n",$ensembleavedist);
	foreach $pair (@pr){
	  ($atomid1,$atomid2)=split(/\|/,$PAIR[$pair]);
          ($rnum1,$atom1)=split(/\_/,$atomid1);
          ($rnum2,$atom2)=split(/\_/,$atomid2);
          # +1 added for output disre list to start with 1
	  printf UF (" %4d   %3d %3s %4s   %3d %3s %4s   %5.2f\n",$unvioln+1,$rnum1,$RESNAME{$rnum1},$atom1,$rnum2,$RESNAME{$rnum2},$atom2,$DIST[$drn]);
	}
	$unvioln++;
    }
    printf OF (" EAD %5.3f VIOL %5.3f\n",$ensembleavedist,$violation);

}
if ($violated > 0){
    $aveviol/=$violated;
}
printf OF ("# There are %d violated restraints. Average violation is %5.3f A, maximum violation is %5.3f A",$violated,$aveviol,$maxviol);
if ($maxviol){print OF " (restraint $INDEX[$maxvioldr])"}
print OF "\n";

printf STDERR ("# There are %d violated restraints. Average violation is %5.3f A, maximum violation is %5.3f A",$violated,$aveviol,$maxviol);
if ($maxviol){print STDERR " (restraint $INDEX[$maxvioldr])"}
print STDERR "\n";

close(OF);
close(UF);
close(VF);
