#!/usr/bin/perl

$|=1;
use Getopt::Std;

$opt_p="eiwit.pdb";

getopts('p:l:hv') || die "ERROR in options, exiting\n";

print STDERR "
############################################################
# str2disre.pl -p<pdb_structure> [-l<convert.log>] [-v]    #
#                                < str_file > disre_file   #
#                                                          #
# An attempt to convert an str file (distance restraints   #
# v2 as denoted in PDB) into a format where all            #
# atoms involved are explicitly named, not unlike          #
# GROMACS's distance restraint section in the topology     #
# file.                                                    #
# The pdb file is used the identify all existing atoms     #
# that might match the restrained ones.                    #
# Problems are reported in the log file.                   #
############################################################

" if ($opt_v || $opt_h);


exit if $opt_h;
if (! -e $opt_p){die "Input pdb file $opt_p does not exist, exiting.\n"}
$opt_l=$opt_p;
$opt_l=~s/.pdb|.ent/_str2disre.log/;
open (LF,">$opt_l") || die "Cannot create log file $opt_l\n";

# Read distance restraint file (str format)
$rn=0;$allrn=0;
$readrst=0;
$readloop=0;
@loopheader=();
$rlist=0;
$my_rid=0; # used for all restraints processed as they might come from multiple lists
while(<>){
    chomp;

    if (($readres) && ($_ =~ /loop_/) && ($rn == 0)){$readloop=1}
    if ($_=~/save_ *$/){$readres=0;$rn=0;$PREVRID=0;}
    if ($_=~/save.*istanc.*traint/){$readres=1;$rlist++}
    if ($_=~/stop_/){$readloop=0;@loopheader=()}
    #print ">$_ | $readres $readloop\n";
    if (($readloop) && ($_=~/^ +\_[A-Z]/i)){
	#print "LOOP header:> $_\n";
	$_=~s/ //g;
	push(@loopheader,$_);
    }
    # Hopefully this is a restraint line
    #elsif (($readloop) && ($_=~/[0-9]\.[0-9].*[0-9]*H.*[0-9]*H/) && ($_ !~ / [NO] /)){
    #elsif (($readloop) && ($_=~/[0-9]\.[0-9]/) && ($_ !~ / [NO] /)){
    elsif (($readloop) && ($_=~/[0-9]\.[0-9]/)){
	if ($rn==0){
	    &processheader;
	}
	$_=~s/^ +//;
	@RL=split(/ +/,$_);

        if (($ATP1_C ne "") && ($ATP2_C ne "") && (($RL[$ATP1_C] ne "H") || ($RL[$ATP2_C] ne "H"))){
           print STDERR "Nem NOE: $_\n";
           next;
        }

	$RID=$RL[$RSID_C];
	if ($RID > $PREVRID){
           #print "DEBUG> $RID > $PREVRID\n";
           $rn++;$allrn++;$my_rid++;
	   $ORIG_RID[$my_rid]=$RID;
	}
	#print "DEBUG> $_\nRID=$RID PREVRID=$PREVRID my_rid=$my_rid >> $RL[$RNO1_C] $RL[$RNM1_C] $RL[$ATM1_C] $RL[$RNO2_C] $RL[$RNM2_C] $RL[$ATM2_C] $RL[$DIST_C] \n\n";
	if ($RESTRAINT[$my_rid] =~ /[A-Z0-9]/){
	    $RESTRAINT[$my_rid].="\n"."$RL[$RNO1_C] $RL[$RNM1_C] $RL[$ATM1_C] $RL[$RNO2_C] $RL[$RNM2_C] $RL[$ATM2_C] $RL[$DIST_C]";
	    #$RESTRAINT[$my_rid].=sprintf("\n  %3d %3s %4s   %3d %3s %4s %5.2f",$RL[$RNO1_C],$RL[$RNM1_C],$RL[$ATM1_C],$RL[$RNO2_C],$RL[$RNM2_C],$RL[$ATM2_C],$RL[$DIST_C]);
	}
	else{
	    $RESTRAINT[$my_rid]="$RL[$RNO1_C] $RL[$RNM1_C] $RL[$ATM1_C] $RL[$RNO2_C] $RL[$RNM2_C] $RL[$ATM2_C] $RL[$DIST_C]";
	    #$RESTRAINT[$my_rid]=sprintf("  %3d %3s %4s   %3d %3s %4s %5.2f",$RL[$RNO1_C],$RL[$RNM1_C],$RL[$ATM1_C],$RL[$RNO2_C],$RL[$RNM2_C],$RL[$ATM2_C],$RL[$DIST_C]);
	}

	if ($RESTRAINT[$my_rid] !~ /^ *[0-9]+ +[A-Z]{3} +[0-9A-Z]+ *[0-9]+ +[A-Z]{3} +[0-9A-Z]+ +[0-9]\.[0-9]/){
	   print LF "WARNING: wrong format? >$_\nRID=$RID my_rid =$my_rid \n$RESTRAINT[$my_rid]\n\n";
	   print STDERR "WARNING: wrong format? >$_\nRID=$RID my_rid =$my_rid \n$RESTRAINT[$my_rid]\n\n" if $opt_v;
	}

	$PREVRID=$RID;
    }


}



#exit;

# reading PDB file
open (PF,"$opt_p") || die "Cannot open pdb file $opt_p\n";

while(<PF>){
    chomp;
    # END or ENDMDL: finish reading 
    # (only first model is read from multi-model files)
    if ($_=~/^END/){last;}
    elsif ($_=~/^ATOM/){
	$atomnum=substr($_,6,5);$resnum=~s/ //g;
	$resname=substr($_,17,3);
	$resnum=substr($_,22,5);$resnum=~s/ //g;
	$atomname=substr($_,12,4);$atomname=~s/ //g;
	$foundres{$resnum}=1;
	$atomid=$resnum."_".$atomname;
	$ATOMNUM{$atomid}=$atomnum;
	$RESNAME{$resnum}=$resname;
	print STDERR "Stored: \#$atomid\# :: $ATOMNUM{$atomid} ($RESNAME{$resnum})\n" if $opt_v;
    }
}


print LF "Read $allrn distance restraints from $rlist list(s)\n";
print STDERR "Read $allrn distance restraints from $rlist list(s)\n" if $opt_v;
for ($rid=1; $rid<= $allrn; $rid++) {
    @A1=();
    @A2=();
    foreach $rst (split(/\n/,$RESTRAINT[$rid])){
	printf ("#%4d $rst\n",$ORIG_RID[$rid]);
	$rst=~s/^ +//;
        ($resnum1,$restype1,$atom1,$resnum2,$restype2,$atom2,$dist)=split(/ +/,$rst);
	@A1=resolveatoms($atom1,$resnum1,$restype1);
	@A2=resolveatoms($atom2,$resnum2,$restype2);
	foreach $aid1 (@A1){
	    $a1=$aid1;$a1=~s/^.*\_//;
	    foreach $aid2 (@A2){
		$a2=$aid2;$a2=~s/^.*\_//;
		printf(" %4d   %3d %3s %4s   %3d %3s %4s %5.2f\n",$rid,$resnum1,$restype1,$a1,$resnum2,$restype2,$a2,$dist);
	    }
	}
    }
}


# shall simply split the restraint lines and for the same restraint
# call resolveatom for atom1 and atom2 for each line and sum them up...



# Iterating over restraints and parsing atoms
for ($n=0; $n < $dn; $n++){

    @A1=();
    @A2=();
    @LINESWITHOUTDIST=();

    &processresline($RESLINE[$n]);

    $maxdist=sprintf("%5.2f",$dist[0]+$dist[2]);

    %WRITTEN=();
    print "$rlo\n" if (@LINESWITHOUTDIST);
    foreach $lwd (@LINESWITHOUTDIST){
	$lwd=~s/ DIST /$maxdist/;
	print $lwd unless $WRITTEN{$lwd};
	$WRITTEN{$lwd}=1;
    }

}

# Reporting atoms that could not be converted
if ($NFN){
 print STDERR "
List of $NFN unresolved atom(s) found
=====================================
" if $opt_v;
 print LF "
List of $NFN unresolved atom(s) found
=====================================
";

 foreach $notfound (sort {$rnfn{$a} <=> $rnfn{$b}} keys %NOTFOUND){
    ($rn,$an)=split(/\_/,$notfound);
    print STDERR "$RESNUM{$notfound} $rn $RESNAME{$rn} $an\n" if $opt_v;
    print LF "$RESNUM{$notfound} $rn $RESNAME{$rn} $an\n";
 }
}#_if NFN

exit;

#---------------================================----------#


sub processheader{

    $n=0;
    foreach $hdk (@loopheader){

	#print "hdk: $hdk | $n\n";
# _Gen_dist_constraint.Distance_upper_bound_val
	if ($hdk=~/istan.*upper.*_val/){$DIST_C=$n}
# _Gen_dist_constraint.PDB_residue_name_1
	elsif ($hdk=~/esid.*_na.*1/){$RNM1_C=$n}
	elsif ($hdk=~/esid.*_na.*2/){$RNM2_C=$n}
# _Gen_dist_constraint.PDB_residue_no_1
	elsif ($hdk=~/esid.*_n.*1/){$RNO1_C=$n}
	elsif ($hdk=~/esid.*_n.*2/){$RNO2_C=$n}
# _Gen_dist_constraint.PDB_atom_name_1
	elsif ($hdk=~/PDB.*tom_na.*1/){$ATM1_C=$n}
	elsif ($hdk=~/PDB.*tom_na.*2/){$ATM2_C=$n}
# _Gen_dist_constraint.ID 
	elsif ($hdk=~/raint\.ID/){$RSID_C=$n}
# _Gen_dist_constraint.Atom_type_1
	elsif ($hdk=~/raint\.Atom_type_1/){$ATP1_C=$n}
	elsif ($hdk=~/raint\.Atom_type_2/){$ATP2_C=$n}
	$n++;
    }
    #print "DIST_C $DIST_C RNM1 $RNM1_C RNO_1 $RNO1_C ATM_1 $ATM1_C RSID $RSID_C\n";


}

sub processresline{
    $rl=shift;
    $rlo=$rl;
    $rlo=~s/ \# /\n\#/g;$rlo=~s/^\n//;
    # Adding spaces to help splitting
    $rl=~s/\(/ ( /g;
    $rl=~s/\)/ ) /g;
    $rl=~s/\" +\"/"VOID"/g;

    @RP=split(/ +/,$rl);
    $PAREN=0;
    $ATOMBLOCK=0;
    $LASTKEY="";
    $DISTNUM=0;
    %a1found=();%a2found=();


    foreach $rk (@RP){
	next if ($rk eq "#");
	if ($rk eq "("){
	    if ($PAREN==0){$ATOMBLOCK++}
	    $PAREN++;
	}
	elsif ($rk eq ")"){
	    $PAREN--;
	}
	elsif ($rk =~ /^or$/i){
	    if ($PAREN == 0){
		#print "paren 0 record lines\n";
		&recordlineswithoutdist;
		#@LINESWITHOUTDIST=();
		$ATOMBLOCK=0;
		@A1=();@A2=();
	    }
	}
	elsif ($rk =~ /^seg/i){
	    $LASTKEY="SEGID";
	}
	elsif ($rk =~ /^res/i){
	    $LASTKEY="RESID";
	}
	elsif ($rk =~ /^nam/i){
	    $LASTKEY="NAME";
	}
	elsif ($rk =~ /^weigh/i){
	    $LASTKEY="WEIGHT";
	}
	if ($rk =~ /^[0-9]+$/){
	    if ($LASTKEY eq "RESID"){
		$resnum=$rk;
	    }
	}
	if ($rk =~ /^[A-Z0-9]*H[A-Z0-9\#\*\%\+]*$/){
	    if ($LASTKEY eq "NAME"){
		$atomname=$rk;
		$UNRESOLVED="";
		#rint "$n -> $resnum $atomname\n";
		if ($ATOMBLOCK == 1){
		    @TA1=resolveatoms($atomname,$resnum);
		    push(@A1,@TA1);
		}
		elsif ($ATOMBLOCK == 2){
		    @TA2=resolveatoms($atomname,$resnum);
		    push(@A2,@TA2);
		}
	    }
	}
	if (($rk =~ /^[0-9\.]+$/) && ($PAREN == 0)){
	    if ($LASTKEY eq "WEIGHT"){
		$weight=$rk;
	    }
	    else{
		$dist[$DISTNUM]=$rk;
		$DISTNUM++;
	    }
	}

    }#_foreach $rk

    #print "A1 $#A1 A2 $#A2\n";
    if ((@A1 > 0) && (@A2 > 0)){
	&recordlineswithoutdist;
    }

 

}#_sub processresline

sub resolveatoms{

    @ALIST=();
    $atomname=shift;
    $resnum=shift;
    $restype=shift;

    # Handling Ile MG/MD and similar 'lazy' notations...
    if ($restype eq "ILE"){
	if ($atomname=~/^[QM]G$/){$atomname.="2"}
	elsif ($atomname=~/^[QM]D$/){$atomname.="1"}
    }
    elsif ($restype eq "THR"){
	if ($atomname=~/^[QM]G$/){$atomname.="2"}
    }

    # pseudoatom handling
    if ($atomname eq "QR"){
	foreach $a ("HD1","HD2","HE1","HE2","HZ"){
	    $aa=$a;$aa=~s/(.*)([12])$/$2$1/;
	    $aid=$resnum."_$a";
	    $aida=$resnum."_$aa";
	    if ($ATOMNUM{$aid} =~ /[0-9]/){
		push(@ALIST,$aid);
		print "Found: $aid : $ATOMNUM{$aid} \n" if $opt_v;
	    }
            # preventing HZ from getting duplicate :-)
	    elsif(($aida ne $aid) && ($ATOMNUM{$aida} =~ /[0-9]/)){
		push(@ALIST,$aida);
		print "Found: $aida :  $ATOMNUM{$aida}\n" if $opt_v;
	    }
	}
    }
    elsif ($atomname =~ /^QQ/){
	$anr=$atomname;$anr=~s/QQ/H/;
	foreach $a (11,12,13,21,22,23){
            # names like HD21 etc.
	    $aid1a=$resnum."_".$anr.$a;
            # if not found, try names like 1HD2
            ($a1,$a2)=split(//,$a);
	    $aid1b=$resnum."_".$a2.$anr.$a1;
	    print "aid1a \#$aid1a\# aid1b \#$aid1b\#\n" if $opt_v;
	    if ($ATOMNUM{$aid1a} =~ /[0-9]/){
		push(@ALIST,$aid1a);
		print "Found: $aid1a : $ATOMNUM{$aid1a} \n" if $opt_v;
	    }
	    elsif($ATOMNUM{$aid1b} =~ /[0-9]/){
		push(@ALIST,$aid1b);
		print "Found: $aid1b :  $ATOMNUM{$aid1b}\n" if $opt_v;
	    }
	}#_foreach possible atom name endings
    }# QQ	
    elsif ($atomname =~ /^[QM]/){
	$anr=$atomname;$anr=~s/[QM]/H/;
	foreach $a (1,2,3){
	    $aid1a=$resnum."_".$anr.$a;
	    $aid1b=$resnum."_".$a.$anr;
	    print "aid1a \#$aid1a\# aid1b \#$aid1b\#\n" if $opt_v;
	    if ($ATOMNUM{$aid1a} =~ /[0-9]/){
		push(@ALIST,$aid1a);
		print "Found: $aid1a : $ATOMNUM{$aid1a} \n" if $opt_v;
	    }
	    elsif($ATOMNUM{$aid1b} =~ /[0-9]/){
		push(@ALIST,$aid1b);
		print "Found: $aid1b :  $ATOMNUM{$aid1b}\n" if $opt_v;
	    }
	}
    }
    else{
	$aid=$resnum."_".$atomname;
	if ($ATOMNUM{$aid}){
	    @ALIST=$aid;
	}
        # A common problem...
	elsif($atomname eq "HN"){
	    $aida=$resnum."_H";
	    if ($ATOMNUM{$aida}){
		@ALIST=$aida;
                $UNRESOLVED=" # $resnum HN -> H"; 
	    }
	}
    }
    if (@ALIST == 0){
	#print STDERR "Could not resolve atom $resnum $atomname\n";
	$UNRESOLVED.=" # $resnum $atomname NOT FOUND";
	$aid=$resnum."_".$atomname;
	push(@ALIST,$aid);
	if ($NOTFOUND{$aid} == 0){
	    $NOTFOUND{$aid}=1;
	    $NFN++;
	    $rnfn{$aid}=$resnum;
	}
    }
    return @ALIST;

}#sub resolveatoms


sub resolveatoms_tbl2disre{

    @ALIST=();
    my $atomname=shift;
    my $resnum=shift;

    if ($atomname =~ /[\#\+\*\%]$/){
	$anr=$atomname;$anr=~s/[\#\+\*\%]//;
	#print "$atomname => $anr\n";
        # all possible options in one cycle
	foreach $a (1,2,3,11,12,13,21,22,23){
	    $aid1a=$resnum."_".$anr.$a;
	    $aid1b=$anr.$a;
	    $aid1b=~s/^(.+)([0-9])$/$2$1/;
	    $aid1b=$resnum."_".$aid1b;
	    print "aid1a \#$aid1a\# aid1b \#$aid1b\#\n" if $opt_v;
	    if ($ATOMNUM{$aid1a} =~ /[0-9]/){
		push(@ALIST,$aid1a);
		print "Found: $aid1a : $ATOMNUM{$aid1a} \n" if $opt_v;
	    }
	    elsif($ATOMNUM{$aid1b} =~ /[0-9]/){
	    	push(@ALIST,$aid1b);
	        print "Found: $aid1b :  $ATOMNUM{$aid1b}\n" if $opt_v;
	    }
	}
    }
    else{
	$aid=$resnum."_".$atomname;
	if ($ATOMNUM{$aid}){
	    @ALIST=$aid;
	}
        # A common problem...
	elsif($atomname eq "HN"){
	    $aida=$resnum."_H";
	    if ($ATOMNUM{$aida}){
		@ALIST=$aida;
                $UNRESOLVED=" # $resnum HN -> H"; 
	    }
	}
    }
    if (@ALIST == 0){
	#print STDERR "Could not resolve atom $resnum $atomname\n";
	$UNRESOLVED.=" # $resnum $atomname NOT FOUND";
	$aid=$resnum."_".$atomname;
	push(@ALIST,$aid);
	if ($NOTFOUND{$aid} == 0){
	    $NOTFOUND{$aid}=1;
	    $NFN++;
	    $rnfn{$aid}=$resnum;
	}
    }
    return @ALIST;

}#sub resolveatoms


sub recordlineswithoutdist{

    foreach $an1 (@A1){
        foreach $an2 (@A2){

            ($rn1o,$an1o)=split(/\_/,$an1);
            ($rn2o,$an2o)=split(/\_/,$an2);
	    #print "$an1, $an2\n";
            #print "$n $resnum1[$n] $resname1[$n] $an1o $resnum2[$n] $resname2[$n] $an2o $distance[$n]\n";
	    $lwd=sprintf (" %4d   %3d %3s %4s   %3d %3s %4s   DIST $UNRESOLVED\n",$n,$rn1o,$RESNAME{$rn1o},$an1o,$rn2o,$RESNAME{$rn2o},$an2o);
	    push(@LINESWITHOUTDIST,$lwd);
        }
    }

}

