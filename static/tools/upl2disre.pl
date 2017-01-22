#!/usr/bin/perl

print STDERR "
############################################################
# upl2disre.pl -p<pdb_structure> < upl_file > disre_file   #
#                                                          #
# An attempt to convert an upl file into a format where    #
# all atoms involved are explicitly named, not unlike      #
# GROMACS's distance restraint section in the topology     #
# file.                                                    #
# The pdb file is used the identify all existing atoms     #
# that might match the restrained ones.                    #
############################################################

";

$|=1;
use Getopt::Std;



getopts('p:hv') || die "ERROR in options, exiting\n";
exit if $opt_h;

# Read distance restraint file (upl format)
$dn=0;
while(<>){
    chomp;
    $_=~s/^ +//;
    # Some trivial conversions before doing any checks...
    $_=~s/CYSS/CYS /g;
    # distance restraint line
    if ($_=~/^[0-9]+ +[A-Z]{3} +[A-Z0-9]+ +[0-9]+ +[A-Z]{3} +[A-Z0-9]+ +[0-9\.]+/){
	($resnum1[$dn],$resname1[$dn],$atomname1[$dn],$resnum2[$dn],$resname2[$dn],$atomname2[$dn],$distance[$dn])=split(/ +/,$_);
	$upline[$dn]=$_;
	$dn++;
    }
}
print STDERR "Read $dn distance restraints\n";

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
	print "Stored: \#$atomid\# :: $ATOMNUM{$atomid} ($RESNAME{$resnum})\n" if $opt_v;
    }
}


# Iterating over restraints and parsing atoms
for ($n=0; $n < $dn; $n++){

    $found_atom1=0;
    $found_atom2=0;
    $ERROR="";
    # residue name check
    if ($resname1[$n] ne $RESNAME{$resnum1[$n]}){$ERROR.="Residue name mismatch: residue $resnum1[$n]: $resname1[$n] (restraint) vs. $RESNAME{$resnum1[$n]} (pdb)\n"}
    if ($resname2[$n] ne $RESNAME{$resnum2[$n]}){$ERROR.="Residue name mismatch: residue $resnum1[$n]: $resname2[$n] (restraint) vs. $RESNAME{$resnum2[$n]} (pdb)\n"}
    die($ERROR) unless ($ERROR eq "");

    $UNRESOLVED="";
    @A1=resolveatoms($atomname1[$n],$resnum1[$n]);
    @A2=resolveatoms($atomname2[$n],$resnum2[$n]);

print "# $upline[$n]\n"; 
    foreach $an1 (@A1){
	foreach $an2 (@A2){
	    $an1o=$an1;$an1o=~s/^.*\_//;
	    $an2o=$an2;$an2o=~s/^.*\_//;
	    #print "$n $resnum1[$n] $resname1[$n] $an1o $resnum2[$n] $resname2[$n] $an2o $distance[$n]\n";
	    printf (" %4d   %3d %3s %4s   %3d %3s %4s   %5.2f $UNRESOLVED\n",$n,$resnum1[$n],$resname1[$n],$an1o,$resnum2[$n],$resname2[$n],$an2o,$distance[$n]);
	}
    }

}

# Reporting atoms that could not be converted
if ($NFN){
 print STDERR "
List of $NFN unresolved atom(s) found
=====================================
";

 foreach $notfound (sort {$rnfn{$a} <=> $rnfn{$b}} keys %NOTFOUND){
    ($rn,$an)=split(/\_/,$notfound);
    print STDERR "$RESNUM{$notfound} $rn $RESNAME{$rn} $an\n";
 }
}#_if NFN

exit;

#---------------================================----------#

sub resolveatoms{

    @ALIST=();
    $atomname=shift;
    $resnum=shift;

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
