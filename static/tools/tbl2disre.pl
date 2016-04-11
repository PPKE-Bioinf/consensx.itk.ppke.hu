#!/usr/bin/perl

print STDERR "
############################################################
# tbl2disre.pl -p<pdb_structure> < tbl_file > disre_file   #
#                                                          #
# An attempt to convert a tbl file into a format where all #
# atoms involved are explicitly named, not unlike          #
# GROMACS's distance restraint section in the topology     #
# file.                                                    #
# The pdb file is used the identify all existing atoms     #
# that might match the restrained ones.                    #
############################################################

";

$|=1;
use Getopt::Std;



getopts('p:hv') || die "ERROR in options, exiting\n";


# Read distance restraint file (tbl format)
$dn=0;
while(<>){
    chomp;

    # omitting everything after comments
    $_=~s/\!.*$//;

    if ($_ =~ /^ *$/){$getres=0}
    # restraint line
    if ($_=~/^ *ass/i){
	$dn++;
	$getres=1;
	$RESLINE[$dn-1]=" # ".$_;
    }
    elsif($getres){
	$RESLINE[$dn-1].=" # ".$_;
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
";

 foreach $notfound (sort {$rnfn{$a} <=> $rnfn{$b}} keys %NOTFOUND){
    ($rn,$an)=split(/\_/,$notfound);
    print STDERR "$RESNUM{$notfound} $rn $RESNAME{$rn} $an\n";
 }
}#_if NFN

exit;

#---------------================================----------#


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

