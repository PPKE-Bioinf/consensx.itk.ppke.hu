#!/usr/bin/perl

print STDERR"
###################################################
# atomconverter.pl -f<from_fromat> -t<to_format>  #
#                  [-q] [-o] [-n] [-H]            #
#                  < in_noe_list > out_noe_list   #
#                                                 #
# Converts atom names from one format to another. #
# Known formats: BMRB, UCSF, XPLOR, MSI, PDB,     #
# SYBYL, MIDAS, DIANA, GROMACS, MOLMOL            #
# Conversion is based on the file:                #
# http://www.bmrb.wisc.edu/ref_info/atom_nom.tbl  #
# GROMACS and MOLMOL formats added later          #
# Option -n causes to rename atoms as 2HB -> HB2  #
# in the PDB file BEFORE attempting format        #
# parsing.                                        #
# -H forces H -> NH conversion BEFORE parsing     #
# -q causes to skip pseudoatoms (default for      #
#    MOLMOL as input format)                      #
# By default, the program will convert terminal   #
# oxygen atoms to PDB format, use option -o to    #
# convert these atoms to the specified format.    #
#                                                 #
###################################################

";

&Init;
use Getopt::Std;
getopts('f:t:oqnHhv') || die "Valid options are: -f<format>, -t<format>, -h[elp], -v[erbose]\n";
if ($opt_h){
    print "Conversion table used:\n$conversion_table\n";
    exit;
}

if (($opt_f !~ /[a-z]/i) || ($opt_t !~ /[a-z]/i)){
    die "Please specify valid 'from' and 'to' formats\n";
}
if ($opt_f !~/^BMRB|UCSF|XPLOR|MSI|PDB|SYBYL|MIDAS|DIANA|GROMACS|MOLMOL$/i){
    die "Invalid 'from' format, please choose from known formats only\n";
}
if ($opt_t !~/^BMRB|UCSF|XPLOR|MSI|PDB|SYBYL|MIDAS|DIANA|GROMACS|MOLMOL$/i){
    die "Invalid 'to' format, please choose from known formats only\n";
}
if ($opt_f eq $opt_t){
    die "'from' and 'to' formats are the same, nothing to do.\n";
}


# converting to uppercase
$opt_f=uc($opt_f);
$opt_t=uc($opt_t);

if ($opt_f eq "MOLMOL"){$opt_q=1}

# reading conversion table
foreach $cl (split(/\n/,$conversion_table)){
    # skip comment or blank lines
    next if (($cl =~ /^\#/) || ($cl !~/[a-z]/i));
    # header line
    if ($cl =~ /^A\tAAA/){
	$n=0;
	foreach $colname (split(/\t/,$cl)){
	    $colname=~s/ //g;
	    if ($colname eq $opt_f){$FROM=$n}
	    if ($colname eq $opt_t){$TO=$n}
	    $n++;
	}
    }
    # atom name line
    else{
	@ac=split(/\t/,$cl);
        if ($ac[$FROM] =~ /[A-Z0-9]/){
	    $from_id="$ac[$FROM] $ac[1]";
	    #print "ILYEN VOLT: \#$from_id\#  ";
	    $from_id=~s/^ +//;$from_id=~s/ +$//;$from_id=~s/  / /g;
	    #print "ILYEN LETT: \#$from_id\# \n";
            $to_id="$ac[$TO]";
	    if (length($to_id) < 4){$to_id=" ".$to_id}
	    while (length($to_id) < 4){$to_id=$to_id." "}
	    $CONVERT{$from_id}=$to_id;
	}
        #print " $from_id -> $to_id\n"; # DEBUG

    }
}

# reading PDB file
while(<>){
    # simplest...
    #$_=~s/ ([A-Z]{3}) +([0-9A-Z]+) / $1 $2 /g;
    if ($_=~/^ATOM/){
     $_=~s/ASP\-/ASP /;
     $_=~s/GLU\-/GLU /;
     $_=~s/ARG\+/ARG /;
     $_=~s/LYS\+/LYS /;
     $_=~s/HIS\+/HIS /;
     $inline=$_;chomp($inline);
     $from_id=substr($inline,12,8);

     if (($opt_H) && ($from_id =~/ H /)){
	 $ofi=$from_id;
	 $from_id=~s/ H / HN/;
	 #$_=~s/$ofi/$from_id/;
     }
     #$from_id=~s/^ +//;$from_id=~s/ +$//; # failsafe way to remove initial and ending spaces
     # Skip MOLMOL pseudoatoms
     if (($opt_q) && ($from_id =~ /Q/)){next}
     #elsif (($opt_f =~/MOLMOL/i) && ($from_id =~ /^[0-9][A-Z]/)){$mfrom_id=$from_id;$mfrom_id=~s/  /   /;print "$from_id > $mfrom_id\n"}

     # convert only if target format specified (leave as it is if not)
     $foundfromid=0;

     #print "ILYEN VOLT: \#$from_id\#  ";
     $from_id=~s/^ +//;$from_id=~s/ +$//;$from_id=~s/  +/ /g;
     #print "ILYEN LETT: \#$from_id\#  \n";

     if ($CONVERT{$from_id} =~ /[0-9A-Z]/){
       $foundfromid=1;
       substr($_,12,4)=$CONVERT{$from_id};
       print ">> CONVERT! $from_id -> $CONVERT{$from_id}\n$w\n$_\n" if $opt_v;
     }
     # trying to throw the first number back, maybe that helps
     elsif($opt_n){
	 $from_id=~s/([0-9])([A-Z0-9]+)/$2$1/;
	 if ($CONVERT{$from_id} =~ /[0-9A-Z]/){
	     $foundfromid=1;
	     substr($_,12,4)=$CONVERT{$from_id};
	     print ">> CONVERT! $from_id -> $CONVERT{$from_id}\n$w\n$_\n" if $opt_v;
	 }
     }
     else{
      $res=$from_id;$res=~s/^.*([A-Z]{3})$/$1/;
      $altfromid=$from_id;$altfromid=~s/[A-Z]{3}$/XXX/;
      #print STDERR "$altfromid $CONVERT{$altfromid}\n";
      if ($CONVERT{$altfromid} =~ /[0-9A-Z]/){
	 $foundfromid=1; 
       $alttoid=$CONVERT{$altfromid};
       # handling terminal atoms: default is PDB convention
       # as consensx will get stuck because shiftx expects O and OXT
       if (!$opt_o){
	   $alttoid=~s/OT1 XXX/O   XXX/;
	   $alttoid=~s/OC1 XXX/O   XXX/;
	   $alttoid=~s/O1  XXX/O   XXX/;
	   $alttoid=~s/O'  XXX/O   XXX/;
	   $alttoid=~s/OT2 XXX/OXT XXX/;
	   $alttoid=~s/OC2 XXX/OXT XXX/;
	   $alttoid=~s/O2  XXX/OXT XXX/;
	   $alttoid=~s/O'' XXX/OXT XXX/;
       }
       $alttoid=~s/[A-Z]{3}$/$res/;
       substr($_,12,4)=$alttoid;
       #$_=~s/ $altfromid / $alttoid /;
       print ">> CONVERT! $altfromid -> $alttoid\n" if $opt_v;
       #print STDERR "$altfromid $CONVERT{$altfromid} $alttoid\n$_\n";
      }
     }
     if (!$foundfromid){
       print STDERR "WARNING: $from_id could not be converted, please check format options. Each unconverted atom type is reported only once.\n" unless $REPORTED{$from_id};
       $REPORTED{$from_id}=1;
     }
    }
    print "$_";
}


sub Init{
# Conversion table put here to avoid dependence on external files.
# This table is a slightly modified version of 
# http://www.bmrb.wisc.edu/ref_info/atom_nom.tbl
$conversion_table="
# Correlation of hydrogen atom naming systems, including diastereotopic
# protons. The original version of this table was created by Charles 
# Hoogstraten.
#
# BMRB  = System in use at BioMagResBank (IUPAC/IUB Biochemistry 9, 3471-3479
#	  [1970]).
#
# SC    = Stereochemical designations 
#
# UCSF  = Mardigras-type software (peptide protonated with newhyd utility).
#
# XPLOR = Peptide protonated with XPLOR 3.1.  Atom nomenclature is derived
#	  from the X-PLOR topology file topallhdg.pro.  
#
# MSI   = Artificial peptide created with InsightII.  For the side chain
#	  protons attached to nitrogen in ASN, GLN, and ARG, the atom
#	  nomenclature does not reflect the potential stereoisomerism of the
#	  planar amide and guanidinium groups.  The correlation with Z and E
#	  nomenclature listed here simply reflects the state of the artificial
#	  peptide generated as an example.  The CG1 and CG2 atoms for VAL in a
#	  peptide generated by InsightII are not labeled according to IUPAC
#	  rules, while the CD1 and CD2 atoms for LEU are.
#
# PDB   = PDB nomenclature (Taken from PDB entry 6I1B REVDAT 15-OCT-92.)
#
# SYBYL = The atom nomenclature was taken from the xxx.res files supplied with
#	  the software package Sybyl version 6.2 from Tripos, Inc.  
#
# MIDAS = MidasPlus from the Computer Graphics Laboratory at UCSF.  The atom 
#	  nomenclature has been taken from the XXX.ins files supplied with the
#	  software.  The prochiral atoms have not been correlated with the
#	  BMRB assignments at this time.  Hydrogens are not included in the
#	  XXX.ins template files.
#
# Note-1: The prochiral methyl group names may reflect convention of code
#	  generating heavy atom names if protons are added later.
#
# Note-2: '*' The stereochemical assignments for the named atoms have not been
#	  determined for these software systems.
#
# Note-3: The Z and E nomenclature is defined in the paper by Blackwood, J.E.,
#	  Gladys, C.L., Loening, K.L., Petrarca, A.E., and Rush, J.E., 
#	  'Unambiguous Specification of Stereoisomerism about a Double Bond,'
#	  J. Amer. Chem. Soc. 90, 509-510 (1968).
#
# Note-3: For the terminal amine and carboxyl atoms, 'X' has been used as a
#	  dummy value for the amino acid type.
#
# Note-4: The terminal secondary amine protons for PRO have been included
#	  with the other PRO atoms.
#
# Note-5: Fields in the table are separated by tabs.
#
# Note-6: Please report errors, updates, or extensions to Eldon Ulrich
#	  (elu@nmrfam.wisc.edu)

# A.A.		BMRB	SC	PDB	UCSF	MSI	XPLOR	SYBYL*	MIDAS*	DIANA
# ----		----	----	----	----	----	-----	-----	-----	-----

# Header line for determining column for each format

A	AAA	BMRB	SC	PDB	UCSF	MSI	XPLOR	SYBYL	MIDAS	DIANA	GROMACS	MOLMOL

# Terminal amine protons

X	XXX	H1		1H	HN1	HN1	HT1	HNCAP			H1	1H
X	XXX	H2		2H	HN2	HN2	HT2				H2	2H
X	XXX	H3		3H	HN3	HN3	HT3				H3	3H

# Terminal carboxyl atoms

X	XXX	H''						HOCAP		
X	XXX	O'		O		O	OT1	O	O		OC1	O1
X	XXX	O''		OXT		OXT	OT2	OXT	OXT		OC2	O2

# Residue atoms

A	ALA	H		H	HN	HN	HN	H		HN	H	H
A	ALA	HA		HA	HA	HA	HA	HA		HA	HA	HA
A	ALA	HB1		1HB	HB1	HB1	HB1	HB1		HB1	HB1	1HB
A	ALA	HB2		2HB	HB2	HB2	HB2	HB2		HB2	HB2	2HB
A	ALA	HB3		3HB	HB3	HB3	HB3	HB3		HB3	HB3	3HB
A	ALA	C		C		C	C	C	C	C	C	C
A	ALA	CA		CA		CA	CA	CA	CA	CA	CA	CA
A	ALA	CB		CB		CB	CB	CB	CB	CB	CB	CB
A	ALA	N		N		N	N	N	N	N	N	N
A	ALA	O		O		O	O	O	O	O	O	O

R	ARG	H		H	HN	HN	HN	H		HN	H	H
R	ARG	HA		HA	HA	HA	HA	HA		HA	HA	HA
R	ARG	HB2	pro-R	1HB	HB1	HB1	HB2	HB2		HB2	HB2	2HB
R	ARG	HB3	pro-S	2HB	HB2	HB2	HB1	HB1		HB3	HB1	1HB
R	ARG	HG2	pro-S	1HG	HG1	HG1	HG2	HG2		HG2	HG2	2HG
R	ARG	HG3	pro-R	2HG	HG2	HG2	HG1	HG1		HG3	HG1	1HG
R	ARG	HD2	pro-S	1HD	HD1	HD1	HD2	HD2		HD2	HD2	2HD
R	ARG	HD3	pro-R	2HD	HD2	HD2	HD1	HD1		HD3	HD1	1HD
R	ARG	HE		HE	HNE	HE	HE	HE		HE	HE	HE
R	ARG	HH11	Z	1HH1	HN11	HH11	HH11	HH11		HH11	1HH1	1HH1
R	ARG	HH12	E	2HH1	HN12	HH12	HH12	HH12		HH12	2HH1	2HH1
R	ARG	HH21	Z	1HH2	HN21	HH22	HH21	HH21		HH21	1HH2	1HH2
R	ARG	HH22	E	2HH2	HN22	HH21	HH22	HH22		HH22	2HH2	2HH2
R	ARG	C		C		C	C	C	C	C	C	C
R	ARG	CA		CA		CA	CA	CA	CA	CA	CA	CA
R	ARG	CB		CB		CB	CB	CB	CB	CB	CB	CB
R	ARG	CG		CG		CG	CG	CG	CG	CG	CG	CG
R	ARG	CD		CD		CD	CD	CD	CD	CD	CD	CD
R	ARG	CZ		CZ		CZ	CZ	CZ	CZ	CZ	CZ	CZ
R	ARG	N		N		N	N	N	N	N	N	N
R	ARG	NE		NE		NE	NE	NE	NE	NE	NE	NE
R	ARG	NH1	Z	NH1		NH1	NH1	NH1	NH1	NH1	NH1	NH1
R	ARG	NH2	E	NH2		NH2	NH2	NH2	NH2	NH2	NH2	NH2
R	ARG	O		O		O	O	O	O	O	O	O

D	ASP	H		H	HN	HN	HN	H		HN	H	H
D	ASP	HA		HA	HA	HA	HA	HA		HA	HA	HA
D	ASP	HB2	pro-S	1HB	HB1	HB1	HB2	HB2		HB2	HB2	2HB
D	ASP	HB3	pro-R	2HB	HB2	HB2	HB1	HB1		HB3	HB1	1HB
D	ASP	HD2				HD2				HD2
D	ASP	C		C		C	C	C	C	C	C	C
D	ASP	CA		CA		CA	CA	CA	CA	CA	CA	CA
D	ASP	CB		CB		CB	CB	CB	CB	CB	CB	CB
D	ASP	CG		CG		CG	CG	CG	CG	CG	CG	CG
D	ASP	N		N		N	N	N	N	N	N	N
D	ASP	O		O		O	O	O	O	O	O	O
D	ASP	OD1		OD1		OD1	OD1	OD1	OD1	OD1	OD1	OD1
D	ASP	OD2		OD2		OD2	OD2	OD2	OD2	OD2	OD2	OD2

N	ASN	H		H	HN	HN	HN	H		HN	H	H
N	ASN	HA		HA	HA	HA	HA	HA		HA	HA	HA
N	ASN	HB2	pro-S	1HB	HB1	HB1	HB2	HB2		HB2	HB2	2HB
N	ASN	HB3	pro-R	2HB	HB2	HB2	HB1	HB1		HB3	HB1	1HB
N	ASN	HD21	E	1HD2	HN21	HD21	HD21	HD21		HD21	1HD2	1HD2
N	ASN	HD22	Z	2HD2	HN22	HD22	HD22	HD22		HD22	2HD2	2HD2
N	ASN	C		C		C	C	C	C	C	C	C
N	ASN	CA		CA		CA	CA	CA	CA	CA	CA	CA
N	ASN	CB		CB		CB	CB	CB	CB	CB	CB	CB
N	ASN	CG		CG		CG	CG	CG	CG	CG	CG	CG
N	ASN	N		N		N	N	N	N	N	N	N
N	ASN	ND2		ND2		ND2	ND2	ND2	ND2	ND2	ND2	ND2
N	ASN	O		O		O	O	O	O	O	O	O
N	ASN	OD1		OD1		OD1	OD1	OD1	OD1	OD1	OD1	OD1

C	CYS	H		H	HN	HN	HN	H		HN	H	H
C	CYS	HA		HA	HA	HA	HA	HA		HA	HA	HA
C	CYS	HB2	pro-S	1HB	HB1	HB1	HB2	HB2		HB2	HB2	2HB
C	CYS	HB3	pro-R	2HB	HB2	HB2	HB1	HB1		HB3	HB1	1HB
C	CYS	HG		HG	HSG	HG	HG	HG		HG	HG	HG
C	CYS	C		C		C	C	C	C	C	C	C
C	CYS	CA		CA		CA	CA	CA	CA	CA	CA	CA
C	CYS	CB		CB		CB	CB	CB	CB	CB	CB	CB
C	CYS	N		N		N	N	N	N	N	N	N
C	CYS	O		O		O	O	O	O	O	O	O
C	CYS	SG		SG		SG	SG	SG	SG	SG	SG	SG

E	GLU	H		H	HN	HN	HN	H		HN	H	H
E	GLU	HA		HA	HA	HA	HA	HA		HA	HA	HA
E	GLU	HB2	pro-R	1HB	HB1	HB1	HB2	HB2		HB2	HB2	2HB
E	GLU	HB3	pro-S	2HB	HB2	HB2	HB1	HB1		HB3	HB1	1HB
E	GLU	HG2	pro-S	1HG	HG1	HG1	HG2	HG2		HG2	HG2	2HG
E	GLU	HG3	pro-R	2HG	HG2	HG2	HG1	HG1		HG3	HG1	1HG
E	GLU	HE2				HE2				HE2
E	GLU	C		C		C	C	C	C	C	C	C
E	GLU	CA		CA		CA	CA	CA	CA	CA	CA	CA
E	GLU	CB		CB		CB	CB	CB	CB	CB	CB	CB
E	GLU	CG		CG		CG	CG	CG	CG	CG	CG	CG
E	GLU	CD		CD		CD	CD	CD	CD	CD	CD	CD
E	GLU	N		N		N	N	N	N	N	N	N
E	GLU	O		O		O	O	O	O	O	O	O
E	GLU	OE1		OE1		OE1	OE1	OE1	OE1	OE1	OE1	OE1
E	GLU	OE2		OE2		OE2	OE2	OE2	OE2	OE2	OE2	OE2

Q	GLN	H		H	HN	HN	HN	H		HN	H	H
Q	GLN	HA		HA	HA	HA	HA	HA		HA	HA	HA
Q	GLN	HB2	pro-R	1HB	HB1	HB1	HB2	HB2		HB2	HB2	2HB
Q	GLN	HB3	pro-S	2HB	HB2	HB2	HB1	HB1		HB3	HB1	1HB
Q	GLN	HG2	pro-S	1HG	HG1	HG1	HG2	HG2		HG2	HG2	2HG
Q	GLN	HG3	pro-R	2HG	HG2	HG2	HG1	HG1		HG3	HG1	1HG
Q	GLN	HE21	E	1HE2	HN21	HE21	HE21	HE21		HE21	1HE2	1HE2
Q	GLN	HE22	Z	2HE2	HN22	HE22	HE22	HE22		HE22	2HE2	2HE2
Q	GLN	C		C		C	C	C	C	C	C	C
Q	GLN	CA		CA		CA	CA	CA	CA	CA	CA	CA
Q	GLN	CB		CB		CB	CB	CB	CB	CB	CB	CB
Q	GLN	CG		CG		CG	CG	CG	CG	CG	CG	CG
Q	GLN	CD		CD		CD	CD	CD	CD	CD	CD	CD
Q	GLN	N		N		N	N	N	N	N	N	N
Q	GLN	NE2		NE2		NE2	NE2	NE2	NE2	NE2	NE2	NE2
Q	GLN	O		O		O	O	O	O	O	O	O
Q	GLN	OE1		OE1		OE1	OE1	OE1	OE1	OE1	OE1	OE1

G	GLY	H		H	HN	HN	HN	H		HN	H	H
G	GLY	HA2	pro-R	1HA	HA1	HA1	HA2	HA2		HA1	HA2	2HA
G	GLY	HA3	pro-S	2HA	HA2	HA2	HA1	HA1		HA2	HA1	1HA
G	GLY	C		C		C	C	C	C	C	C	C
G	GLY	CA		CA		CA	CA	CA	CA	CA	CA	CA
G	GLY	N		N		N	N	N	N	N	N	N
G	GLY	O		O		O	O	O	O	O	O	O

H	HIS	H		H	HN	HN	HN	H		HN	H	H
H	HIS	HA		HA	HA	HA	HA	HA		HA	HA	HA
H	HIS	HB2	pro-S	1HB	HB1	HB1	HB2	HB2		HB2	HB2	2HB
H	HIS	HB3	pro-R	2HB	HB2	HB2	HB1	HB1		HB3	HB1	1HB
H	HIS	HD1		HD1	HND1	HD1	HD1	HD1		HD1	HD1	1HD
H	HIS	HD2		HD2	HD2	HD2	HD2	HD2		HD2	HD2	2HD
H	HIS	HE1		HE1	HE1	HE1	HE1	HE1		HE1	HE1	1HE
H	HIS	HE2			HNE2	HE2				HE2	HE2	2HE
H	HIS	C		C		C	C	C	C	C	C	C
H	HIS	CA		CA		CA	CA	CA	CA	CA	CA	CA
H	HIS	CB		CB		CB	CB	CB	CB	CB	CB	CB
H	HIS	CG		CG		CG	CG	CG	CG	CG	CG	CG
H	HIS	CD2		CD2		CD2	CD2	CD2	CD2	CD2	CD2	CD2
H	HIS	CE1		CE1		CE1	CE1	CE1	CE1	CE1	CE1	CE1
H	HIS	N		N		N	N	N	N	N	N	N
H	HIS	ND1		ND1		ND1	ND1	ND1	ND1	ND1	ND1	ND1
H	HIS	NE2		NE2		NE2	NE2	NE2	NE2	NE2	NE2	NE2
H	HIS	O		O		O	O	O	O	O	O	O

I	ILE	H		H	HN	HN	HN	H		HN	H	H
I	ILE	HA		HA	HA	HA	HA	HA		HA	HA	HA
I	ILE	HB		HB	HB	HB	HB	HB		HB	HB	HB
I	ILE	HG12	pro-R	1HG1	HG11	HG11	HG12	HG12		HG12	2HG1	2HG1
I	ILE	HG13	pro-S	2HG1	HG12	HG12	HG11	HG11		HG13	1HG1	1HG1
I	ILE	HG21		1HG2	HG21	HG21	HG21	HG21		HG21	1HG2	1HG2
I	ILE	HG22		2HG2	HG22	HG22	HG22	HG22		HG22	2HG2	2HG2
I	ILE	HG23		3HG2	HG23	HG23	HG23	HG23		HG23	3HG2	3HG2
I	ILE	HD11		1HD1	HD11	HD11	HD11	HD11		HD11	HD1	1HD 
I	ILE	HD12		2HD1	HD12	HD12	HD12	HD12		HD12	HD2	2HD 
I	ILE	HD13		3HD1	HD13	HD13	HD13	HD13		HD13	HD3	3HD 
I	ILE	C		C		C	C	C	C	C	C	C
I	ILE	CA		CA		CA	CA	CA	CA	CA	CA	CA
I	ILE	CB		CB		CB	CB	CB	CB	CB	CB	CB
I	ILE	CG1		CG1		CG1	CG1	CG1	CG1	CG1	CG1	CG1
I	ILE	CG2		CG2		CG2	CG2	CG2	CG2	CG2	CG2	CG2
I	ILE	CD1		CD1		CD1	CD1	CD1	CD1	CD1	CD1	CD
I	ILE	N		N		N	N	N	N	N	N	N
I	ILE	O		O		O	O	O	O	O	O	O

L	LEU	H		H	HN	HN	HN	H		HN	H	H
L	LEU	HA		HA	HA	HA	HA	HA		HA	HA	HA
L	LEU	HB2	pro-R	1HB	HB1	HB1	HB2	HB2		HB2	HB2	2HB
L	LEU	HB3	pro-S	2HB	HB2	HB2	HB1	HB1		HB3	HB1	1HB
L	LEU	HG		HG	HG	HG	HG	HG		HG	HG	HG
L	LEU	HD11		1HD1	HD11	HD11	HD11	HD11		HD11	1HD1	1HD1
L	LEU	HD12		2HD1	HD12	HD12	HD12	HD12		HD12	2HD1	2HD1
L	LEU	HD13		3HD1	HD13	HD13	HD13	HD13		HD13	3HD1	3HD1
L	LEU	HD21		1HD2	HD21	HD21	HD21	HD21		HD21	1HD2	1HD2
L	LEU	HD22		2HD2	HD22	HD22	HD22	HD22		HD22	2HD2	2HD2
L	LEU	HD23		3HD2	HD23	HD23	HD23	HD23		HD23	3HD2	3HD2
L	LEU	C		C		C	C	C	C	C	C	C
L	LEU	CA		CA		CA	CA	CA	CA	CA	CA	CA
L	LEU	CB		CB		CB	CB	CB	CB	CB	CB	CB
L	LEU	CG		CG		CG	CG	CG	CG	CG	CG	CG
L	LEU	CD1	pro-R	CD1		CD1	CD1	CD1	CD1	CD1	CD1	CD1
L	LEU	CD2	pro-S	CD2		CD2	CD2	CD2	CD2	CD2	CD2	CD2
L	LEU	N		N		N	N	N	N	N	N	N
L	LEU	O		O		O	O	O	O	O	O	O

K	LYS	H		H	HN	HN	HN	H		HN	H	H
K	LYS	HA		HA	HA	HA	HA	HA		HA	HA	HA
K	LYS	HB2	pro-R	1HB	HB1	HB1	HB2	HB2		HB2	HB2	2HB
K	LYS	HB3	pro-S	2HB	HB2	HB2	HB1	HB1		HB3	HB1	1HB
K	LYS	HG2	pro-R	1HG	HG1	HG1	HG2	HG2		HG2	HG2	2HG
K	LYS	HG3	pro-S	2HG	HG2	HG2	HG1	HG1		HG3	HG1	1HG
K	LYS	HD2	pro-S	1HD	HD1	HD1	HD2	HD2		HD2	HD2	2HD
K	LYS	HD3	pro-R	2HD	HD2	HD2	HD1	HD1		HD3	HD1	1HD
K	LYS	HE2	pro-S	1HE	HE1	HE1	HE2	HE2		HE2	HE2	2HE
K	LYS	HE3	pro-R	2HE	HE2	HE2	HE1	HE1		HE3	HE1	1HE
K	LYS	HZ1		1HZ	HNZ1	HZ1	HZ1	HZ1		HZ1	HZ1	1HZ
K	LYS	HZ2		2HZ	HNZ2	HZ2	HZ2	HZ2		HZ2	HZ2	2HZ
K	LYS	HZ3		3HZ	HNZ3	HZ3	HZ3	HZ3		HZ3	HZ3	3HZ
K	LYS	C		C		C	C	C	C	C	C	C
K	LYS	CA		CA		CA	CA	CA	CA	CA	CA	CA
K	LYS	CB		CB		CB	CB	CB	CB	CB	CB	CB
K	LYS	CG		CG		CG	CG	CG	CG	CG	CG	CG
K	LYS	CD		CD		CD	CD	CD	CD	CD	CD	CD
K	LYS	CE		CE		CE	CE	CE	CE	CE	CE	CE
K	LYS	N		N		N	N	N	N	N	N	N
K	LYS	NZ		NZ		NZ	NZ	NZ	NZ	NZ	NZ	NZ
K	LYS	O		O		O	O	O	O	O	O	O

M	MET	H		H	HN	HN	HN	H		HN	H	H
M	MET	HA		HA	HA	HA	HA	HA		HA	HA	HA
M	MET	HB2	pro-S	1HB	HB1	HB1	HB2	HB2		HB2	HB2	2HB
M	MET	HB3	pro-R	2HB	HB2	HB2	HB1	HB1		HB3	HB1	1HB
M	MET	HG2	pro-S	1HG	HG1	HG1	HG2	HG2		HG2	HG2	2HG
M	MET	HG3	pro-R	2HG	HG2	HG2	HG1	HG1		HG3	HG1	1HG
M	MET	HE1		1HE	HE1	HE1	HE1	HE1		HE1	HE1	1HE
M	MET	HE2		2HE	HE2	HE2	HE2	HE2		HE2	HE2	2HE
M	MET	HE3		3HE	HE3	HE3	HE3	HE3		HE3	HE3	3HE
M	MET	C		C		C	C	C	C	C	C	C
M	MET	CA		CA		CA	CA	CA	CA	CA	CA	CA
M	MET	CB		CB		CB	CB	CB	CB	CB	CB	CB
M	MET	CG		CG		CG	CG	CG	CG	CG	CG	CG
M	MET	CE		CE		CE	CE	CE	CE	CE	CE	CE
M	MET	N		N		N	N	N	N	N	N	N
M	MET	O		O		O	O	O	O	O	O	O
M	MET	SD		SD		SD	SD	SD	SD	SD	SD	SD

F	PHE	H		H	HN	HN	HN	H		HN	H	H
F	PHE	HA		HA	HA	HA	HA	HA		HA	HA	HA
F	PHE	HB2	pro-R	1HB	HB1	HB1	HB2	HB2		HB2	HB2	2HB
F	PHE	HB3	pro-S	2HB	HB2	HB2	HB1	HB1		HB3	HB1	1HB
F	PHE	HD1		HD1	HD1	HD1	HD1	HD1		HD1	HD1	1HD
F	PHE	HD2		HD2	HD2	HD2	HD2	HD2		HD2	HD2	2HD
F	PHE	HE1		HE1	HE1	HE1	HE1	HE1		HE1	HE1	1HE
F	PHE	HE2		HE2	HE2	HE2	HE2	HE2		HE2	HE2	2HE
F	PHE	HZ		HZ	HZ	HZ	HZ	HZ		HZ	HZ	HZ
F	PHE	C		C		C	C	C	C	C	C	C
F	PHE	CA		CA		CA	CA	CA	CA	CA	CA	CA
F	PHE	CB		CB		CB	CB	CB	CB	CB	CB	CB
F	PHE	CG		CG		CG	CG	CG	CG	CG	CG	CG
F	PHE	CD1		CD1		CD1	CD1	CD1	CD1	CD1	CD1	CD1
F	PHE	CD2		CD2		CD2	CD2	CD2	CD2	CD2	CD2	CD2
F	PHE	CE1		CE1		CE1	CE1	CE1	CE1	CE1	CE1	CE1
F	PHE	CE2		CE2		CE2	CE2	CE2	CE2	CE2	CE2	CE2
F	PHE	CZ		CZ		CZ	CZ	CZ	CZ	CZ	CZ	CZ
F	PHE	N		N		N	N	N	N	N	N	N
F	PHE	O		O		O	O	O	O	O	O	O

P	PRO	H2	pro-R	H2		HN2	HT2					
P	PRO	H3	pro-S	H1		HN1	HT1					
P	PRO	HA		HA	HA	HA	HA	HA		HA	HA	HA
P	PRO	HB2	pro-R	1HB	HB1	HB1	HB2	HB2		HB2	HB2	2HB
P	PRO	HB3	pro-S	2HB	HB2	HB2	HB1	HB1		HB3	HB1	1HB
P	PRO	HG2	pro-S	1HG	HG1	HG1	HG2	HG2		HG2	HG2	2HG
P	PRO	HG3	pro-R	2HG	HG2	HG2	HG1	HG1		HG3	HG1	1HG
P	PRO	HD2	pro-S	1HD	HD2	HD1	HD2	HD2		HD2	HD2	2HD
P	PRO	HD3	pro-R	2HD	HD1	HD2	HD1	HD1		HD3	HD1	1HD
P	PRO	C		C		C	C	C	C	C	C	C
P	PRO	CA		CA		CA	CA	CA	CA	CA	CA	CA
P	PRO	CB		CB		CB	CB	CB	CB	CB	CB	CB
P	PRO	CG		CG		CG	CG	CG	CG	CG	CG	CG
P	PRO	CD		CD		CD	CD	CD	CD	CD	CD	CD
P	PRO	N		N		N	N	N	N	N	N	N
P	PRO	O		O		O	O	O	O	O	O	O

S	SER	H		H	HN	HN	HN	H		HN	H	H
S	SER	HA		HA	HA	HA	HA	HA		HA	HA	HA
S	SER	HB2	pro-S	1HB	HB1	HB1	HB2	HB2		HB2	HB2	2HB
S	SER	HB3	pro-R	2HB	HB2	HB2	HB1	HB1		HB3	HB1	1HB
S	SER	HG		HG	HOG	HG	HG	HG		HG	HG	HG
S	SER	C		C		C	C	C	C	C	C	C
S	SER	CA		CA		CA	CA	CA	CA	CA	CA	CA
S	SER	CB		CB		CB	CB	CB	CB	CB	CB	CB
S	SER	N		N		N	N	N	N	N	N	N
S	SER	O		O		O	O	O	O	O	O	O
S	SER	OG		OG		OG	OG	OG	OG	OG	OG	OG

T	THR	H		H	HN	HN	HN	H		HN	H	H
T	THR	HA		HA	HA	HA	HA	HA		HA	HA	HA
T	THR	HB		HB	HB	HB	HB	HB		HB	HB	HB
T	THR	HG1		HG1	HOG1	HG1	HG1	HG1		HG1	HG1	1HG
T	THR	HG21		1HG2	HG21	HG21	HG21	HG21		HG21	1HG2	1HG2
T	THR	HG22		2HG2	HG22	HG22	HG22	HG22		HG22	2HG2	2HG2
T	THR	HG23		3HG2	HG23	HG23	HG23	HG23		HG23	3HG2	3HG2
T	THR	C		C		C	C	C	C	C	C	C
T	THR	CA		CA		CA	CA	CA	CA	CA	CA	CA
T	THR	CB		CB		CB	CB	CB	CB	CB	CB	CB
T	THR	CG2		CG2		CG2	CG2	CG2	CG2	CG2	CG2	CG2
T	THR	N		N		N	N	N	N	N	N	N
T	THR	O		O		O	O	O	O	O	O	O
T	THR	OG1		OG1		OG1	OG1	OG1	OG1	OG1	OG1	OG1

W	TRP	H		H	HN	HN	HN	H		HN	H	H
W	TRP	HA		HA	HA	HA	HA	HA		HA	HA	HA
W	TRP	HB2	pro-R	1HB	HB1	HB1	HB2	HB2		HB2	HB2	2HB
W	TRP	HB3	pro-S	2HB	HB2	HB2	HB1	HB1		HB3	HB1	1HB
W	TRP	HD1		HD1	HD1	HD1	HD1	HD1		HD1	HD1	HD
W	TRP	HE1		HE1	HNE1	HE1	HE1	HE1		HE1	HE1	1HE
W	TRP	HE3		HE3	HE3	HE3	HE3	HE3		HE3	HE3	3HE
W	TRP	HZ2		HZ2	HZ2	HZ2	HZ2	HZ2		HZ2	HZ2	2HZ
W	TRP	HZ3		HZ3	HZ3	HZ3	HZ3	HZ3		HZ3	HZ3	3HZ
W	TRP	HH2		HH2	HH2	HH2	HH2	HH2		HH2	HH2	HH
W	TRP	C		C		C	C	C	C	C	C	C
W	TRP	CA		CA		CA	CA	CA	CA	CA	CA	CA
W	TRP	CB		CB		CB	CB	CB	CB	CB	CB	CB
W	TRP	CG		CG		CG	CG	CG	CG	CG	CG	CG
W	TRP	CD1		CD1		CD1	CD1	CD1	CD1	CD1	CD1	CD1
W	TRP	CD2		CD2		CD2	CD2	CD2	CD2	CD2	CD2	CD2
W	TRP	CE2		CE2		CE2	CE2	CE2	CE2	CE2	CE2	CE2
W	TRP	CE3		CE3		CE3	CE3	CE3	CE3	CE3	CE3	CE3
W	TRP	CZ2		CZ2		CZ2	CZ2	CZ2	CZ2	CZ2	CZ2	CZ2
W	TRP	CZ3		CZ3		CZ3	CZ3	CZ3	CZ3	CZ3	CZ3	CZ3
W	TRP	CH2		CH2		CH2	CH2	CH2	CH2	CH2	CH2	CH2
W	TRP	N		N		N	N	N	N	N	N	N
W	TRP	NE1		NE1		NE1	NE1	NE1	NE1	NE1	NE1	NE1
W	TRP	O		O		O	O	O	O	O	O	O

Y	TYR	H		H	HN	HN	HN	H		HN	H	H
Y	TYR	HA		HA	HA	HA	HA	HA		HA	HA	HA
Y	TYR	HB2	pro-R	1HB	HB1	HB1	HB2	HB2		HB2	HB2	2HB
Y	TYR	HB3	pro-S	2HB	HB2	HB2	HB1	HB1		HB3	HB1	1HB
Y	TYR	HD1		HD1	HD1	HD1	HD1	HD1		HD1	HD1	1HD
Y	TYR	HD2		HD2	HD2	HD2	HD2	HD2		HD2	HD2	2HD
Y	TYR	HE1		HE1	HE1	HE1	HE1	HE1		HE1	HE1	1HE
Y	TYR	HE2		HE2	HE2	HE2	HE2	HE2		HE2	HE2	2HE
Y	TYR	HH		HH	HOH	HH	HH	HH		HH	HH	HH
Y	TYR	C		C		C	C	C	C	C	C	C
Y	TYR	CA		CA		CA	CA	CA	CA	CA	CA	CA
Y	TYR	CB		CB		CB	CB	CB	CB	CB	CB	CB
Y	TYR	CG		CG		CG	CG	CG	CG	CG	CG	CG
Y	TYR	CD1		CD1		CD1	CD1	CD1	CD1	CD1	CD1	CD1
Y	TYR	CD2		CD2		CD2	CD2	CD2	CD2	CD2	CD2	CD2
Y	TYR	CE1		CE1		CE1	CE1	CE1	CE1	CE1	CE1	CE1
Y	TYR	CE2		CE2		CE2	CE2	CE2	CE2	CE2	CE2	CE2
Y	TYR	CZ		CZ		CZ	CZ	CZ	CZ	CZ	CZ	CZ
Y	TYR	N		N		N	N	N	N	N	N	N
Y	TYR	O		O		O	O	O	O	O	O	O
Y	TYR	OH		OH		OH	OH	OH	OH	OH	OH	OH

V	VAL	H		H	HN	HN	HN	H		HN	H	H
V	VAL	HA		HA	HA	HA	HA	HA		HA	HA	HA
V	VAL	HB		HB	HB	HB	HB	HB		HB	HB	HB
V	VAL	HG11		1HG1	HG11	HG21	HG11	HG11		HG11	1HG1	1HG1
V	VAL	HG12		2HG1	HG12	HG22	HG12	HG12		HG12	2HG1	2HG1
V	VAL	HG13		3HG1	HG13	HG23	HG13	HG13		HG13	3HG1	3HG1
V	VAL	HG21		1HG2	HG21	HG11	HG21	HG21		HG21	1HG2	1HG2
V	VAL	HG22		2HG2	HG22	HG12	HG22	HG22		HG22	2HG2	2HG2
V	VAL	HG23		3HG2	HG23	HG13	HG23	HG23		HG23	3HG2	3HG2
V	VAL	C		C		C	C	C	C	C	C	C
V	VAL	CA		CA		CA	CA	CA	CA	CA	CA	CA
V	VAL	CB		CB		CB	CB	CB	CB	CB	CB	CB
V	VAL	CG1	pro-R	CG1		CG2	CG1	CG1	CG1	CG1	CG1	CG1
V	VAL	CG2	pro-S	CG2		CG1	CG2	CG2	CG2	CG2	CG2	CG2
V	VAL	N		N		N	N	N	N	N	N	N
V	VAL	O		O		O	O	O	O	O	O	O
";
}#_sub Init
