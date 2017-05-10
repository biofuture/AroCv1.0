#/usr/bin/env perl -w
use strict;

##This script is developed to fastly and accurately calculate the average rrn operon copy (AroC) number in 16s amplicion sequencing data or metagenomics sequencing data

##1. for 16s amplicon data, the input is the OTUs mothur table, the aroc is calculated by community structure and CopyRighter copy number database 
##2. for metagenomics data, the input is the orignial high quality fastq reads, there are two methods for metagenomics data. This first one is to extract 16S sequences from metagenomics data and then using these 16S information to get the aroc with method metioned in 1. The second method is to utilize universal single-copy marker gene to get the averagly coverage of these marker genes and then estimate the total possible cell number, then we using the total number of 16S to normalize aganist the cell number to get the aroc
#developed by Xiaotao-Jiang on 2016-10-24 at the University of Hong Kong 
#Email: biofuture.jiang@gmail.com

#------------------------------------------------------------------------------------------------------
use Getopt::Std;
use File::Basename;
use FindBin qw($Bin);

my $AroCDir;

BEGIN {
    my @dirset = split(/\//,$Bin);
    #pop(@dirset); pop(@dirset);
    $AroCDir = join("/", @dirset);
    unshift @INC, "$AroCDir";
}
#print "$AroCDir\n";

my $usage =<<USG;
    
    Author Xiao-tao JIANG
    Email: biofuture.jiang\@gmail.com
    Date 24/10/2016
		
===>For 16S amplicon microbial communities
    perl $0 -t <otu-table-with-tax> -o <output.aroc.txt>  -h [help information]

    -t the otu table with columns are samples and rows are OTUs, the end column should be greengen format taxonomy annotation for that OTUs
    -o the text format file contains aroc for each sample in otu table
 
===>For metagenomics sequencing data	
    perl $0 -1 <1.fq> -2 <2.fq> -U <single.fq>  -m <s/m> -o <aroc.txt> -a [bowtie2/usearch] -b [bowtie2/diamond] -n [threads number] -h [help information] 	

    -m s/m; s represnts by extracting 16S from metagenomics data and then using the Copyrighter copy number database to calculate the aroc; m represents using universal single-copy phylogenetic marker gene to calculate the aroc (default m)  
    -1 pair-end fastq/fasta 1
    -2 pair-end fastq.fasta 2
    -U single-end fastq/fasta		

    For metagenomics data, the work is composed of two parts, the first part is searching reads from 16SrDNA genes, to count the total full length 16S numbers or to extract the hyper variable regions to g et the community composition and abundance information. The second part is to extract those reads orignial from 40 universal single copy phylogenetic marker genes (COGs, ref: ) or 30 essential single copy genes (KOs, ref: ). 

    Alignment methods for extracting 16S including: usearch and bowtie2 
    Methods to calculate the essential single copy genes: bowtie2, or diamond (blastx searching) 
    This version -a and -b are blocked.      
	
EXample: 
    perl $0 -1 test_1.fq -2 test_2.fq -m m -o test.aroc.txt -n 8
   
    perl $0 -1 test_1.fq -2 test_2.fq -U test.single.fq -m m -o test.aroc.txt -n 8 -f fq	
	
    By extracting 16S hyper variable regions and CopyRighter database to estimate the ARCN 
    perl $0 -1 test_1.fq -2 test_2.fq -m s -o test.aroc.txt -n 8 
USG

our($opt_1, $opt_2, $opt_U, $opt_t, $opt_n, $opt_m, $opt_o, $opt_f, $opt_a, $opt_b, $opt_s, $opt_h)="";
getopts('1:2:U:n:m:t:o:f:a:b:s:h');

my $date = localtime;
$opt_n ||= 8;
$opt_1 ||= "";
$opt_2 ||= "";
$opt_U ||= "";
$opt_a ||= "usearch";
$opt_b ||= "diamond";
$opt_f ||= "fq";
$opt_s ||=100;
my $evu ||=1e-4;
my $eval ||=6;

if($opt_h){

	print "$usage\n";
	exit;
}

##---------------------------------------------------------------------------------------
##Database prepare for the pipelie 

my $REFHVR6 = "$AroCDir/db/RefHVR.V6.udb"; #usearch index for V6, nuc sequences 
my $ID2TAXV6 = "$AroCDir/db/RefHVR.V6.taxonomy.txt";
my $ggnr85 = "$AroCDir/db/gg85.udb"; ##usearch index, nuc sequences 
my $ggnr85index = "$AroCDir/db/template_85_16S"; ##greengene 85% idnetity OTUs 16S bowtie2 index, nuc sequences 

##This is the 39 COGs (Orignial from SpecieI) from 5850 complete genomes from NCBI bacteria database
my $duscmgbowtie2index = "$AroCDir/db/uscmg.bowtie2index"; ##bowttie2 index  for 39 COGs nucletide sequences 
my $cogsnamelist = "$AroCDir/db/all_COGs_nuc.name.map.list"; ##gene -> COGs name list 

##This is the 30 essential single copy KOs across Bacteria, Archaea and Fungi with at least 99% coverage
my $dbindex = "$AroCDir/db/KO30_DIAMOND.dmnd"; #diamond index for 30 KOs aa sequences 
my $cogslis = "$AroCDir/db/all_KO30_name.list";#gene->KO name list 

##import Copy number database into hash table 
my $cprdb = "$AroCDir/db/Copy_db.copyrighter.txt";
die "$!\n" unless open(CPDB, "$cprdb");
my %copyn;
<CPDB>;
while(<CPDB>){
        chomp;
        my @tem = split(/\t/, $_);
        $copyn{$tem[0]} = $tem[1];
}
close CPDB;

my %id2tax;
die "$!\n"unless open(ID2TAX, "$ID2TAXV6");
while(<ID2TAX>){
        chomp;
        my @tem = split("\t", $_);
        $id2tax{$tem[0]} = $tem[1];
}
close ID2TAX;


my @blank = ("k__", "p__", "c__", "o__", "f__", "g__", "s__"); ##this is for lca
##------------------------------------------------------------------------------------------

if($opt_t){
	#die "perl $0 <otu.matrx> <cnv.out>\n" unless (@ARGV == 2);
	otutablearoc($opt_t, $opt_o);

}elsif($opt_1 || $opt_2 || $opt_U){ ## single, pair-end and pair-end with single are recived as parameters 
	##Extract all 16S reads from metagenomics samples and estimate the total number of full length 16S rRNA gene in this metagenomis sample
	
	my $out16sreads = "out.reads.16s";
	my $usv6 = "out.us.v6";
	my $comfile = "out.v6.community";
	my $outcov = "out.coverage.txt";
	
	if($opt_U && $opt_1 eq ""){

		if($opt_m eq "s"){

		##extract all hyper variable region reads 
		##using hyper variable regions to get the microbial community information; applying copyrigher information to calculate the aroc 
		
		extract_hvr($opt_U, $usv6);
		hyper_community($usv6, $comfile);
		aroc_16scommunity($comfile, $opt_o);


		}elsif($opt_m eq "m"){

			#Mapping to a collection of universal single-copy phylogenetic marker genes
			#using bowtie2 to map reads agianst the seletect universal pMGs	
			#Get the average coverage for all these universal marker genes, taken this value as cell number estimation and calculate the aroc by dividing the cell number with 16S number 

			if($opt_a eq "usearch"){
				extract_16s_usearch($opt_U, $out16sreads);
			
			}elsif($opt_a eq "bowtie2"){
				extrac_16s_bowtie2($opt_U, $out16sreads, $opt_f);
			}else{
				die "Wrong line 134";
			}
			
			#diamond 
			map2scpMGs($opt_U, $outcov, $opt_f);			
			get_aroc($outcov, $out16sreads, $opt_o);


		}else{
			print "None support format input for -m\n";
			exit;
		}

	}elsif($opt_1 && $opt_2 && $opt_U eq ""){
		if($opt_m eq "s"){

		##extract all hyper variable region reads 
		##using hyper variable regions to get the microbial community information; applying copyrigher information to calculate the aroc 
		#hyper_community($out16sreads, $community);
		#aroc_16scommunity($community, $outcopy);
		extract_hvr($opt_1, $opt_2, $usv6);
		hyper_community($usv6, $comfile);
		aroc_16scommunity($comfile, $opt_o);


		}elsif($opt_m eq "m"){

			#Mapping to a collection of universal single-copy phylogenetic marker genes
			#using bowtie2 to map reads agianst the seletect universal pMGs	
			#Get the average coverage for all these universal marker genes, taken this value as cell number estimation and calculate the aroc by dividing the cell number with 16S number 
			
                        if($opt_a eq "usearch"){

				extract_16s_usearch($opt_1, $opt_2, $out16sreads);

			}elsif($opt_a eq "bowtie2"){

				extract_16s_bowtie2($opt_1, $opt_2, $out16sreads, $opt_f);

			}else{
				die "Wrong line 171\n";
			}
			map2scpMGs($opt_1, $opt_2,$outcov, $opt_f);
			get_aroc($outcov, $out16sreads, $opt_o);


		}else{
			print "None support format input for -m\n";
			exit;
		}


	
	}elsif($opt_1 && $opt_2 && $opt_U){
		if($opt_m eq "s"){

		##extract all hyper variable region reads 
		##using hyper variable regions to get the microbial community information; applying copyrigher information to calculate the aroc 
		#hyper_community($out16sreads, $community);
		#aroc_16scommunity($community, $outcopy);
		extract_hvr($opt_1, $opt_2,$opt_U, $usv6);
		hyper_community($usv6, $comfile);
		aroc_16scommunity($comfile, $opt_o);


		}elsif($opt_m eq "m"){

			#Mapping to a collection of universal single-copy phylogenetic marker genes
			#using bowtie2 to map reads agianst the seletect universal pMGs	
			#Get the average coverage for all these universal marker genes, taken this value as cell number estimation and calculate the aroc by dividing the cell number with 16S number 
			if($opt_a eq "usearch"){
				extract_16s_usearch($opt_1, $opt_2, $opt_U,$out16sreads);
			}elsif($opt_a eq "bowtie2"){
				extract_16s_bowtie2($opt_1, $opt_2, $opt_U,$out16sreads, $opt_f);
			}else{
				die "Wrong line 203\n";
			}
			map2scpMGs($opt_1, $opt_2, $opt_U, $outcov, $opt_f);
			get_aroc($outcov, $out16sreads, $opt_o);


		}else{
			print "None support format input for -m\n";
			exit;
		}

	}


}else{
	print $usage;
	exit;	
}

#--------------------------------------------------------------------------------------------
##This function extract 16S reads from metagenomics samples and retun the reads file 

##extract_16s diamond 

sub extract_16s_usearch{

	my @tem = @_;	
	
	if($#tem == 1){
		##single end sequences  
		`$AroCDir/bin/usearch -ublast $tem[0] -db $ggnr85 -evalue $evu -id 0.85 -accel 0.5 -matched $tem[1] -threads $opt_n -strand both  -maxaccepts 1`;

	}elsif($#tem == 2 ){
		##pair end sequences 
		my $temr1 = "temp.1.16s.reads";
		my $temr2 = "temp.2.16s.reads";
		`$AroCDir/bin/usearch -ublast $tem[0]  -db $ggnr85 -evalue $evu -id 0.85 -accel 0.5 -matched $temr1 -threads $opt_n -strand both  -maxaccepts 1`;
		`$AroCDir/bin/usearch -ublast $tem[1]  -db $ggnr85 -evalue $evu -id 0.85 -accel 0.5 -matched $temr2 -threads $opt_n -strand both  -maxaccepts 1`;
		`cat $temr1 $temr2 > $tem[2]`;

	}elsif($#tem == 3){
		##both pair end and signle end sequences 
		my $temr1 = "temp.1.16s.reads";
		my $temr2 = "temp.2.16s.reads";
		my $temr3 = "temp.3.16s.reads";
		
		`$AroCDir/bin/usearch -ublast $tem[0] -db $ggnr85 -evalue $evu -id 0.85 -accel 0.5 -matched $temr1  -threads $opt_n -strand both  -maxaccepts 1`;
		`$AroCDir/bin/usearch -ublast $tem[1]  -db $ggnr85 -evalue $evu -id 0.85 -accel 0.5 -matched $temr2 -threads $opt_n -strand both  -maxaccepts 1`;
		`$AroCDir/bin/usearch -ublast $tem[2]  -db $ggnr85 -evalue $evu -id 0.85 -accel 0.5 -matched $temr3 -threads $opt_n -strand both  -maxaccepts 1`;

		`cat $temr1 $temr2 $temr3 > $tem[3]`;
	}else{
		print "Wrong number of sequences for metagenomics sequences\no";
		exit;
	}
}#extract_16S_usearch

##extract 16s bowtie2
sub extract_16s_bowtie2{

	my @tem = @_;

	if($#tem == 2){
		##single end sequences  
		my $tems = "temp.align.single";
		my $temp = "temp.align.pair";
		my $osam = "temp.16s.sam";
		my $temp1 = "temp.pair.1.align";
		my $temp2 = "temp.pair.2.align";
		if($tem[2] eq "fq"){
			`$AroCDir/bin/bowtie2 -x $ggnr85index  -U $tem[0] -S $osam -p $opt_n --no-unal --al $tems  --sensitive`;
			fastq2fasta($temp1, $temp2, $tems, $tem[1]);

		}elsif($tem[2] eq "fa"){
			`$AroCDir/bin/bowtie2 -x $ggnr85index  -U $tem[0] -S $osam -p $opt_n --no-unal --al $tems  --sensitive -f`;
			`cat $tems > $tem[1]`;
		}else{
			die "Wrong $tem[2]\n";
		}

	}elsif($#tem == 3 ){
		##pair end sequences 
		my $tems = "temp.single.align";
		my $temp = "temp.pair.align";
		my $temp1 = "temp.pair.1.align";
		my $temp2 = "temp.pair.2.align";
		my $osam = "temp.16s.sam";

		if($tem[3] eq "fq"){
			#die "$AroCDir/bin/bowtie2 -x $ggnr85index  -1 $tem[0] -2 $tem[1] -S $osam -p $opt_n --no-unal --al $tems --al-conc $temp --sensitive\n";
			`$AroCDir/bin/bowtie2 -x $ggnr85index  -1 $tem[0] -2 $tem[1] -S $osam -p $opt_n --no-unal --al $tems --al-conc $temp --sensitive`;		
	
			fastq2fasta($temp1, $temp2, $tems, $tem[2]);
		}elsif($tem[3] eq "fa"){
			`$AroCDir/bin/bowtie2 -x $ggnr85index  -1 $tem[0] -2 $tem[1] -S $osam -p $opt_n --no-unal --al $tems --al-conc $temp --sensitive -f`;			
			`cat $temp1 $temp2 $tems > $tem[2]`;
		}else{
			die "Wrong $tem[3]\n";
		}

	}elsif($#tem == 4){
		##both pair end and signle end sequences 
		my $tems = "temp.single.align";
		my $temp = "temp.pair.align";
		my $temp1 = "temp.pair.1.align";
		my $temp2 = "temp.pair.2.align";
		my $osam = "temp.16s.sam";
		
		`$AroCDir/bin/bowtie2 -x $ggnr85index  -1 $tem[0] -2 $tem[1] -U $tem[2] -S $osam -p $opt_n --no-unal --al $tems --al-conc $temp --sensitive -f`;		
		fastq2fasta($temp1, $temp2, $tems, $tem[3]);

	}else{
		print "Wrong number of sequences for metagenomics sequences\n";
		exit;
	}
	
}

sub fastq2fasta{
	##This script transform fastq format 16S reads into 16S fasta file
	
	my @files = @_;

	die "$!\n" unless open(OUTFA, ">$files[3]");
	die "$!\n" unless open(P1, "$files[0]");
	die "$!\n" unless open(P2, "$files[1]");
	die "$!\n" unless open(SING, "$files[2]");

	while(my $n = <P1>){
		my $seq = <P1>; <P1>; <P1>;
		$n =~ s/^@/^>/;
		print OUTFA "$n$seq";
	}
	close P1;
	
	while(my $n = <P2>){
		my $seq = <P2>; <P2>; <P2>;
		$n =~ s/^@/^>/;
		print OUTFA "$n$seq";
	}
	close P2;

	while(my $n = <SING>){
		my $seq = <SING>; <SING>; <SING>;
		$n =~ s/^@/^>/;
		print OUTFA "$n$seq";
	}	
	close SING;

	close OUTFA;

}#fastq2fasta

sub map2scpMGs{

	my @tem = @_;

	my $tab = "tab.dimand.out";
	
	if($#tem == 2){
		#map2scpMGs($opt_1, $opt_2,$outcov, $opt_f);
		`$AroCDir/bin/diamond blastx -q $tem[0] -d $dbindex -o $tab -f tab --max-hits 1 -p $opt_n  --max-target-seqs 1 -e $eval`;
		
	
	}elsif($#tem == 3){
		my $temp1 = "temp.1";
		my $temp2 = "temp.2";
		`$AroCDir/bin/diamond blastx -q $tem[0] -d $dbindex -o $temp1 -f tab --max-hits 1 -p $opt_n  --max-target-seqs 1 -e $eval`;
		`$AroCDir/bin/diamond blastx -q $tem[1] -d $dbindex -o $temp2 -f tab --max-hits 1 -p $opt_n  --max-target-seqs 1 -e $eval`;
		`cat $temp1 $temp2 > $tab`;			
		
	}elsif($#tem == 4){
		
		my $temp1 = "temp.1";
		my $temp2 = "temp.2";
		my $temp3 = "temp.3";

		`$AroCDir/bin/diamond blastx -q $tem[0] -d $dbindex -o $temp1 -f tab --max-hits 1 -p $opt_n  --max-target-seqs 1 -e $eval`;
		`$AroCDir/bin/diamond blastx -q $tem[1] -d $dbindex -o $temp2 -f tab --max-hits 1 -p $opt_n  --max-target-seqs 1 -e $eval`;
		`$AroCDir/bin/diamond blastx -q $tem[2] -d $dbindex -o $temp3 -f tab --max-hits 1 -p $opt_n  --max-target-seqs 1 -e $eval`;
	
		`cat $temp1 $temp2 $temp3 > $tab`;
	}

	##get the coverage of each KO

        my %seq2OGs;
	my %seqlen;
        die "$!\n" unless open(OGMAP, "$cogslis");
        while(<OGMAP>){
                chomp;
                my @tem = split(/\t/, $_);
                $seq2OGs{$tem[0]} = $tem[1];
		$seqlen{$tem[0]} = $tem[2];
        }
        close OGMAP;

	die "$!\n" unless open(TAB, "$tab");
	my %koscov;
	while(<TAB>){
		chomp;
		my @tem = split(/\t/, $_);
		if(exists $seq2OGs{$tem[1]}){
			$koscov{$tem[1]} += $tem[3];				
		}else{
			$koscov{$tem[1]} = $tem[3];
		}
	}	
	close TAB;

	die "$!\n" unless open(OUTCOV, ">$tem[-2]");
	print OUTCOV "Name,Average.coverage,Reference.length\n";
	
	for my $koseq (keys %koscov){
		my $ave = $koscov{$koseq} / $seqlen{$koseq};
		print OUTCOV "$koseq,$ave,$seqlen{$koseq}\n";
	}

	close OUTCOV;
}#map2scpMGs


sub get_aroc{

	my ($cogcov, $out16s, $output) = @_;
	
	##For each COGs get the average coverage firstly 
	##standard deviation and mean value across all COGs
	##calculate the number of full length 16S in this metagenomics sample
	##caculate the AROC by 16s/cell 
	my %seq2OGs;
	die "$!\n" unless open(OGMAP, "$cogslis");
	while(<OGMAP>){
		chomp;
		my @tem = split(/\t/, $_);
		$seq2OGs{$tem[0]} = $tem[1];
	}
	close OGMAP;
	
	die "$!\n" unless open(OGCOV, "$cogcov");
	my %ogAb;
	my %ognum;
	<OGCOV>;
	while(<OGCOV>){
		chomp;
		my @tem = split(/\,/, $_);
		if(exists $seq2OGs{$tem[0]}){
			if(exists $ogAb{$seq2OGs{$tem[0]}}){
				$ogAb{$seq2OGs{$tem[0]}} += $tem[1];	 		
				$ognum{$seq2OGs{$tem[0]}} ++;	

			}else{
				$ogAb{$seq2OGs{$tem[0]}} = $tem[1];	 		
				$ognum{$seq2OGs{$tem[0]}} = 1;	
			}
		}else{
			die "Impossible $_\n";
		}
	}
	close OGCOV;

	die "$!\n" unless open(E16S, "$out16s");
	#my %uqseq;
	my $nu16s = 0;
	while(<E16S>){
		#chomp;
		#my @tem = split /\t/;
		#$uqseq{$tem[0]} = 1;
		if(/^>/){
			$nu16s++;
		}
	}
	close E16S;

	#my @total16sreads = keys %uqseq;
	#my $num16sreads = $#total16sreads + 1;	
	#$num16sreads = $num16sreads * 100/1432; 
	my $num16sreads = $nu16s * $opt_s / 1520;
	die "$!" unless open(TOUT, ">$output");
	print TOUT "Estimated_16S_number\t$nu16s\n";

	my $overallcov = 0;
	my $totalOGs = 0;
	for my $ogid (keys %ogAb){
		my $avc = $ogAb{$ogid};	
		$overallcov += $avc;
		$totalOGs ++;
		print TOUT "$ogid\t$avc\t$ognum{$ogid}\n";
	}

	my $finalaroc = $num16sreads * $totalOGs/ $overallcov;
	my $aveko = $overallcov / $totalOGs;
	print TOUT "AverageKO_Estimated_CellNumber\t$aveko\n";
	print TOUT "ARCN\t$finalaroc\n";

}#get_aroc


#-------------------------------------------------------------------------
##This function calculate aroc from otu table format 16S data 
sub otutablearoc{

	my ($tab, $out) = @_;

	die "$!\n" unless open(I, "$tab");
	die "$!\n" unless open(T, ">$out");

	<I>;
	my $samplename = <I>; chomp($samplename);
	my @samples = split(/\t/, $samplename);
	shift @samples; pop @samples;

	my %sotuab;
	my %sotucn;
	my %stotal;
	##inital sample abundance
	for(my $i=0; $i<=$#samples; $i++){
		$stotal{$samples[$i]} = 0;
	}


	while(<I>){
		chomp;
		my @tem = split /\t/;
		my $otuid = shift @tem;
		my $otutax = pop @tem; 
		my $cn = 2; ##Copy number default	
		##average copy number of this otu
			if(exists $copyn{$otutax}){
				$cn = $copyn{$otutax};
			}else{
				my @tttmp = split("; ", $otutax);	
				$cn = copyupper(@tttmp);

				if($cn == 0){
					next;
				}
			}
		for(my $i=0; $i <= $#tem; $i++){
			$sotuab{$samples[$i]}{$otuid} = $tem[$i];
			$sotucn{$samples[$i]}{$otuid} = $cn;	
			$stotal{$samples[$i]} += $tem[$i];
		}

	}
	close I;

	##Output the ACN
	for(my $i=0; $i<=$#samples; $i++){
		my $scn = 0;
		for my $id (keys %{$sotuab{$samples[$i]}}){
			next if($stotal{$samples[$i]} == 0);
			$scn += ($sotuab{$samples[$i]}{$id} / $stotal{$samples[$i]} ) * $sotucn{$samples[$i]}{$id};	
		}
		print T "$samples[$i]\t$scn\n";
	}
}

##sub function copyupper, while there is not a copy number for a certain taxonomy rank, get the ancestor level copy number until there is one, otherwise negnect this OTUs 
sub copyupper {
        my @all = @_;
        my @tem = ("k__", "p__", "c__", "o__", "f__", "g__", "s__");
        my $cop = 1;

        for(my $i = 0; $i <= 6; $i++){
                $tem[$i] = $all[$i];
                my $test = join("; ", @tem);
                if(exists $copyn{$test}){
                        $cop = $copyn{$test};
                }else{
			##If the taxa is not in CopyRighter database, even in domain rank, this OTU is juest ignore
                        return 0;
                }
        }
        #print "$cop\n";
        return $cop;
}

#-------------------------------------------------------------------------------------------------------
#      extract_hvr($opt_U, $$us16s2v6);
#      hyper_community($out16sreads, $community);
#      aroc_16scommunity($community, $opt_o);
#------------------------------------------------------------------------------------------------------
sub extract_hvr{

	my @tem = @_;	
	##using usearch to extract all the potential 16S rRNA V6 sequences from metagenomics dataset.
		
	if($#tem == 1){
		`$AroCDir/bin/usearch -ublast $tem[0]  -db $REFHVR6 -evalue 1e-5 -accel 0.5 -blast6out $tem[1] -threads $opt_n -strand both  -maxaccepts 10 `;
		
	}elsif($#tem == 2){
		my $o1 = "temp1.16sv6";
		my $o2 = "temp2.16sv6";

		`$AroCDir/bin/usearch -ublast $tem[0]  -db $REFHVR6 -evalue 1e-5 -accel 0.5 -blast6out $o1 -threads $opt_n -strand both  -maxaccepts 10 `;
		`$AroCDir/bin/usearch -ublast $tem[1]  -db $REFHVR6 -evalue 1e-5 -accel 0.5 -blast6out $o2 -threads $opt_n -strand both  -maxaccepts 10 `;
		`cat $o1 $o2 > $tem[2]`;

	}elsif($#tem == 3){	
		my $o1 = "temp1.16sv6";
		my $o2 = "temp2.16sv6";
		my $o3 = "temp3.16sv6";

		`$AroCDir/bin/usearch -ublast $tem[0]  -db $REFHVR6 -evalue 1e-5 -accel 0.5 -blast6out $o1 -threads $opt_n -strand both  -maxaccepts 10 `;
		`$AroCDir/bin/usearch -ublast $tem[1]  -db $REFHVR6 -evalue 1e-5 -accel 0.5 -blast6out $o2 -threads $opt_n -strand both  -maxaccepts 10 `;
		`$AroCDir/bin/usearch -ublast $tem[2]  -db $REFHVR6 -evalue 1e-5 -accel 0.5 -blast6out $o3 -threads $opt_n -strand both  -maxaccepts 10 `;

		`cat $o1 $o2 $o3 > $tem[3]`;
	}
}#extract_hvr

sub hyper_community{
	
	##16s blast6out, output community 
	my @tem = @_;	
	#There are possible multi-alignment for one reads aganist the REFHVR database 	
	my %hits;
	my %community;
	die "$!\n" unless open(TEM, "$tem[0]");
	while(<TEM>){
		chomp;
		my @ts = split(/\t/, $_);
		$hits{$ts[0]}{$ts[1]} = $ts[3];
	}
	close TEM;
	
	for my $s (keys %hits){

		my @seqtax;
        	my $averlen = 0;
		my $count = 0;
		for my $tid (keys %{$hits{$s}}){	
			push @seqtax, $id2tax{$tid};
			$averlen += $hits{$s}{$tid};
			$count ++;
		}
		##For abundance
		$averlen = $averlen / $count;
		##for taxonomy LCA
		my $taxs = lca(@seqtax);

		if(exists $community{$taxs}){
			$community{$taxs} += $averlen/66.2;
			#$community{$taxs} += $maxle/66.2;
		}else{
			$community{$taxs} = $averlen/66.2;
		}

	}
	die "$!\n" unless open(TOUT, ">$tem[1]");
	for my $key (sort keys %community){
        	print TOUT "$key\t$community{$key}\n";
	}
	close TOUT;
}##hyper_community

sub lca{
        my @all = @_;
        my @total = ();
        my @common = ();

        ##split each taxonomy by "; " 


        my $i = 0;
        my $j = 0;
        ##store in matrix and transform column to row
        for($i = 0; $i <= $#all; $i++){
                my @tem = split("; ", $all[$i]);
                die "$all[$i]" if($#tem != 6);
                for($j = 0; $j <= $#tem; $j++){
                        $total[$j][$i] = $tem[$j];
                }
        }

        ##compare and generate consistent taxonomy assignment by over 2/3 of all the matched
        for(my $k = 0; $k < $j; $k++){
                my $flag = 0;

                my %temp = ();

                for(my $l = 0; $l < $i; $l++){
                        if(exists $temp{$total[$k][$l]}){
                              $temp{$total[$k][$l]} ++;
                        }else{
                                $temp{$total[$k][$l]} = 1;
                        }
                }

                my @temn = sort {$temp{$b} <=> $temp{$a}} keys %temp;
                if($temp{$temn[0]} >= 2*($k+1)/3){
                        $common[$k] = $temn[0];

                }else{
                        $common[$k] = $blank[$k];
                }
        }

        return join("; ", @common);
}##Using LCA to generate the taxonomy from aligned sequences 

sub aroc_16scommunity{

	##community file / output aroc.txt 
	my @tem = @_;

	my $total = 0;
	my %weight;
	die "$!\n" unless open(I, "$tem[0]");
	while(<I>){
		chomp;
		my @ts = split(/\t/, $_);
		$total += $ts[1];
		if(exists $copyn{$ts[0]}){

			my $acnum = $ts[1] * $copyn{$ts[0]};
			#print T "$ts[0]\t$acnum\n"; 
			$weight{$ts[0]} = $acnum;

		}else{
			my @tttmp = split("; ", $ts[0]);
			my $val = copyupper(@tttmp);
			#my $acnum = $ts[1] * 2.45;
			my $acnum = $ts[1] * $val;
			next if($val == 0);
			$weight{$ts[0]} = $acnum;
		}
	}
	close I;

	my $normacopy =0;
	for my $k (keys %weight){
		$normacopy += ($weight{$k} / $total);
	}
	
	die "$!\n" unless open(T, ">$tem[1]");
	print T "ARCN\t$normacopy\n";
	close T;
}#get aroc	
#------------------------------------------------------------------------------------------------------
__END__
1;
