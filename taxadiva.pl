#!/usr/bin/perl -w
#
# Authors        : Lavanya Rishishwar, Chris Gaby
# Creation Date  : 23rd Aug 2015
# Last Modified  : 26th Dec 2018
# Version        : 0.12.1
my $version = "0.12.1";
#
#############################################################
use strict;
use Getopt::Long;
use threads;
use threads::shared;
#############################################################

#   Pretty usage
my $programHead = "TaxADivA - TAXonomic Assignment and DIVersity Assessment for amplicon reads\n";
my $baseUsage = "           [-o <STRING. output dir and PREFIX to store results. All results will be stored as PREFIX inside the directory PREFIX. Default: input filename>]
           [-d <STRING. Database. Default: db.fasta>] [-t <STRING. Taxonomy file. Default: tax.tsv>]
           [-p <INT. Primer length to be trimmed. Default: no trimming>]
           [-r <INT. Primer length to be trimmed from the RIGHT. Default: no trimming>]
           [-l <INT. Primer length to be trimmed from the LEFT. Default: no trimming>]
           [-k <FLAG. If specified, the trimmed primers will be retained as separate files.  Default: Don't retain trimmed primers.>]
           [-y <FLAG. Performs MED analysis>]
           [-j <INT. Number of threads. Default: 1>]
           [-g <INT. Depth cutoff for considering file. Default: 5000>]
           [-u <STRING. VSEARCH program path. Default: VSEARCH>]
           [--pear <STRING. PEAR paramaters in double quotes.  Default: \"-v 50 -m 450 -n 350 -p 1.0 -j <threads>\". Validity of the arguments not checked.>]
           [--med <STRING. Oligotyping paramaters in double quotes.  Default: \"\". Validity of the arguments not checked.>]
           [--med-metadata <STRING. MED metadata file to be used for the decompose command.  Specified by the -E option in the decompose command.  Default: decompose_map3.tab.>]
           [--keepc4 <FLAG. TaxADivA will not filter out cluster IV sequences>]
		   
           [-h <FLAG. Prints this help>]
           [-v <FLAG.  Prints the current version of the script.>]
           [--version <FLAG.  Prints the current version of the script.>]";
my $usage = "\n=====================================\n$programHead\n=====================================\n$0  [-1 <STRING. forward read file>] [-2 <STRING. reverse read file>]\n$baseUsage\n\n$0  [-s <STRING. file with set of forward and reverse files>]\n$baseUsage \n=====================================\nExample usage: $0 -d db.fasta -t tax.db.tsv -j 18 -s list.txt -o output1 -m \"-v 50 -m 450 -n 350 -p 1.0 -j 4\"\n";

my $shortUsage = "\nRun $0 -h to display the full help page.\n";

#############################################################
#	Global Variables
my $r1;
my $r2;
my $l;
my $help = 0;
my $b = "beta.tsv";
my $outDir = "out"; # Output file
my $db      = "db.fasta"; # Database
my $taxFile = "tax.tsv"; # Database
my $threads = 10;
my $pLen;
my $prLen;
my $plLen;
my $keepTrimmedPrimers = 0;
my $vsearchProgram = "vsearch";
my $pearParameters = "-v 50 -m 450 -n 350 -p 1.0 -j $threads";
my $medParameters = "";
my $goodDepth = 5000;
my $oligotyping = 0;
my $medMetadataFile = "decompose_map3.tab";
my $maxHsps = "";
my $versionPrint = 0;
my $cluster4filter = 0;
my %tax; my %tax2Species; my %family; my %genus; my %class; my %uniqueClass; my %taxid; my %spfFormat;
#############################################################

#	Argument Check
my $args  = GetOptions ("1=s"     => \$r1,
                        "2=s"     => \$r2,
                        "s=s"     => \$l,
                        "b=s"     => \$b,
                        "h"       => \$help,
                        "d=s"     => \$db,
                        "u=s"     => \$vsearchProgram,
                        "t=s"     => \$taxFile,
                        "j=i"     => \$threads,
                        "p=i"     => \$pLen,
                        "r=i"     => \$prLen,
                        "l=i"     => \$plLen,
                        "k+"      => \$keepTrimmedPrimers,
                        "g=i"     => \$goodDepth,
                        "y+"      => \$oligotyping,
                        "o=s"     => \$outDir,
                        "pear=s"  => \$pearParameters,
                        "med=s"   => \$medParameters,
                        "keepc4+"   => \$cluster4filter,
                        "med-metadata=s"   => \$medMetadataFile,
                        "version+"   => \$versionPrint);
                        
#############################################################

	
if($help == 1){
	print STDERR $programHead;
	print STDERR $usage;
	exit;
}

if($versionPrint > 0){
	print STDERR "Current version of the program TaxADivA: $version\n";
	exit;
}
# Placed here to speed up the help or version print.
share(%tax); share(%tax2Species); share(%family); share(%genus); share(%class); share(%uniqueClass); share(%taxid);
#############################################################

print STDERR "Checking for the provided arguments...";
                        
if(! defined $l && !( defined $r1 && defined $r2)){
	print STDERR "\nERROR (Line ".__LINE__."): list file (option -s) or forward (option -1) and reverse (option -2) reads are mandatory options.  Please specify them!\n";
	print STDERR "$shortUsage\n";
	exit;
}
if(defined $l){
	if(! -e $l){
		print STDERR "\nERROR (Line ".__LINE__."): Please make sure that list file (option -s) exist or the location provided is correct!\n";
		print STDERR "$shortUsage\n";
		exit;
	}
} else {
	if(! -e $r1 || ! -e $r2){
		print STDERR "\nERROR (Line ".__LINE__."): Please make sure that the forward (option -1) and reverse (option -2) reads exist or the location provided is correct!\n";
		print STDERR "$shortUsage\n";
		exit;
	}
}
if(! -e $taxFile){
	print STDERR "\nERROR (Line ".__LINE__."): Please make sure that the taxonomy (option -t) file exist or the location provided is correct!\n";
	print STDERR "$shortUsage\n";
	exit;
}
if(! -e $db){
	print STDERR "\nERROR (Line ".__LINE__."): Please make sure that the database (option -d) file exist or the location provided is correct!\n";
	print STDERR "$shortUsage\n";
	exit;
}
if( ! -e "$db.nhr" || ! -e "$db.nin" || ! -e "$db.nsq" ){
	`makeblastdb -in $db -out $db -dbtype nucl`;
}

if(length($outDir) == 0 || $outDir =~ m{[\\:*?"<>|]}){
	print STDERR "\nERROR (Line ".__LINE__."): Invalid output prefix.  Please provide a prefix with acceptable characters (alphanumeric, _, .)\n";
	print STDERR "$shortUsage\n";
	exit;
}

if ($oligotyping > 0 && ! -e $medMetadataFile){
	print STDERR "\nERROR (Line ".__LINE__."): Please make sure that MED metadata (option --med-metadata) file exist or the location provided is correct!\n";
	print STDERR "$shortUsage\n";
	exit;
}

print STDERR "everything looks fine.\n";


if(-e $outDir && -d $outDir){
	print STDERR "\nWARNING: Output directory already exists! Content will be overwritten!\n";
}
`mkdir -p $outDir`;


#############################################################
#	Prerequisite check

print STDERR "Checking for dependencies...";

my $programPath;
# PEAR
$programPath = `which pear`;
die "\nERROR (Line ".__LINE__."): PEAR not found! Please make sure PEAR is installed and is in the \$PATH variable.\nTo install PEAR, download the binaries/source from http://sco.h-its.org/exelixis/web/software/pear/ and place it in your \$PATH\n" if(length($programPath) == 0);
$programPath = "";

# BLAST
$programPath = `which blastn`;
die "\nERROR (Line ".__LINE__."): BLASTN not found! Please make sure BLASTN is installed and is in the \$PATH variable.\nTo install BLASTN, download the binaries/source from http://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download and place it in your \$PATH\n" if(length($programPath) == 0);
# Not used for now
#my $blastVersion = `blastn -version | head -1 | sed 's/blastn: //'`;
#chomp $blastVersion;
#if($blastVersion =~ m/^2.2/){
	# The BLAST version is 2.2.XX
#	$blastVersion =~ s/2.2.//;
#	$blastVersion =~ s/\+$//;
#	if($blastVersion < 31){
#		die "An older version of BLAST is detected (2.2.$blastVersion).  Please upgrade to BLAST+ version 2.2.31+.  The BLAST latest binaries can be obtained from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/\n";
#	}
#	$maxHsps = "-max_hsps 1";
#} elsif ($blastVersion =~ m/^2.3/) {
#	$maxHsps = "-max_hsps 1";
#} elsif ($blastVersion =~ m/^2.4/){
#	$maxHsps = "-max_hsps 1";
#} else {
#	print STDERR "WARNING: BLAST version of $blastVersion was found, this has not been tested with the script and may throw error.  The script has been tested on 2.2.31+, 2.3 or 2.4 versions.  Proceeding nevertheless.\n";
#}
$programPath = "";

# VSEARCH
$programPath = `which $vsearchProgram`;
die "\nERROR (Line ".__LINE__."): VSEARCH not found! Please make sure VSEARCH is installed and is in the \$PATH variable.\nTo install VSEARCH, download the source from https://github.com/torognes/vsearch and place it in your \$PATH\n" if(length($programPath) == 0);
$programPath = "";

# KRONA
$programPath = `which ktImportText`;
die "\nERROR (Line ".__LINE__."): Krona not found! Please make sure Krona is installed and is in the \$PATH variable.\nTo install Krona, download the binaries/source from https://github.com/marbl/Krona/wiki and place it in your \$PATH\n" if(length($programPath) == 0);
$programPath = "";



# Optional program check # 1
if(defined $pLen || defined $prLen || defined $plLen){
	# Prinseq-lite.pl
	$programPath = `which prinseq-lite.pl`;
	die "\nERROR (Line ".__LINE__."): Prinseq-lite not found! Please make sure Prinseq-lite is installed and is in the \$PATH variable.\nTo install prinseq-lite, download the source from http://prinseq.sourceforge.net/ and place it in your \$PATH\n" if(length($programPath) == 0);
	$programPath = "";
}
if(defined $pLen){
	$prLen = $pLen;
	$plLen = $pLen;
} else {
	$prLen //= 0;
	$plLen //= 0;
}

# Optional program check # 2
if($oligotyping == 1){
	# Minimum Entropy Decomposition (MED)
	$programPath = `which decompose`;
	die "\nERROR (Line ".__LINE__."): MED not found! Please make sure MED is installed and is in the \$PATH variable.\nTo install MED, refer to the instructions provided here: http://merenlab.org/projects/oligotyping/\n" if(length($programPath) == 0);
	$programPath = "";
	
	$programPath = `which o-pad-with-gaps`;
	die "\nERROR (Line ".__LINE__."): o-pad-with-gaps not found! Please make sure oligotype is installed and is in the \$PATH variable.\nTo install oligotyping, refer to the instructions provided here: http://merenlab.org/projects/oligotyping/\n" if(length($programPath) == 0);
	$programPath = "";
	
	$programPath = `which entropy-analysis`;
	die "\nERROR (Line ".__LINE__."): entropy-analysis not found! Please make sure oligotype is installed and is in the \$PATH variable.\nTo install oligotyping, refer to the instructions provided here: http://merenlab.org/projects/oligotyping/\n" if(length($programPath) == 0);
	$programPath = "";
}

print STDERR "found everything I need.\n";


#############################################################
#	Function definitions
#############################################################

# Function: acceptable#####################
#	Input: text
#   Small function that filters out things we do not want to go on the output
#	Output: filtered result
###########################################
sub acceptable{
	if($_[0] =~ m/(unclassified)|(unidentified)|(enrichment)|(miscellaneous)|(marine)|(environmental)|(^\s*\.\s*$)/){
		return 0;
	} else {
		return 1;
	}
}

# Function: readTaxonomy####################
#	Input: none
#
#	Output: none; initializes hashes
###########################################
sub readTaxonomy{
	print STDERR "Commencing preprocessing...";
	open TAX, "<$taxFile" or die "Cannot open input file $taxFile: $!\n";
	while(<TAX>){
		chomp $_;
		my ($accn, $cluster, $taxid, $tax) = split(/\t/, $_);
		my ($kingdom, $phylum, $class, $order, $family, $genus, $species, $strain) = split(/;/, $tax);
		$taxid{$accn} = $taxid;
		
		$tax =~ s/;/\t/g;
		
		die "\nERROR (Line ".__LINE__."): Issue at stage 1 with $_\nNew tax line is $tax\n" if($tax =~ m/^\s*Bacteria\s*$/);
		$tax{$accn} = "$cluster\t$tax";
		
		if(defined $genus && length($genus) > 2){
			$tax =~ s/\t[^\t]+?$//;
		}
		die "\nERROR (Line ".__LINE__."): Issue at stage 2 with $_\nNew tax line is $tax\n" if($tax =~ m/^\s*Bacteria\s*$/);
		$genus{$accn} = "$cluster\t$tax";
		
		if(defined $family && length($family) > 2){
			$tax =~ s/\t[^\t]+$//;
		}
		die "\nERROR (Line ".__LINE__."): Issue at stage 3 with $_\nNew tax line is $tax\n" if($tax =~ m/^\s*Bacteria\s*$/);
		#$tax = join("\t", $kingdom, $phylum, $class, $order, $family);
		$family{$accn} = "$cluster\t$tax";
		
		$class =~ s/^\s+//;
		$class =~ s/\s+$//;
		
		$class{$accn} = $class;
		$uniqueClass{$class} = 1;
		
		$phylum = "unclassified" unless(acceptable($phylum));
		$class = "unclassified" unless(acceptable($class));
		$order = "unclassified" unless(acceptable($order));
		$family = "unclassified" unless(acceptable($family));
		$genus = "unclassified" unless(acceptable($genus));
		$species = "unclassified" unless(acceptable($species));
		
		$spfFormat{$accn}{"full"} = "k__".$kingdom."\tp__".$phylum."\tc__".$class."\to__".$order."\tf__".$family."\tg__".$genus."\ts__".$species;
		$spfFormat{$accn}{"kingdom"} = "k__".$kingdom;
		$spfFormat{$accn}{"phylum"} = "p__".$phylum;
		$spfFormat{$accn}{"class"} = "c__".$class;
		$spfFormat{$accn}{"order"} = "o__".$order;
		$spfFormat{$accn}{"family"} = "f__".$family;
		$spfFormat{$accn}{"genus"} = "g__".$genus;
		$spfFormat{$accn}{"species"} = "s__".$species;
	}
	close TAX;
	print STDERR "Done\n";
}
# Function: preprocess ####################
#	Input: read1, read2, output prefix 
#		and primer length to be trimmed off
#
#	Output: 0 if preprocessing failed, 1 if
#			passed
#
#	This will generate $out.reads.temp
###########################################
sub preprocess{
	my ($r1, $r2, $out, $pRight, $pLeft) = @_;	
	print STDERR "[$out][Step 0]\tLooking if the input file has sufficient depth for the analysis...\n";
	
	my $count = `wc -l $r1 | awk '{print \$1}'`;
	$count /= 4;
	
	if($count < $goodDepth){
		print STDERR "[$out][Step 0]\tSorry, less than $goodDepth reads found for $out. Skipping this file\n";
		return 0;
	}
	
	print STDERR "[$out][Step 0]\t$out passed the criteria\n";
	
	print STDERR "[$out][Step 1]\tMerging reads using PEAR (approx ~ 2 min)...\n";
	
	`pear -f $r1 -r $r2 -o $outDir/$out-pear $pearParameters`; # Need to revisit this command in future
	
	print STDERR "[$out][Step 1]\tConverting FASTQ to FASTA...\n";
	`awk 'BEGIN{c=0} {c++; if(c==1){print ">"\$0}; if(c==2){print}; if(c==4){c=0}}' $outDir/$out-pear.assembled.fastq > $outDir/$out-pear.fa`;
	
	if($pRight+$pLeft > 0){
		print STDERR "[$out][Step 1]\tTrimming primers...\n";
		`prinseq-lite.pl -fasta $outDir/$out-pear.fa -out_good $outDir/$out -trim_left $pLeft -trim_right $pRight 1>&2 2>> $outDir/$out.log`;
		if($keepTrimmedPrimers > 0){
			open FILE, "<$outDir/$out-pear.fa" or die "Cannot open PEAR output $outDir/$out-pear.fa: $!\n";
			open LEFT, ">$outDir/$out.trimmedLeft.fa" or die "Cannot create output file $outDir/$out.trimmedLeft.fa: $!\n";
			open RIGHT, ">$outDir/$out.trimmedRight.fa" or die "Cannot create output file $outDir/$out.trimmedRight.fa: $!\n";
			
			while(<FILE>){
				my $desc = $_;
				my $seq = <FILE>;
				chomp $seq;
				my $left = substr($seq, 0, $pLeft);
				my $right = substr($seq, -1*$pRight);
				print LEFT "$desc\n$left\n";
				print RIGHT "$desc\n$right\n";
			}
			close FILE;
		}
		`mv $outDir/$out.fasta $outDir/$out.reads.temp`;
		print STDERR "[$out][Step 1]\tTrimming Done\n";
	} else {
		`mv $outDir/$out-pear.fa $outDir/$out.reads.temp`;
	}
	
	`rm $outDir/$out-pear*`;
	
	return 1;
}

# Function: mergeReads ####################
#	Input: merged files that passed 
#			preprocessing
#
#	Output: Single merged file name
#
#	This will generate $out.read.temp
###########################################
sub mergeReads{
	my ($prefix, @files) = @_;
	my $out = $prefix.".reads.temp";
	print STDERR "[$prefix][Step 1]\tMerging processed reads (total of ".@files." files)...";
	open OUT, ">$outDir/$out" or die "\nERROR (Line ".__LINE__."): Cannot create output file $outDir/$out:$!\n";
	
	foreach my $file (@files){
		
		open FILE, "<$outDir/$file.reads.temp" or die "\nERROR (Line ".__LINE__."): Cannot open input file $outDir/$file.reads.temp:$!\n";
		my $prefix = $file;
		$prefix =~ s/.reads.temp//;
		
		while(<FILE>){
			if($_ =~ m/^>/){
				$_ =~ s/^>/>$prefix:/;
			}
			print OUT $_;
		}
		
		close FILE;
	}
	close OUT;
	print STDERR "done\n";
}

# Function: taxAssign #####################
#	Input: File prefix to be processed,
#	blast db prefix
#
#	Output: None
#
#	This generates the following files:
#	* {prefix}.reads.fasta - chimera free merged reads
#	* {prefix}.cent.fa - centroids of the clusters
#	* {prefix}.clusters.tsv - tabular clusters files
#	* {prefix}.blast.tsv - blast output file
#	* {prefix}.alpha.txt - alpha diversity file
#	* {prefix}.abundance.txt - class abundance file
#	* {prefix}.krona.txt - krona input file
###########################################
sub taxAssign{
	my ($out, $db) = @_;
	
	# Chimera removal and clustering
	print STDERR "[$out][Step 2]\tDereplicating sequences...";
	`$vsearchProgram -derep_fulllength $outDir/$out.reads.temp --output $outDir/$out.reads.derep -sizeout -threads $threads 1>> $outDir/$out.log 2>> $outDir/$out.log`;
	#`rm $out.reads.temp`;
	
	print STDERR "Done\n[$out][Step 3]\tClustering sequences...";
	`$vsearchProgram -cluster_fast $outDir/$out.reads.derep -centroids $outDir/$out.otus.fa -uc $outDir/$out.up -minsize 2 --id 0.9 -strand both --sizeout 1>> $outDir/$out.log 2>> $outDir/$out.log`;
	
	print STDERR "Done\n[$out][Step 3]\tRemoving chimeras...";
	`$vsearchProgram -uchime_ref $outDir/$out.otus.fa -db $db -strand plus -minh 1.0 -nonchimeras $outDir/$out.nochim.fa -uchimeout $outDir/$out.uchime -uchimealns $outDir/$out.aln 1>> $outDir/$out.log 2>> $outDir/$out.log`;
	

	my $total = 0;
	my %otuSizes;
	my %otu2Keep;
	my $otu2Keep = `grep '^>' $outDir/$out.nochim.fa`;
	my @otu2Keep = split(/\n/, $otu2Keep);
	foreach(@otu2Keep){
		chomp $_;
		$_ =~ s/^>//;
		$otu2Keep{$_} = 1;
	}
	
	open FILE, "<$outDir/$out.up" or die "\nERROR (Line ".__LINE__."): Cannot open input file $outDir/$out.clusters.tsv: $!\n";
	while(<FILE>){
		chomp $_;
		my ($type, @vals) = split(/\t/, $_);
		next if($type ne "S");
		
		pop(@vals);
		my $otu = pop(@vals);
		
		next if(! defined $otu2Keep{$otu}); #Throw out chimeric OTUs
		
		my $size = $otu;
		$size =~ s/.*;size=//;
		$size =~ s/;\s*$//;
		
		$total += $size;
		$otuSizes{$otu} += $size;
	}
	close FILE;
	
	# Perform the BLAST search
	print STDERR "Done\n[$out][Step 4]\tLooking for homologs...";
	my $blOut = `blastn -query $outDir/$out.nochim.fa -db $db -outfmt "6 qseqid sseqid evalue pident" -evalue 0.01 -max_target_seqs 1  -max_hsps 1 -num_threads $threads | tee $outDir/$out.blast.tsv`;
	print STDERR "Done\n[$out][Step 5]\tProcessing results...";
	
	# These are emperically calculated threshold values
	my $family  = 75;
	my $genus   = 88.1;
	my $species = 91.9;
	
	# These variables will store the processed results
	my %counts; # stores the counts of the assignment
	my %classCounts; # this guy stores how many times we have seen a class
	#my %alpha; # stores each species name for alpha diversity # Not used
	my $assigned = 0; # this guy stores how many assignments have been made
	
	
	my @blast = split(/\n/, $blOut);
	foreach (@blast){
		my ($q, $s, $e, $p) = split(/\t/, $_);
		
		if($p > $family){
			my ($accn, undef) = split(/;/, $s);
			my $size = $otuSizes{$q};
			
			# Throw away paralogs and account for that change in the total
			if($cluster4filter == 0 && $s =~ /cluster_IV/){
				next;
			}
			#$total += $size;
			
			# If the sequence has passed the family threshold, than it can be assigned to the same class
			$classCounts{$class{$accn}} += $size;
			$assigned += $size;
			
			if($p > $species){ # Species cutoff
				die "\nERROR (Line ".__LINE__."): Cannot file $accn from $q!exiting\n" if(! defined $tax{$accn}); # failsafe
				$counts{$tax{$accn}} += $size;
				#$alpha{$tax2Species{$accn}}  += $size;
			} elsif($p > $genus) { # Genus cutoff
				die "\nERROR (Line ".__LINE__."): Cannot file $accn from $q!exiting\n" if(! defined $genus{$accn}); # failsafe
				$counts{$genus{$accn}} += $size;
			} else { # Everything else
				die "\nERROR (Line ".__LINE__."): Cannot file $accn from $q!exiting\n" if(! defined $family{$accn}); # failsafe
				$counts{$family{$accn}} += $size;
			}
		}
	}
	
	print STDERR "Done\n[$out][Step 6]\tOutputting results to files...";
	
	
#	#This files stores the alpha diversity
#	open OUT, ">$out.alpha.txt" or die "\nERROR (Line ".__LINE__."): Cannot create output file $out.alpha.txt: $!\n";
#	print OUT "Species\tCounts\n";
#	foreach my $key (keys %alpha){
#		next if($key eq "NA");
#		print OUT "$key\t$alpha{$key}\n";
#	}
#	close OUT;
	
	# Class abundance plot file
	open OUT, ">$outDir/$out.abundance.txt" or die "\nERROR (Line ".__LINE__."): Cannot create output file $outDir/$out.abundance.txt: $!\n";
	my $percent = 0;
	my @classes = sort keys %uniqueClass;
	print OUT "Class\t$out\n";
	foreach my $key (@classes){
		$classCounts{$key} //= 0;
		my $this = 100*$classCounts{$key}/$total;
		print OUT "$key\t$this\n";
		$percent += $this;
	}
	print OUT "Unknown OTUs\t".(100-$percent)."\n";
	close OUT;
	
	# This file will be used to generate the krona plot
	open OUT, ">$outDir/$out.krona.txt" or die "\nERROR (Line ".__LINE__."): Cannot create output file $outDir/$out.krona.txt: $!\n";
	foreach my $key (keys %counts){
		my $text = $key;
		$text =~ s/\//&#47;/g;
		print OUT "$counts{$key}\t$text\n";
	}
	print OUT ($total-$assigned)."\tUnassigned OTU\n";
	close OUT;
	print STDERR "Done\n\tMaking plots...";
	`ktImportText -n nifH $outDir/$out.krona.txt -o $outDir/$out.krona.html`;
}
# Function: taxAssignBatch ################
#	Input: File prefix to be processed,
#	blast db prefix
#
#	Output: None
#
#	This generates the following files:
#	* {prefix}.reads.fasta - chimera free merged reads
#	* {prefix}.cent.fa - centroids of the clusters
#	* {prefix}.clusters.tsv - tabular clusters files
#	* {prefix}.blast.tsv - blast output file
#	* {prefix}.alpha.txt - alpha diversity file
#	* {prefix}.abundance.txt - class abundance file
#	* {prefix}.krona.txt - krona input file
#	* {prefix}.otu_table.txt - OTU table that can be converted to BIOM format
###########################################
sub taxAssignBatch{
	my ($out, $db) = @_;
	
	# Chimera removal and clustering
	print STDERR "[$out][Step 2]\tDereplicating sequences...";
	`$vsearchProgram -derep_fulllength $outDir/$out.reads.temp --output $outDir/$out.reads.derep -sizeout -threads $threads 1>> $outDir/$out.log 2>> $outDir/$out.log`;
	#`rm $out.reads.temp`;
	
	print STDERR "Done\n[$out][Step 3]\tClustering sequences...";
	`$vsearchProgram -cluster_fast $outDir/$out.reads.derep -centroids $outDir/$out.otus.fa -uc $outDir/$out.up -minsize 2 --id 0.9 -strand both --sizeout 1>> $outDir/$out.log 2>> $outDir/$out.log`;
	
	print STDERR "Done\n[$out][Step 3]\tRemoving chimeras...";
	`$vsearchProgram -uchime_ref $outDir/$out.otus.fa -db $db -strand plus -minh 1.0 -nonchimeras $outDir/$out.nochim.fa -uchimeout $outDir/$out.uchime -uchimealns $outDir/$out.aln 1>> $outDir/$out.log 2>> $outDir/$out.log`;
	
	# Perform the BLAST search
	print STDERR "Done\n[$out][Step 4]\tLooking for homologs...";
	my $blOut = `blastn -query $outDir/$out.nochim.fa -db $db -outfmt "6 qseqid sseqid evalue pident" -evalue 0.01 -max_target_seqs 1 -max_hsps 1 -num_threads $threads | tee $outDir/$out.blast.tsv`;
	print STDERR "Done\n[$out][Step 5]\tProcessing results...";
	
	# These are emperically calculated threshold values
	my $family  = 75;
	my $genus   = 88.1;
	my $species = 91.9;
	
	my $totalSeqs = `grep -c '>' $outDir/$out.nochim.fa`;
	chomp $totalSeqs;
	
	# These variables will store the processed results
	my %assignments; # stores mapping of centroids to taxonomy
	
	my @blast = split(/\n/, $blOut);
	my $seqsProcessed = 0;
	
	
	## 
	## This piece of code stores the taxonomy assignment of
	## each centroid
	##
	for(my $i=0; $i < @blast; $i++){
		my ($centroid, $s, $e, $p) = split(/\t/, $blast[$i]);
		
	
		$centroid =~ s/;size.*$//;
		if($cluster4filter == 0 && $s =~ /cluster_IV/){
			$assignments{$centroid}{"class"} = "cluster_IV";
			next;
		}
		my ($accn, undef) = split(/;/, $s);
		
		$assignments{$centroid}{"stamp"} = join("\t", $spfFormat{$accn}{"kingdom"}, $spfFormat{$accn}{"phylum"}, $spfFormat{$accn}{"class"}, $spfFormat{$accn}{"order"});
		$assignments{$centroid}{"qiime"} = join("; ", $spfFormat{$accn}{"kingdom"}, $spfFormat{$accn}{"phylum"}, $spfFormat{$accn}{"class"}, $spfFormat{$accn}{"order"});
		$seqsProcessed++;
		
		if($p > $family){
			#print "$centroid assigned $class{$accn}\n";
			$assignments{$centroid}{"class"} = $class{$accn};
			$assignments{$centroid}{"accn"} = $accn;
			$assignments{$centroid}{"stamp"} .= "\t".$spfFormat{$accn}{"family"};
			$assignments{$centroid}{"qiime"} .= "; ".$spfFormat{$accn}{"family"};
			
			if($p > $species){
				die "\nERROR (Line ".__LINE__."): Cannot find $accn from $centroid! Do all the sequences in the database exist in the taxonomy file (or is the database build an different version from taxonomy file)? Exiting\n" if(! defined $tax{$accn});
				$assignments{$centroid}{"assign"} = $tax{$accn};
				#$assignments{$centroid}{"alpha"} = $tax2Species{$accn};
				$assignments{$centroid}{"stamp"} .= "\t".$spfFormat{$accn}{"genus"}."\t".$spfFormat{$accn}{"species"};
				$assignments{$centroid}{"qiime"} .= "; ".$spfFormat{$accn}{"genus"}."; ".$spfFormat{$accn}{"species"};
				
			} elsif($p > $genus) {
				die "\nERROR (Line ".__LINE__."): Cannot find $accn from $centroid! Do all the sequences in the database exist in the taxonomy file (or is the database build an different version from taxonomy file)? Exiting\n" if(! defined $genus{$accn});
				$assignments{$centroid}{"assign"} = $genus{$accn};
				$assignments{$centroid}{"stamp"} .= "\t".$spfFormat{$accn}{"genus"}."\tUnclassified";
				$assignments{$centroid}{"qiime"} .= ";".$spfFormat{$accn}{"genus"}."; s__";
			} else {
				die "\nERROR (Line ".__LINE__."): Cannot find $accn from $centroid! Do all the sequences in the database exist in the taxonomy file (or is the database build an different version from taxonomy file)? Exiting\n" if(! defined $family{$accn});
				$assignments{$centroid}{"assign"} = $family{$accn};
				$assignments{$centroid}{"stamp"} .= "\tUnclassified\tUnclassified";
				$assignments{$centroid}{"qiime"} .= "; g__; s__";
			}
		} else {
			$assignments{$centroid}{"stamp"} .= "\tUnclassified\tUnclassified\tUnclassified";
			$assignments{$centroid}{"qiime"} .= "; f__; g__; s__";
		}
		
		# print "Looking at $centroid which matches $s at $p -> assigned to ".$assignments{$centroid}{"qiime"}."\n";
	}
	
	# These variables will store the processed results
	my %clusterCount; # stores the counts of the assignment
	my %classCounts; # this guy stores how many times we have seen a class
	my %assigned; # this guy stores how many assignments have been made
	my %total; # this guy stores how many assignments have been made
	my %otus;
	my $otuNum = 1;
	
	# These variables only store things for printing purposes
	my %classesSeen;
	
	# This variable stores the count for creating a biom file
	# bins are the samples
	
	my %biom;
	my %biomSamples;
	
	my %otu2Keep;
	my $otu2Keep = `grep '^>' $outDir/$out.nochim.fa`;
	my @otu2Keep = split(/\n/, $otu2Keep);
	foreach(@otu2Keep){
		chomp $_;
		$_ =~ s/^>//;
		$_ =~ s/;size.*//;
		$otu2Keep{$_} = 1;
	}
	
	my $unassignedSeqs = 0;
	open FILE, "<$outDir/$out.up" or die "\nERROR (Line ".__LINE__."): Cannot open input file $outDir/$out.up: $!\n";
	while(<FILE>){
		chomp $_;
		my ($type, undef, undef, undef, undef, undef, undef, undef, $label1, $label2) = split(/\t/, $_);
		# columns are: (1) record type, (2) cluster number, (3) centroid length / cluster size, 
		#              (4) % identity, (5) match orientation, (6) *, (7) *, 
		#              (8) alignment in CIGAR, (9) label of query (H) or centroid (S), 
		#              (10) label of centroid (H) or * (S)
		next if($type eq "C");
		
		# $otu is the centroid label, which is used to figure out
		# what taxonomic assignment was made
		my $otu = "";
		my $sample = "";
		my $size = 0; # this guy contains dereplication information
		if($type eq "S"){
			$otu = $label1;
			($sample, undef) = split(/:/, $otu);
			
			$size = $otu;
		} else {
			$otu = $label2;
			($sample, undef) = split(/:/, $label1);
			
			$size = $label1;
		}
		$otu =~ s/;size=.*$//;
		$size =~ s/.*;size=//;
		$size =~ s/;\s*$//;
		
		next if(! defined $otu2Keep{$otu}); #Throw out chimeric OTUs
		
		if(defined $assignments{$otu}{"class"}){
			next if($cluster4filter == 0 && $assignments{$otu}{"class"} eq "cluster_IV"); # Precautionary check, should never be true
			
			if(defined $assignments{$otu}{"class"}){
				$classCounts{$sample}{$assignments{$otu}{"class"}} += $size;
				$classesSeen{$assignments{$otu}{"class"}} = 1;
			}
			
			$assigned{$sample} += $size;
			
		} else {
			$assignments{$otu}{"assign"} = "OTU$otuNum";
			$otuNum++;
		}
		
		$biom{$otu}{$sample} += $size;
		$biomSamples{$sample} = 1;
		
		$clusterCount{$sample}{$assignments{$otu}{"assign"}} += $size;
		$total{$sample} += $size;
	}
	close FILE;
	
	print STDERR "Done\n[$out][Step 6]\tOutputting results to files...";
	
	# Class relative abundance files
	open ABND, ">$outDir/$out.abundance.txt" or die "\nERROR (Line ".__LINE__."): Cannot create output file $outDir/$out.abundance.txt: $!\n";
	foreach my $key (keys %classesSeen){
		next unless(acceptable($key));
		print ABND "\t$key";
	}
	print ABND "\tUnassignedOTUs\n";
	foreach my $sample (keys %classCounts){
		print ABND "$sample";
		my $percent = 0;
		foreach my $tax (keys %classesSeen){
			next unless(acceptable($tax));
			$classCounts{$sample}{$tax} //= 0;
			my $this = $classCounts{$sample}{$tax}*100/$total{$sample};
			print ABND "\t$this";
			$percent += $this;
		}
		print ABND "\t".(100-$percent)."\n";
	}
	close ABND;
	print STDERR "\n\tCreated $out.abundance.txt";
		
	# Set of krona files
	foreach my $sample (keys %clusterCount){
		open KRONA, ">$outDir/$sample.krona.txt" or die "\nERROR (Line ".__LINE__."): Cannot create output file $outDir/$sample.krona.txt: $!\n";
		foreach my $key (keys %{$clusterCount{$sample}}){
			my $text = $key;
			$text =~ s/\//&#47;/g;
			next if($text =~ m/^OTU/);
			print KRONA "$clusterCount{$sample}{$key}\t$text\n";
		}
		print KRONA ($total{$sample}-$assigned{$sample})."\tUnassigned OTU\n";
		close KRONA;
		`ktImportText -n nifH $outDir/$sample.krona.txt -o $outDir/$sample.krona.html`;
	}
	print STDERR "\n\tCreated Krona files (Total files = ".(scalar keys %clusterCount).")";
	
	
	# This is where the BIOM file is generated for input for QIIME
	
	my @biomSamples = keys %biomSamples;
	my @otus = keys %biom;
	
	open BIOM, ">$outDir/$out.otu_table.txt" or die "\nERROR (Line ".__LINE__."): Cannot create output file $outDir/$out.otu_table.txt: $!\n";
	print BIOM "# Constructed from biom file\n";
	print BIOM "#OTU ID";
	foreach (@biomSamples){
		print BIOM "\t$_";
	}
	print BIOM "\ttaxonomy\n";
	for(my $i = 0; $i < @otus; $i++){
		
		print BIOM "denovo".($i+1);
		foreach my $sample (@biomSamples){
			$biom{$otus[$i]}{$sample} //= 0;
			print BIOM "\t$biom{$otus[$i]}{$sample}";
		}
		if(! defined $assignments{$otus[$i]}{"qiime"}){
			print BIOM "\tk__; p__; c__; o__; f__; g__; s__\n";
		} else {
			#The following code would've printed the taxonomy ID
			#print BIOM "\t".$taxid{$assignments{$otus[$i]}{"accn"}} if (defined $assignments{$otus[$i]}{"accn"});
			#The new one will print the taxonomy hierarchy
			# if(defined $spfFormat{$assignments{$otus[$i]}{"accn"}}{"full"}){
				# my $name = $spfFormat{$assignments{$otus[$i]}{"accn"}}{"full"};
				# $name =~ s/\t/;/g;
				# print BIOM "\t$name";
			# }
			# print BIOM "\n";
			print BIOM "\t".$assignments{$otus[$i]}{"qiime"}."\n";
		}
	}
	close BIOM;
	print STDERR "\n\tCreated $out.otu_table.txt";
	
	
	# This is where the STAMPS spf file is generated
	
	open SPF, ">$outDir/$out.spf" or die "\nERROR (Line ".__LINE__."): Cannot create output file $outDir/$out.spf: $!\n";
	print SPF "Domain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies";
	foreach (@biomSamples){
		print SPF "\t$_";
	}
	print SPF "\n";
	my %unclass;
	my $unclassTag = "";
	my %spfTag;
	my %spfBinCounts;
	for(my $i = 0; $i < @otus; $i++){
		if(! defined $assignments{$otus[$i]}{"stamp"}){
			foreach my $sample (@biomSamples){
				$biom{$otus[$i]}{$sample} //= 0;
				$unclass{$sample} += $biom{$otus[$i]}{$sample};
			}
		} else {
			$spfTag{$assignments{$otus[$i]}{"stamp"}} = 1;
			foreach my $sample (@biomSamples){
				$biom{$otus[$i]}{$sample} //= 0;
				$spfBinCounts{$assignments{$otus[$i]}{"stamp"}}{$sample} += $biom{$otus[$i]}{$sample};
			}
		}
	}
	foreach my $tag (keys %spfTag){
		print SPF $tag;
		foreach my $sample (@biomSamples){
			print SPF "\t$spfBinCounts{$tag}{$sample}";
		}
		print SPF "\n";
	}
	print SPF "Unclassified\tUnclassified\tUnclassified\tUnclassified\tUnclassified\tUnclassified\tUnclassified"; 
	foreach my $sample (@biomSamples){
		print SPF "\t$unclass{$sample}";
	}
	print SPF "\n";
	close SPF;
	print STDERR "\n\tCreated $out.spf";
	
	print STDERR "\n\tDone creating files!\n";
}

# Function: oligotype #####################
#	Input: merged files from mergeReads, 
#	metadata file
#
#	Output: MED files
#
#	Will also re-label the NETWORK.gexf and TOPOLOGY.gefx
#	based on the labelling from the database and the 
#	representative sequences from NODE-REPRESENTATIVES.fasta
#	file 
#
###########################################
sub oligotype{
	my ($prefix) = @_;
	
	print STDERR "[$prefix][Step 7]\tRunning MED analysis...\n";
	print STDERR "[$prefix][Step 7]\tReformatting sequence description to work with MED...";
	
	# Create a subdirectory MED inside the output prefix for storing the output
	`mkdir -p $outDir/MED`;
	
	# this is the input file for the oligotyping analysis
	my $in = $prefix.".reads.med";
	
	# this is the original merged input file that needs to be edited
	# the descriptor lines are required to be changed
	open FILE, "<$outDir/$prefix.reads.temp" or die "\nERROR (Line ".__LINE__."): Cannot open input file $outDir/$prefix.reads.temp: $!\n";
	open OUT, ">$outDir/MED/$in" or die "\nERROR (Line ".__LINE__."): Cannot open output file $outDir/MED/$in: $!\n";
	my $seqindex = 1;
	my $last = "";
	while(<FILE>){
		if($_ =~ m/^>/){
			$_ =~ s/^>//;
			$_ =~ s/:.*$//;
			$_ =~ s/[^a-zA-Z0-9]//g;
			if($last ne $_){
				$seqindex = 1;
				$last = $_;
			}
			$_ = ">$_"."_$seqindex\n";
			$seqindex++;
		}
		print OUT $_;
	}
	close OUT;
	close FILE;
	
	
	
	if(length($medParameters) != 0){
		my @parts = split(/\s+/, $medParameters);
		$medParameters = "";
		for(my $i=0; $i < @parts; $i+=2){
			if($parts[$i] !~ m/-o/){
				$medParameters .= "$parts[$i] $parts[$i+1] ";
			}
		}
	}
	
	print STDERR "Done!\n\tPadding sequences with gaps...";
	
	`o-pad-with-gaps $outDir/MED/$in -o $outDir/MED/$in.withPadgaps`;
	
	print STDERR "Done!\n";
	print STDERR "\tStarting decompose...";
	
	`decompose -o $outDir/MED $medParameters -E $medMetadataFile $outDir/MED/$in.withPadgaps`; #JCGaby removed --gen-html on July 7 2016
	
	die "\nError: The output from oligotyping was not found! Exiting.\n" if(! -e "$outDir/MED/NODE-REPRESENTATIVES.fasta");
	#my $key = "";
	
	print STDERR "\n\t\tFinish running decompose, initiating BLAST...";
	
	# BLAST the representative sequences against the database to get the name
	my $results = `blastn -query $outDir/MED/NODE-REPRESENTATIVES.fasta -db $db -outfmt "6 qseqid sseqid evalue pident" -evalue 0.01 -max_target_seqs 1 -max_hsps 1 -num_threads $threads 2> $outDir/MED/blast.err | tee $outDir/MED/oligoBlast.tsv`;
	
	open ERR, "<$outDir/MED/blast.err" or die "\nERROR (Line ".__LINE__."): Cannot open input file $prefix/MED/blast.err: $!\n";
	while(<ERR>){
		next if($_ =~ m/Hyphens are invalid and will be ignored/);
		print STDERR $_;
	}
	close ERR;
	
	my @results = split(/\n/, $results);
	my %mapping;
	foreach my $this (@results){
		chomp $this;
		my ($qseqid, $sseqid, $eval, $pident) = split(/\t/, $this);
		$qseqid =~ s/\|.*//;
		$sseqid =~ s/.*;//;
		$sseqid =~ s/_/ /;
		$mapping{$qseqid} = "$sseqid ($pident%)";
	}
	print STDERR "done\n";
	
	print STDERR "\n\t\tFinished BLAST, starting output file processing...";
	
	open FILE, "<$outDir/MED/NETWORK.gexf" or die "Cannot open input file $outDir/MED/NETWORK.gexf: $!\n";
	my @file = <FILE>;
	my $file = join("\n", @file);
	foreach my $key (keys %mapping){
		$file =~ s/$key/$mapping{$key}/g;
	}
	close FILE;
	open OUT, ">$outDir/MED/labelled.NETWORK.gexf" or die "Cannot create output file $outDir/MED/labelled.NETWORK.gexf: $!\n";
	print OUT $file;
	close OUT;
	@file = ();
	
	open FILE, "<$outDir/MED/TOPOLOGY.gexf" or die "Cannot open input file $outDir/MED/TOPOLOGY.gexf: $!\n";
	@file = <FILE>;
	$file = join("\n", @file);
	foreach my $key (keys %mapping){
		$file =~ s/$key/$mapping{$key}/g;
	}
	close FILE;
	open OUT, ">$outDir/MED/labelled.TOPOLOGY.gexf" or die "Cannot create output file $outDir/MED/labelled.NETWORK.gexf: $!\n";
	print OUT $file;
	close OUT;
	print STDERR "done\n";
	print STDERR "[$prefix][Step 7]\tCompleted\n";
}

# Function: singleSampleAnalysis ##########
#	Input: Do analysis of a single sample
#
#	Output: none
#
# 	This function basically gives sequence of
#	instructions (calls function) to do the analysis
###########################################
sub singleSampleAnalysis{
	my $start_run = time();
	my $date = localtime();
	
	my ($out) = @_;
	
	print STDERR "[$out][Step 0]\tStarted task at time ".$date."\n";
	return if(!preprocess($r1, $r2, $out, $prLen, $plLen));
	taxAssign($out, $db);
	
	my $end_run = time();
	$date = localtime();
	print STDERR "[$out][Done]\tTask took: ".($end_run - $start_run)." sec.  Current time: ".$date."\n";
}

# Function: multiSampleAnalysis ####################
#	Input: Do analysis of a single run (multiple samples)
#
#	Output: none
#
# 	This function basically gives sequence of
#	instructions (calls function) to do the analysis
###########################################
sub multiSampleAnalysis{
	
	my @files;
	share(@files);
	my %seen;
	my ($prefix) = @_;
	
	# merge the files
	open FILE, "<$l" or die "\nERROR (Line ".__LINE__."): Cannot open input file $l: $!\n";
	while(<FILE>){
		chomp $_;
		$_ =~ s/^\s+//;
		
		next if($_ =~ m/^\s*#/ || length($_) < 5);
		
		my ($out, $r1, $r2) = split(/\t/, $_);
		next if(defined $seen{$out});
		$seen{$out} = 1;
		
		my $start_run = time();
		
		$out =~ s/.*\///;
		
		my $t = async{push(@files, $out) if(preprocess($r1, $r2, $out, $prLen, $plLen))};
		
		my @joinable = threads->list(threads::joinable);
		$_->join() for(@joinable);
		my @running = threads->list(threads::running);
		while(@running == $threads){
			@running = threads->list(threads::running);
			sleep(1);
		}
	}
	close FILE;
	
	my @running = threads->list(threads::running);
	while(@running > 0){
		my @joinable = threads->list(threads::joinable);
		$_->join() for(@joinable);
		@running = threads->list(threads::running);
		sleep(1);
	}

	mergeReads($prefix, @files);
	taxAssignBatch($prefix, $db);
	oligotype($prefix) if($oligotyping == 1);
	
	print STDERR "[$prefix] All done... bye bye!\n";
}

#############################################################
readTaxonomy();
my $prefix = $outDir;
$prefix =~ s/.*\///;
if(defined $l) {
	multiSampleAnalysis($prefix);
} else {
	singleSampleAnalysis($prefix);
}
