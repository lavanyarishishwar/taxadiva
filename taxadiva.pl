#!/usr/bin/perl -w
#
# Authors        : Lavanya Rishishwar, Chris Gaby
# Creation Date  : 23rd Aug 2015
# Last Modified  : 15th July 2016
# Version        : 0.10
#
#############################################################
use strict;
use Getopt::Long;
use threads;
use threads::shared;
#############################################################

#   Pretty usage
my $programHead = "TaxADivA - TAXonomic Assignment and DIVersity Assessment for amplicon reads\n";
my $baseUsage = "           [-o <STRING. output prefix to store results. Default: input f>]
           [-d <STRING. Database. Default: db.fasta>] [-t <STRING. Taxonomy file. Default: tax.tsv>]
           [-p <INT. Primer length to be trimmed. Default: no trimming>]
           [-r <INT. Primer length to be trimmed from the RIGHT. Default: no trimming>]
           [-l <INT. Primer length to be trimmed from the LEFT. Default: no trimming>]
           [-k <FLAG. If specified, the trimmed primers will be retained as separate files.  Default: Don't retain trimmed primers.>]
           [-y <FLAG. Performs MED analysis>]
           [-j <INT. Number of threads. Default: 1>]
           [-g <INT. Depth cutoff for considering file. Default: 5000>]
           [--pear <STRING. PEAR paramaters in double quotes.  Default: \"-v 50 -m 450 -n 350 -p 1.0 -j <threads>\">]
           [--med <STRING. Oligotyping paramaters in double quotes.  Default: \"\">]
           [--med-metadata <STRING. MED metadata file to be used for the decompose command.  Specified by the -E option in the decompose command.  Default: decompose_map3.tab.>]
           [-u <STRING. USEARCH program path. Default: usearch>]
           [-h <FLAG. Prints this help>]";
my $usage = "\n=====================================\n$programHead\n=====================================\n$0  [-1 <STRING. forward read file>] [-2 <STRING. reverse read file>]\n$baseUsage\n\n$0  [-s <STRING. file with set of forward and reverse files>]\n$baseUsage \n=====================================\nExample usage: $0 -d db.fasta -t tax.db.tsv -j 18 -s list.txt -o output1 -m \"-v 50 -m 450 -n 350 -p 1.0 -j 4\"\n";

#############################################################
#	Global Variables
my $r1;
my $r2;
my $l;
my $help = 0;
my $b = "beta.tsv";
my $outFile = "out"; # Output file
my $db      = "db.fasta"; # Database
my $taxFile = "tax.tsv"; # Database
my $threads = 10;
my $pLen;
my $prLen;
my $plLen;
my $keepTrimmedPrimers = 0;
my $usearchProgram = "usearch";
my $pearParameters = "-v 50 -m 450 -n 350 -p 1.0 -j $threads";
my $medParameters = "";
my $goodDepth = 5000;
my $oligotyping = 0;
my $medMetadataFile = "decompose_map3.tab";
my %tax; my %tax2Species; my %family; my %genus; my %class; my %uniqueClass; my %taxid; my %spfFormat;
share(%tax); share(%tax2Species); share(%family); share(%genus); share(%class); share(%uniqueClass); share(%taxid);
#############################################################

#	Argument Check
my $args  = GetOptions ("1=s"     => \$r1,
                        "2=s"     => \$r2,
                        "s=s"     => \$l,
                        "b=s"     => \$b,
                        "h"       => \$help,
                        "d=s"     => \$db,
                        "u=s"     => \$usearchProgram,
                        "t=s"     => \$taxFile,
                        "j=i"     => \$threads,
                        "p=i"     => \$pLen,
                        "r=i"     => \$prLen,
                        "l=i"     => \$plLen,
                        "k+"      => \$keepTrimmedPrimers,
                        "g=i"     => \$goodDepth,
                        "y+"      => \$oligotyping,
                        "o=s"     => \$outFile,
                        "pear=s"  => \$pearParameters,
                        "med=s"   => \$medParameters,
                        "med-metadata=s"   => \$medMetadataFile);
                        
#############################################################

print STDERR "Checking for the provided arguments...";
	
if($help == 1){
	print STDERR $programHead;
	print STDERR $usage;
	exit;
}

                        
if(! defined $l && !( defined $r1 && defined $r2)){
	print STDERR "\nERROR (Line ".__LINE__."): list file (option -s) or forward (option -1) and reverse (option -2) reads are mandatory options.  Please specify them!\n";
	print STDERR "$usage\n";
	exit;
}
if(defined $l){
	if(! -e $l){
		print STDERR "\nERROR (Line ".__LINE__."): Please make sure that list file (option -s) exist or the location provided is correct!\n";
		print STDERR "$usage\n";
		exit;
	}
} else {
	if(! -e $r1 || ! -e $r2){
		print STDERR "\nERROR (Line ".__LINE__."): Please make sure that the forward (option -1) and reverse (option -2) reads exist or the location provided is correct!\n";
		print STDERR "$usage\n";
		exit;
	}
}
if(! -e $taxFile){
	print STDERR "\nERROR (Line ".__LINE__."): Please make sure that the taxonomy (option -t) file exist or the location provided is correct!\n";
	print STDERR "$usage\n";
	exit;
}
if(! -e $db){
	print STDERR "\nERROR (Line ".__LINE__."): Please make sure that the database (option -d) file exist or the location provided is correct!\n";
	print STDERR "$usage\n";
	exit;
}
if( ! -e "$db.nhr" || ! -e "$db.nin" || ! -e "$db.nsq" ){
	`makeblastdb -in $db -out $db -dbtype nucl`;
}

if(length($outFile) == 0 || $outFile =~ m{[\\:*?"<>|]}){
	print STDERR "\nERROR (Line ".__LINE__."): Invalid output prefix.  Please provide a prefix with acceptable characters (alphanumeric, _, .)\n";
	print STDERR "$usage\n";
	exit;
}

if($outFile =~ m/\//){
	my @parts = split(/\//, $outFile);
	pop(@parts);
	my $prefix = join("/", @parts);
	`mkdir -p $prefix` if(! -e $prefix);
}

print STDERR "everything looks fine.\n";

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
my $blastVersion = `blastn -version | head -1 | sed 's/blastn: //'`;
chomp $blastVersion;
if($blastVersion =~ m/^2./){
	# The BLAST version is 2.2.XX
	$blastVersion =~ s/2.2.//;
	$blastVersion =~ s/\+$//;
	if($blastVersion < 31){
		die "An older version of BLAST is detected (2.2.$blastVersion).  Please upgrade to BLAST+ version 2.2.31+.  The BLAST latest binaries can be obtained from ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/\n";
	}
}
$programPath = "";

# USEARCH
$programPath = `which $usearchProgram`;
die "\nERROR (Line ".__LINE__."): USEARCH not found! Please make sure USEARCH is installed and is in the \$PATH variable.\nTo install USEARCH, download the binaries from http://www.drive5.com/usearch/download.html and place it in your \$PATH\n" if(length($programPath) == 0);
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
		
		$spfFormat{$accn} = "k__".$kingdom."\tp__".$phylum."\tc__".$class."\to__".$order."\tf__".$family."\tg__".$genus."\ts__".$species;
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
	
	`pear -f $r1 -r $r2 -o $out-pear $pearParameters`; # Need to revisit this command in future
	
	print STDERR "[$out][Step 1]\tConverting FASTQ to FASTA...\n";
	`awk 'BEGIN{c=0} {c++; if(c==1){print ">"\$0}; if(c==2){print}; if(c==4){c=0}}' $out-pear.assembled.fastq > $out-pear.fa`;
	
	if($pRight+$pLeft > 0){
		print STDERR "[$out][Step 1]\tTrimming primers...\n";
		`prinseq-lite.pl -fasta $out-pear.fa -out_good $out -trim_left $pLeft -trim_right $pRight 1>&2 2>> $out.log`;
		if($keepTrimmedPrimers > 0){
			open FILE, "<$out-pear.fa" or die "Cannot open PEAR output $out-pear.fa: $!\n";
			open LEFT, ">$out.trimmedLeft.fa" or die "Cannot create output file $out.trimmedLeft.fa: $!\n";
			open RIGHT, ">$out.trimmedRight.fa" or die "Cannot create output file $out.trimmedRight.fa: $!\n";
			
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
		`mv $out.fasta $out.reads.temp`;
		print STDERR "[$out][Step 1]\tTrimming Done\n";
	} else {
		`mv $out-pear.fa $out.reads.temp`;
	}
	
	`rm $out-pear*`;
	
	return 1;
}

# Function: mergeReads ####################
#	Input: merged files that passed 
#			preprocessing
#
#	Output: Single merged file name
#
#	This will generate $out$$.read.temp
###########################################
sub mergeReads{
	my ($prefix, @files) = @_;
	my $out = $prefix.".reads.temp";
	print STDERR "[$prefix][Step 1]\tMerging processed reads (total of ".@files." files)...";
	open OUT, ">$out" or die "\nERROR (Line ".__LINE__."): Cannot create output file $out:$!\n";
	
	foreach my $file (@files){
		
		open FILE, "<$file.reads.temp" or die "\nERROR (Line ".__LINE__."): Cannot open input file $file.reads.temp:$!\n";
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
	`$usearchProgram -derep_fulllength $out.reads.temp -fastaout $out.reads.derep -sizeout -threads $threads 1>> $out.log 2>> $out.log`;
	#`rm $out.reads.temp`;
	
	print STDERR "Done\n[$out][Step 3]\tClustering sequences...";
	`$usearchProgram -cluster_otus $out.reads.derep -otus $out.otus.fa -uparseout $out.up -relabel OTU -minsize 2 1>> $out.log 2>> $out.log`;
	
	print STDERR "Done\n[$out][Step 3]\tRemoving chimeras...";
	`$usearchProgram -uchime_ref $out.otus.fa -db $db -strand plus -minh 1.0 -nonchimeras $out.nochim.fa -uchimeout $out.uchime -uchimealns $out.aln 1>> $out.log 2>> $out.log`;
	
	# Old paradigm
	#`usearch -uchime_ref $out.reads.temp -db $db -strand plus -minh 0.28 -nonchimeras $out.reads.fasta 1>> $out.log 2>> $out.log`;
	#`usearch -cluster_fast $out.reads.fasta -id 0.9 -centroids $out.cent.fa -uc $out.clusters.tsv -strand both -sizeout >> $out.log 2>> $out.log`;
	
	# Get the count of total number of sequences processed from the centroid file
	
	my $total = 0;
	my %otuSizes;
	my %otu2Keep;
	my $otu2Keep = `grep '^>' $out.nochim.fa`;
	my @otu2Keep = split(/\n/, $otu2Keep);
	foreach(@otu2Keep){
		chomp $_;
		$_ =~ s/^>//;
		$otu2Keep{$_} = 1;
	}
	
	open FILE, "<$out.up" or die "\nERROR (Line ".__LINE__."): Cannot open input file $out.clusters.tsv: $!\n";
	while(<FILE>){
		chomp $_;
		my ($query, $type, @vars) = split(/\t/, $_);
		my $otu = pop(@vars);
		
		next if(! defined $otu2Keep{$otu}); #Throw out chimeric OTUs
		
		my $size = $query;
		$size =~ s/.*;size=//;
		$size =~ s/;\s*$//;
		
		$total += $size;
		$otuSizes{$otu} += $size;
	}
	close FILE;
	
	# Perform the BLAST search
	print STDERR "Done\n[$out][Step 4]\tLooking for homologs...";
	my $blOut = `blastn -query $out.nochim.fa -db $db -outfmt "6 qseqid sseqid evalue pident" -evalue 0.01 -max_target_seqs 1  -max_hsps 1 -num_threads $threads | tee $out.blast.tsv`;
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
			if($s =~ /cluster_IV/){
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
	open OUT, ">$out.abundance.txt" or die "\nERROR (Line ".__LINE__."): Cannot create output file $out.abundance.txt: $!\n";
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
	open OUT, ">$out.krona.txt" or die "\nERROR (Line ".__LINE__."): Cannot create output file $out.krona.txt: $!\n";
	foreach my $key (keys %counts){
		my $text = $key;
		$text =~ s/\//&#47;/g;
		print OUT "$counts{$key}\t$text\n";
	}
	print OUT ($total-$assigned)."\tUnassigned OTU\n";
	close OUT;
	print STDERR "Done\n\tMaking plots...";
	`ktImportText -n nifH $out.krona.txt -o $out.krona.html`;
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
	`$usearchProgram -derep_fulllength $out.reads.temp -fastaout $out.reads.derep -sizeout -threads $threads 1>> $out.log 2>> $out.log`;
	#`rm $out.reads.temp`;
	
	print STDERR "Done\n[$out][Step 3]\tClustering sequences...";
	`$usearchProgram -cluster_otus $out.reads.derep -otus $out.otus.fa -uparseout $out.up -relabel OTU -minsize 2 1>> $out.log 2>> $out.log`;
	
	print STDERR "Done\n[$out][Step 3]\tRemoving chimeras...";
	`$usearchProgram -uchime_ref $out.otus.fa -db $db -strand plus -minh 1.0 -nonchimeras $out.nochim.fa -uchimeout $out.uchime -uchimealns $out.aln 1>> $out.log 2>> $out.log`;
	
	# Old paradigm
	#`usearch -uchime_ref $out.reads.temp -db $db -strand plus -minh 0.28 -nonchimeras $out.reads.fasta 1>> $out.log 2>> $out.log`;
	#`usearch -cluster_fast $out.reads.fasta -id 0.9 -centroids $out.cent.fa -uc $out.clusters.tsv -strand both -sizeout >> $out.log 2>> $out.log`;
	
	# Perform the BLAST search
	print STDERR "Done\n[$out][Step 4]\tLooking for homologs...";
	my $blOut = `blastn -query $out.nochim.fa -db $db -outfmt "6 qseqid sseqid evalue pident" -evalue 0.01 -max_target_seqs 1 -max_hsps 1 -num_threads $threads | tee $out.blast.tsv`;
	print STDERR "Done\n[$out][Step 5]\tProcessing results...";
	
	# Old Paradigm
	#print STDERR "[$out][Step 2]\tRemoving chimeras...";
	#`usearch -uchime_ref $out.reads.temp -db $db -strand plus -minh 0.28 -nonchimeras $out.reads.fasta 1> $out.log 2>> $out.log`;
	#`rm $out.reads.temp`;
	#print STDERR "Done\n[$out][Step 3]\tClustering sequences...";
	#`usearch -cluster_fast $out.reads.fasta -id 0.9 -centroids $out.cent.fa -uc $out.clusters.tsv -strand both -sizeout > $out.log 2>> $out.log`;
	
	
	# These are emperically calculated threshold values
	my $family  = 75;
	my $genus   = 88.1;
	my $species = 91.9;
	
	# These variables will store the processed results
	my %assignments; # stores mapping of centroids to taxonomy
	
	my @blast = split(/\n/, $blOut);
	
	for(my $i=0; $i < @blast; $i++){
		my ($q, $s, $e, $p) = split(/\t/, $blast[$i]);
		if($p > $family){
			
			$q =~ s/;size.*$//;
			
			if($s =~ /cluster_IV/){
				$assignments{$q}{"class"} = "cluster_IV";
				next;
			}
			
			my ($accn, undef) = split(/;/, $s);
			
			$assignments{$q}{"class"} = $class{$accn};
			$assignments{$q}{"accn"} = $accn;
			
			if($p > $species){
				die "\nERROR (Line ".__LINE__."): Cannot find $accn from $q! Do all the sequences in the database exist in the taxonomy file (or is the database build an different version from taxonomy file)? Exiting\n" if(! defined $tax{$accn});
				$assignments{$q}{"assign"} = $tax{$accn};
				#$assignments{$q}{"alpha"} = $tax2Species{$accn};
			} elsif($p > $genus) {
				die "\nERROR (Line ".__LINE__."): Cannot find $accn from $q! Do all the sequences in the database exist in the taxonomy file (or is the database build an different version from taxonomy file)? Exiting\n" if(! defined $genus{$accn});
				$assignments{$q}{"assign"} = $genus{$accn};
			} else {
				die "\nERROR (Line ".__LINE__."): Cannot find $accn from $q! Do all the sequences in the database exist in the taxonomy file (or is the database build an different version from taxonomy file)? Exiting\n" if(! defined $family{$accn});
				$assignments{$q}{"assign"} = $family{$accn};
			}
		}
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
	my %biom;
	my %biomBins;
	
	my %otu2Keep;
	my $otu2Keep = `grep '^>' $out.nochim.fa`;
	my @otu2Keep = split(/\n/, $otu2Keep);
	foreach(@otu2Keep){
		chomp $_;
		$_ =~ s/^>//;
		$otu2Keep{$_} = 1;
	}
	
	open FILE, "<$out.up" or die "\nERROR (Line ".__LINE__."): Cannot open input file $out.up: $!\n";
	while(<FILE>){
		chomp $_;
		my ($query, $type, @vars) = split(/\t/, $_);
		my $otu = pop(@vars);
		
		next if(! defined $otu2Keep{$otu}); #Throw out chimeric OTUs
		
		my ($bin, undef) = split(/:/, $query);
		my $size = $query;
		$size =~ s/.*;size=//;
		$size =~ s/;\s*$//;
		
		if(defined $assignments{$otu}{"class"}){
			next if($assignments{$otu}{"class"} eq "cluster_IV"); # Precautionary check, should never be true
			
			if(defined $assignments{$otu}{"class"}){
				$classCounts{$bin}{$assignments{$otu}{"class"}} += $size;
				$classesSeen{$assignments{$otu}{"class"}} = 1;
			}
			
			$assigned{$bin} += $size;
		} else {
			$assignments{$otu}{"assign"} = "OTU$otuNum";
			$otuNum++;
		}
		
		$biom{$otu}{$bin} += $size;
		$biomBins{$bin} = 1;
		
		$clusterCount{$bin}{$assignments{$otu}{"assign"}} += $size;
		$total{$bin} += $size;
	}
	close FILE;
	
	print STDERR "Done\n[$out][Step 6]\tOutputting results to files...";
	
	# Class abundance files
	open ABND, ">$out.abundance.txt" or die "\nERROR (Line ".__LINE__."): Cannot create output file $out.abundance.txt: $!\n";
	foreach my $key (keys %classesSeen){
		next unless(acceptable($key));
		print ABND "\t$key";
	}
	print ABND "\tUnassignedOTUs\n";
	foreach my $bin (keys %classCounts){
		print ABND "$bin";
		my $percent = 0;
		foreach my $tax (keys %classesSeen){
			next unless(acceptable($tax));
			$classCounts{$bin}{$tax} //= 0;
			my $this = $classCounts{$bin}{$tax}*100/$total{$bin};
			print ABND "\t$this";
			$percent += $this;
		}
		print ABND "\t".(100-$percent)."\n";
	}
	close ABND;
	print STDERR "\n\tCreated $out.abundance.txt";
		
	# Set of krona files
	foreach my $bin (keys %clusterCount){
		open KRONA, ">$bin.krona.txt" or die "\nERROR (Line ".__LINE__."): Cannot create output file $bin.krona.txt: $!\n";
		foreach my $key (keys %{$clusterCount{$bin}}){
			my $text = $key;
			$text =~ s/\//&#47;/g;
			next if($text =~ m/^OTU/);
			print KRONA "$clusterCount{$bin}{$key}\t$text\n";
		}
		print KRONA ($total{$bin}-$assigned{$bin})."\tUnassigned OTU\n";
		close KRONA;
		`ktImportText -n nifH $bin.krona.txt -o $bin.krona.html`;
	}
	print STDERR "\n\tCreated Krona files (Total files = ".(scalar keys %clusterCount).")";
	
	
	# This is where the BIOM file is generated for input for QIIME
	
	my @biomBins = keys %biomBins;
	my @otus = keys %biom;
	
	open BIOM, ">$out.otu_table.txt" or die "\nERROR (Line ".__LINE__."): Cannot create output file $out.otu_table.txt: $!\n";
	print BIOM "# Constructed from biom file\n";
	print BIOM "#OTU ID";
	foreach (@biomBins){
		print BIOM "\t$_";
	}
	print BIOM "\ttaxonomy\n";
	for(my $i = 0; $i < @otus; $i++){
		print BIOM "denovo".($i+1);
		foreach my $bin (@biomBins){
			$biom{$otus[$i]}{$bin} //= 0;
			print BIOM "\t$biom{$otus[$i]}{$bin}";
		}
		#The following code would've printed the taxonomy ID
		#print BIOM "\t".$taxid{$assignments{$otus[$i]}{"accn"}} if (defined $assignments{$otus[$i]}{"accn"});
		#The new one will print the taxonomy hierarchy
		next if(! defined $assignments{$otus[$i]}{"accn"});
		if(defined $spfFormat{$assignments{$otus[$i]}{"accn"}}){
			my $name = $spfFormat{$assignments{$otus[$i]}{"accn"}};
			$name =~ s/\t/;/g;
			print BIOM "\t$name";
		}
		print BIOM "\n";
	}
	close BIOM;
	print STDERR "\n\tCreated $out.otu_table.txt";
	
	
	# This is where the STAMPS spf file is generated
	
	open SPF, ">$out.spf" or die "\nERROR (Line ".__LINE__."): Cannot create output file $out.spf: $!\n";
	print SPF "Domain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies";
	foreach (@biomBins){
		print SPF "\t$_";
	}
	print SPF "\n";
	for(my $i = 0; $i < @otus; $i++){
		next if(!defined $assignments{$otus[$i]}{"accn"});
		print SPF $spfFormat{$assignments{$otus[$i]}{"accn"}};
		foreach my $bin (@biomBins){
			$biom{$otus[$i]}{$bin} //= 0;
			print SPF "\t$biom{$otus[$i]}{$bin}";
		}
		print SPF "\n";
	}
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
	
	print STDERR "[$prefix][Step 7]\tRunning MED analysis...";
	print STDERR "[$prefix][Step 7]\t\tReformatting sequence description to work with MED...";
	
	# this is the input file for the oligotyping analysis
	my $in = $prefix.".reads.med";
	
	# this is the original merged input file that needs to be edited
	# the descriptor lines are required to be changed
	open FILE, "<$prefix.reads.temp" or die "\nERROR (Line ".__LINE__."): Cannot open input file $prefix.reads.temp: $!\n";
	open OUT, ">$in" or die "\nERROR (Line ".__LINE__."): Cannot open output file $in: $!\n";
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
			$_ = "Sample-".$_."_Read$seqindex\n";
			$seqindex++;
		}
		print OUT $_;
	}
	close OUT;
	close FILE;
	
	
	# Create a subdirectory MED inside the output prefix for storing the output
	`mkdir -p $prefix/MED`;
	
	if(length($medParameters) != 0){
		my @parts = split(/\s+/, $medParameters);
		$medParameters = "";
		for(my $i=0; $i < @parts; $i+=2){
			if($parts[$i] !~ m/-o/){
				$medParameters .= "$parts[$i] $parts[$i+1] ";
			}
		}
	}
	
	`o-pad-with-gaps $in -o $in.withPadgaps`;
	`decompose -o $prefix/MED $medParameters $in.withPadgaps -E $medMetadataFile`; #JCGaby removed --gen-html on July 7 2016
	close OUT;
	
	#my $key = "";
	
	print STDERR "\n\t\tFinish running decompose, initiating BLAST...";
	
	# BLAST the representative sequences against the database to get the name
	my $results = `blastn -query $prefix/MED/NODE-REPRESENTATIVES.fasta -db $db -outfmt "6 qseqid sseqid evalue pident" -evalue 0.01 -max_target_seqs 1 -max_hsps 1 -num_threads $threads | tee $prefix/MED/oligoBlast.tsv`;
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
	
	open FILE, "<$prefix/MED/NETWORK.gexf" or die "Cannot open input file $prefix/MED/NETWORK.gexf: $!\n";
	my @file = <FILE>;
	my $file = join("\n", @file);
	foreach my $key (keys %mapping){
		$file =~ s/$key/$mapping{$key}/g;
	}
	close FILE;
	open OUT, ">$prefix/MED/labelled.NETWORK.gexf" or die "Cannot create output file $prefix/MED/labelled.NETWORK.gexf: $!\n";
	print OUT $file;
	close OUT;
	@file = ();
	
	open FILE, "<$prefix/MED/TOPOLOGY.gexf" or die "Cannot open input file $prefix/MED/TOPOLOGY.gexf: $!\n";
	@file = <FILE>;
	$file = join("\n", @file);
	foreach my $key (keys %mapping){
		$file =~ s/$key/$mapping{$key}/g;
	}
	close FILE;
	open OUT, ">$prefix/MED/labelled.TOPOLOGY.gexf" or die "Cannot create output file $prefix/MED/labelled.NETWORK.gexf: $!\n";
	print OUT $file;
	close OUT;
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
	my ($r1, $r2, $db, $out, $primerTrimLen) = @_;
	
	print STDERR "[$out][Step 0]\tStarted task at time ".$date."\n";
	return if(!preprocess($r1, $r2, $outFile, $pLen));
	taxAssign($outFile, $db);
	
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
	
	# merge the files
	open FILE, "<$l" or die "\nERROR (Line ".__LINE__."): Cannot open input file $l: $!\n";
	while(<FILE>){
		next if($_ =~ m/^#/);
		chomp $_;
		my ($out, $r1, $r2) = split(/\t/, $_);
		next if(defined $seen{$out});
		$seen{$out} = 1;
		
		my $start_run = time();
		
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
	
	mergeReads($outFile, @files);
	taxAssignBatch($outFile, $db);
	oligotype($outFile) if($oligotyping == 1);
	
	print STDERR "[$outFile] All done... bye bye!\n";
}

#############################################################
readTaxonomy();
if(defined $l) {
	multiSampleAnalysis();
} else {
	singleSampleAnalysis();
}
