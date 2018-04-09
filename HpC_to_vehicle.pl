#!/usr/bin/perl

=head1 NAME

	HpC_to_vehicle.pl
	
=head1 SYNOPSIS

	Vehicle analysis with genes :
	perl HpC_to_vehicle.pl -equi equivalence.txt -netin Network.csv -tabin Table.csv -seq nt -ppA 'Analyse/muscle/*' -pg 'directoryGene/*'
	
	Vehicle analysis with proteins :
	perl HpC_to_vehicle.pl -equi equivalence.txt -netin Network.csv -tabin Table.csv -seq aa -ppM 'Analyse/phylip/'

=head1 DESCRIPTION

	The perl script HpC_to_vehicle retrieves the results of the perl script AccNetPhylA (files to visualize the networks, multiple alignment of Muscle, phylogenetic results of PHYLIP, ...) 
	to create HpC (Homologous protein Cluster) and GU (Genomic Unit), to thus make different networks representing the links between "vehicles" of genes or proteins.

	For each vehicle, there are two types of graphs:
		-network between vehicles that can belong to the same GU
		-network between different GU vehicles

	For gene vehicle analysis, a patristic distance (calculated from tree branch lengths describe the amount of genetic change represented by a tree) will be used to build two more graphs.

	In the end, there will be 4 networks resulting from an analysis with the genes and 2 from an analysis with the proteins.

=head1 OPTION

	-Help|help|h, produces this help file. 
	
	-Verbose|verbose|v, Verbose mode. Can be used multiple times for increased verbosity.
	
	-equi, Equivalence file of AccNetPhylA.pl.
	
	-netin, Network filename of AccNetPhylA.pl.
	
	-tabin, Table filename of AccNetPhylA.pl.

	-seq, Proteic or nucleic sequence (nt|aa).
	
	-ppA, Path of the protein alignment files to analyze (with -seq nt).
	
	-ppM, Path of the protein matrix files to analyze (with -seq aa).
	
	-pg, Path of the fasta files of genes to to analyze (with -seq nt).
	
	-clean, Remove and class files. Default(yes).

	-dir, Directory's name. Default(Analyse_vehicle).
	

=head1 AUTHORS

	SIGURET Clea
	
	
=head1 VERSION

	V1.0


=head1 DATE

	Creation : 20.06.2017
	Last modification : 23.03.2018

=cut 

# librairies 

use strict;
use warnings;

use Getopt::Long;
use Pod::Usage;
use FindBin;
use Statistics::R;
use Data::Dumper;
use List::Util qw(min);
use Graph::Undirected;#graph data structures and algorithms (version 0.9704)

# scalars 
my $help;
my $verbose;
my $version;
my $VERSION = 1.0 ;
my $force;
my $exe = $FindBin::RealScript;
my $R = Statistics::R->new();
my $g0 = Graph::Undirected->new;
my $g1 = Graph::Undirected->new;
my $g2 = Graph::Undirected->new;
my $g3 = Graph::Undirected->new;
my $g4 = Graph::Undirected->new;
my $g5 = Graph::Undirected->new;

my $equi;
my $netin;
my $tabin;
my $seqTY;
my $pathProtA;
my $pathProtM;
my $pathGene;
my $clean = "yes";
my $dir = "Analyse_vehicle";
my $NbClustBon = 0;
my $idClust;
my $idPC;
my $idGFnum;
my $idGFcont;
my $tree;
my $actualNumCluster = 0;

# tables
my @linet;
my @litab;
my @AliProt;
my @FastaGene;
my @clMatr;
my @clMatrG;

# hashes 
my %hashEqui;
my %hashEquiG;
my %hashnet;
my %comptarget;
my %hashtab;
my %GU;
my %twin;
my %twtc;
my %clBon;
my %hashProtClust;
my %hashProtList;
my %hashGeneFasta;
my %hashG;
my %hashMCG;
my %hashMCG2;


# function

sub error {
	#~ management of error messages and help page layout, will stop execution
	#~ local arguments passed:1st, error message to output
	my $error = shift;
	my $filename = ($0);
	pod2usage(-message => "$filename (error): $error. Execution halted.", -verbose => 2, -noperldoc => 1);
	exit(2);
}

sub warning {
	#~ management of warnings and execution carry on
	#~ local arguments passed: 1st, warning message to output
	my $message = "@_";
	if ($verbose) {
		my $filename = $0;
		#~ warn("$filename (info): ".$message."\n");
		warn($message."\n");
	}
}

sub version {
  print STDERR "$exe $VERSION\n";
  exit;
}


#~ ×××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××
#~ Get options

MAIN: {
	GetOptions(
			"Help|help|h" 			=> \$help,				#help flag
			"Verbose|verbose|v!" 	=> \$verbose,			#verbose flag
			"version!"				=> \$version,			#version flag
			"force!"				=> \$force,				#force flag
			"equi=s"				=> \$equi,				#Equivalence file of AccNetPhylA.pl
			"netin=s"				=> \$netin,				#Network file of AccNetPhylA.pl
			"tabin=s"				=> \$tabin,				#Table file of AccNetPhylA.pl
			"seq=s"					=> \$seqTY,				#Proteic or nucleic sequence
			"ppA=s"					=> \$pathProtA,			#Path of the protein alignment files to analyze
			"ppM=s"					=> \$pathProtM,			#Path of the protein matrix files to analyze
			"pg=s"					=> \$pathGene,			#Path of the fasta files of genes to to analyze
			"clean=s"				=> \$clean,				#Remove and class files
			"dir=s"					=> \$dir,				#Directory's name
		);
	
	version() if ($version);
	
	if ($help) {
		pod2usage(-verbose => 2, -noperldoc => 1);
		exit;
	}
		
	unless($equi){
		unless($help){print STDERR "You must specified the equivalence file\n";exit;}
	}
	
	unless($netin){
		unless($help){print STDERR "You must specified the network file\n";exit;}
	}
	
	unless($tabin){
		unless($help){print STDERR "You must specified the table file\n";exit;}
	}
	
	unless($seqTY){
		unless($help){print STDERR "You must specified the sequence type : protein (aa) or gene (nt)\n";exit;}
	}
	
	if($seqTY eq "nt"){
		unless($pathProtA){
			unless($help){print STDERR "You must specified the path of the protein alignment files\n";exit;}
		}
				
		unless($pathGene){
			unless($help){print STDERR "You must specified the path of the fasta files of genes\n";exit;}
		}	
	}
	
	if($seqTY eq "aa"){
		unless($pathProtM){
			unless($help){print STDERR "You must specified the path of the protein matrix files\n";exit;}
		}
	}
	
#~ ×××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××××
#~ Main

print "Opening input files :\n";


######### Equivalence file
print "\t- Equivalence file : ";

open(F1, $equi) or die ("Opening Problem $equi");
my %hashPROT;
while(my $li = <F1>){
	chomp $li;
	#Recovery of the identifier of the sequence and belonging to the contig / species
	if($li =~ m/^>([0-9]+\|[0-9]+)\t>(.+)$/){
		$hashPROT{$1} = $2;
		my @liEqui = split("\t",$li);
		if($liEqui[1] =~ m/^>gi\|\S+\|ref\|/){
			if($liEqui[0] =~ m/^>([0-9]+)(\|[0-9]+)$/){
				my $fich = $1;
				my $id1 = $1.$2;
				#Recovery of the species (for bacteria / archaea)
				if($liEqui[1] =~ m/^>gi\|\S+\|ref\|(.+)\|.+\[([A-Za-z]+.+)\]$/){
					my $id2 = $1;
					my $sp = $2;
					my $val = $fich."_".$sp;
					if ($seqTY eq "nt"){
						$hashEquiG{$id1}{$sp}{$id2} = $fich;
					}
					$hashEqui{$id1} = $val;
				}
				
			}
		}
		else{
			if($liEqui[0] =~ m/^>([0-9]+)(\|[0-9]+)$/){
				my $fich = $1;
				my $id1 = $1.$2;
				if($liEqui[1] =~ m/^>ref\|.+\|(.+)\|(.+)$/){
					my $id2 = $2;
					my $sp = $1;
					my $val = $fich."_".$sp;
					if ($seqTY eq "nt"){
						$hashEquiG{$id1}{$sp}{$id2} = $fich;
					}
					$hashEqui{$id1} = $val;
				}
				#Recovery of contig (for prophages or viruses)
				elsif(($liEqui[1] =~ m/^>(\S+).protein_([0-9]+)/) || ($liEqui[1] =~ m/^>(\S+)\.([0-9]+)/)){
					my $val = $fich."_".$1;
					if ($seqTY eq "nt"){
						$hashEquiG{$id1}{$1}{$2} = $fich;
					}
					
					$hashEqui{$id1} = $val;
					
				}
			}
		}
	}
}
close F1;
print "OK\n";

foreach my $i (keys %hashEqui){
	if(($hashEqui{$i} =~ m/^([0-9]+)_gi\|\S+\|ref\|.+\|.+\[(.+)$/) || ($hashEqui{$i} =~ m/^([0-9]+)_gi\|\S+\|ref\|.+\|.+\[(.+)$/)){
		$hashEqui{$i} = "$1_$2";
	}
}
	
######### Table file
print "\t- Table file : ";

open(FT, $tabin) or die ("Opening Problem $tabin");

while(my $li = <FT>){
	chomp($li);
	if($li =~ m/^Cluster/){
		@litab = split ("\t",$li);
		if($litab[2] =~ m/^Twin/){
			$twin{$litab[0]}{$litab[2]}=$litab[3];
			$twtc{$litab[2]}{$litab[0]}=$litab[3];
		}
		else{
			$hashtab{$litab[0]}{$litab[2]}=$litab[3];
		}
	}
	else{
		@litab = split ("\t",$li);
		if($litab[1] eq "GU"){
			$GU{$litab[0]}=$litab[3];
		}
	}
}
close FT;
print "OK\n";


######### Network file
print "\t- Network file : ";

open (FN, $netin) or die ("Opening Problem $netin");

while(my $li = <FN>){
	chomp($li);
	@linet = split ("\t",$li);
	if($linet[0] =~ m/^Cluster/){
		foreach my $tw (keys %{$twin{$linet[0]}}){
			if($twin{$linet[0]}{$tw}){
				$hashnet{$tw}{$linet[1]}=$linet[2];
				$comptarget{$tw}=2;
			}
		}
		foreach my $tg (keys %{$hashtab{$linet[0]}}){
			$hashnet{$linet[0]}{$linet[1]}=$linet[2];
			++$comptarget{$linet[0]};
		}
	}
}
close FN;
print "OK\n\n";


######### Recovery of clusters of interest (connected to at least two GUs)
print "Recovery of clusters of interest (connected to at least two GUs) : ";

foreach my $tw (keys %twtc){
	if($tw=~ m/^Twin/){
		foreach my $tc (keys %{$twtc{$tw}}){
			if($tc =~ m/^(Cluster_[0-9]+_[0-9]+)/){
				my $key = $1;
				$clBon{$key} = 1;
				++$NbClustBon;
			}
		}
	}
}

foreach my $s (%comptarget){
	if($s =~ m/^Cluster/){
		if($comptarget{$s} ne "1"){
			foreach my $ty(keys %{$hashtab{$s}}){
				if($ty=~ m/^Single/){
					if($s =~ m/^(Cluster_[0-9]+_[0-9]+)/){
						my $key = $1;
						$clBon{$key} = 1;
						++$NbClustBon;
					}
				}
			}
		}	
	}
}
######### End recovery
print "done\n\n";

my %hashGeneFista;
#########Vehicle analysis with genes

if($seqTY eq "nt"){
	
	print "Vehicle analysis with genes : \n";
	
	######### Protein alignment files
	@AliProt = `ls $pathProtA`;

	foreach my $ap (@AliProt){
		if($ap =~ /(Cluster_[0-9]+_[0-9]+).aln/){
			$idClust = $1;
			if($clBon{$idClust}){
				
				open(AP, $ap) or die ("Opening Problem $ap");
				while(my $li = <AP>){
					chomp($li);
					if($li =~ /^>([0-9]+\|[0-9]+)$/){
						$idPC = $1;
						push(@{$hashProtList{$idClust}} , $idPC);
					}
					else{
						$hashProtClust{$idClust}{$idPC} .= $li;
					}
				}
				close AP;	
			}
		}
	}

	open(FILEG1, ">bilan.gene.vehicule.txt") or die ("Error\n");
	print FILEG1 "Vehicule1\tVehicule2\tGeneVehicule1\tGeneVehicule2\n";

	open(FILEG2, ">bilan2.gene.vehicule.txt") or die ("Error\n");
	print FILEG2 "Vehicule1\tVehicule2\tGeneVehicule1\tGeneVehicule2\n";

	open(FILEG1M, ">bilan.gene.vehiculeMC.txt") or die ("Error\n");
	print FILEG1M "Vehicule1\tVehicule2\tGeneVehicule1\tGeneVehicule2\n";

	open(FILEG2M, ">bilan2.gene.vehiculeMC.txt") or die ("Error\n");
	print FILEG2M "Vehicule1\tVehicule2\tGeneVehicule1\tGeneVehicule2\n";	

	######### Gene fasta files
	@FastaGene = `ls $pathGene`;
	my %DerProt;
	my $protD = 0;
	foreach my $fg (@FastaGene){
		
		open(FG, $fg) or die ("Opening Problem $fg");
		while(my $li = <FG>){
			chomp($li);
			
			if($fg =~ m/bacteria/ || $fg =~ m/arch/ || $fg =~ m/plasmid/){
				if($li =~ m/>gi\|[0-9]+\|ref\|(.+)\|.+\[([A-Za-z]+.+)\]$/){
					$idGFnum = $1;
					$idGFcont = $2;
					$protD = 0;
					if($DerProt{$idGFnum}{$protD}){
						$protD=1;
					}
				}
				elsif($li =~ m/^>ref\|.+\|(.+)\|(.+)$/){
					$idGFnum = $2;
					$idGFcont = $1;
					$protD = 0;
					if($DerProt{$idGFnum}{$protD}){
						$protD=1;
					}
				}
				else{
					$hashGeneFasta{$idGFcont}{$idGFnum} .= $li;
					if($protD eq "0"){
						$DerProt{$idGFnum}{$protD} .= $li;
					}
				}
			}
			else{
				if($li =~ m/^>/){
					if(($li =~ m/>(\S+)\.gene_([0-9]+)\D#/) || ($li =~ m/>(\S+)\.([0-9]+)\D#/)){
						$idGFcont = $1;
						$idGFnum = $2;
					}
					elsif(($li =~ m/>(\S+)\.gene_([0-9]+)/) || ($li =~ m/>(\S+)\.([0-9]+)/)){
						$idGFcont = $1;
						$idGFnum = $2;
					}
					elsif($li =~ m/>gi\|\S+\|ref\|(.+)\|.+\[([A-Za-z]+.+)\]$/){
						$idGFnum = $1;
						$idGFcont = $2;
					}
					elsif($li =~ m/^>ref\|.+\|(.+)\|(.+)$/){
						$idGFnum = $2;
						$idGFcont = $1;
					}
				}
				else{
					$hashGeneFasta{$idGFcont}{$idGFnum} .= $li;
				}
			}
		}
		close FG;
	}



	foreach my $cl (keys %hashProtClust){
		++$actualNumCluster;
		open(CLG, ">$cl.gene.fasta") or die ("Creative problem $cl.gene.fasta");
		
		for (my $i = 0 ; $i <= scalar @{$hashProtList{$cl}}-1 ; ++$i){
			print "$i\t@{$hashProtList{$cl}}[$i]\n";
			foreach my $c (keys %{$hashEquiG{${$hashProtList{$cl}}[$i]}}){
				foreach my $p (keys %{$hashEquiG{${$hashProtList{$cl}}[$i]}{$c}}){
					if ($hashGeneFasta{$c}{$p}){
						if(${$hashProtList{$cl}}[$i] =~ s/\|/\.\./){
							print CLG ">${$hashProtList{$cl}}[$i]\n$hashGeneFasta{$c}{$p}\n";
						}
						
					}
					else{
						if($DerProt{$p}{"0"}){
							if(${$hashProtList{$cl}}[$i] =~ s/\|/\.\./){
								print CLG ">${$hashProtList{$cl}}[$i]\n$DerProt{$p}{0}\n";
							}
						}
						else{
							`echo $p >> sedmanque`;
						}
					}
				}
			}
		}
		close CLG;
		
		foreach my $ap (@AliProt){
			if($ap =~ /(Cluster_[0-9]+_[0-9]+).aln$/){
				if($cl eq $1){
					
					open(CA,$ap) or die ("Opening problem $ap");
					open(CAP,">tmp.aln") or die ("Creative problem tmp.aln");
					
					while(my $li = <CA>){
						chomp($li);
						if($li =~ s/\|/\.\./){
							print CAP "$li\n";
						}
						else{
							print CAP "$li\n";
						}
					}
					close CAP;
					close CA;
					#chomp($ap);
					#system("sed 's/|/../' $ap > tmp.aln");
					#system("tranalign -asequence $cl.gene.fasta -bsequence tmp.aln -outseq $cl.gene.aln -table 11");
					system("muscle -quiet -in $cl.gene.fasta -out $cl.gene.aln");
				}			
				
			}
		}
		system("trimal -strictplus -phylip -in $cl.gene.aln -out $cl.gene.aln.phylip >/dev/null 2>&1");
		open(TMP,">seq_temp");
		print TMP "$cl.gene.aln.phylip\nY\n";
		close TMP;
		system("phylip dnadist < seq_temp > temp");
		system("mv outfile $cl.gene.aln.phylip.rendu");
		
			
				
		if($cl=~ m/Cluster_.+_([0-9]+)/){
			my $TotGene = $1;
			if($TotGene > 2){
				system("phyml -d nt -i $cl.gene.aln.phylip");
				
				$tree = $cl.".gene.aln.phylip_phyml_tree.txt";
				
				#Script python
				open(PYT,">patri.py") or die ("Error patri.py");
				print PYT "#! /usr/bin/env python\nimport sys\nimport dendropy\nfrom dendropy import treecalc\n";
				print PYT "arg = sys.argv[1]\nif arg.startswith('-H')or arg.startswith('-h'):\n\tprint(\"python patri.py tree.txt\")\n";
				print PYT "else:\n\ttr = dendropy.Tree.get(path=arg,schema=\"newick\")\n\tdist = dendropy.PhylogeneticDistanceMatrix.from_tree(tr)\n";
				print PYT "\tfor i, t1 in enumerate(tr.taxon_namespace):\n\t\tmatr=t1.label\n\t\tfor j, t2 in enumerate(tr.taxon_namespace):\n\t\t\tmatr+= \"\\t\"+str(dist(t1, t2))\n\t\tprint(matr)\n";

				#Step Patristic
				system("python patri.py $tree > $tree.patristic.txt"); 									
				open(MCG, "$tree.patristic.txt") or die ("Opening Problem $tree.patristic.txt");		
				my $actuGene = 0;
				my $matGene = 0;
				my $gID;
				my %hashMCGene;
				my %hashMCDisGene;
				my %hashMCDisGene2;
					
				while(my $li = <MCG>){
					chomp($li);
					#Recuperation du nombre de genes dans la matrice
					if($li=~ m/^[0-9]+\.\.[0-9]+.+$/){
						my @ligne = split("\t",$li);
						foreach my $c (@ligne){
							if($c =~ m/^[0-9]+\.\.[0-9]+$/){
								++$actuGene;
							}
						}
					}
					if($li =~ m/^([0-9]+\.\.[0-9]+)\t(.+)$/){
						++$matGene;
						$gID = $1;
						my $partd = $2;
						if($matGene <= $actuGene){
							#Recuperation ID gene
							if($gID =~ s/\.\./\|/){
								$hashMCGene{$matGene} = $gID;
							}
							#Recuperation distance gene
							my @liDis = split(/\t/,$partd);
							@{$hashMCDisGene{$gID}} = @liDis;
						}
					}
					if(($actuGene > 6) && ($li =~ m/^\t([0-9.]+\t.+)/)){
						my @liDis2 = split(/\t/,$1);
						#Recuperation distance gene (suite)
						foreach my $d2 (@liDis2){
							push(@{$hashMCDisGene{$hashMCGene{$matGene}}},$d2);
						}
					}
				}
				close MCG;
				my $matrGene = 0;
				foreach my $ggm (keys %hashMCDisGene){
					$matrGene = 0;
					foreach my $d (@{$hashMCDisGene{$ggm}}){
						++$matrGene;
						if($matrGene <= $actuGene){
							if($hashMCGene{$matrGene} =~ m/^([0-9]+)\D+.+$/){
								my $verf = $1;
								if($ggm =~ m/([0-9]+)\|[0-9]+/){
									my $verf2 = $1;

									if(($verf ne $verf2) && ($hashMCGene{$matrGene} ne $ggm)){
										push(@{$hashMCDisGene2{$ggm}},$d);
									}
									else{
										my $distance = 1000000;
										push(@{$hashMCDisGene2{$ggm}},$distance);
									}
								}
							}
						}
					}
				}
				
				foreach my $g (keys %hashMCDisGene){
					for(my $i = 0; $i < scalar(@{$hashMCDisGene{$g}}) ; ++$i){
						my $tab = $i + 1;
						if($g eq $hashMCGene{$tab}){
							if((${$hashMCDisGene{$g}}[$i] =~ m/0\.0+$/) || (${$hashMCDisGene{$g}}[$i] eq "0e+00") || (${$hashMCDisGene{$g}}[$i] eq "0")){
								${$hashMCDisGene{$g}}[$i] = 1000000;
							}
						}
					}
					my $distMin = min @{$hashMCDisGene{$g}};

					for(my $i = 0; $i < scalar(@{$hashMCDisGene{$g}}) ; ++$i){
						#Identification des noeuds et aretes
						my $tab = $i + 1;
						if((${$hashMCDisGene{$g}}[$i] eq $distMin) && ($hashMCGene{$tab} ne $g)){
							if($hashEqui{$g} ne $hashEqui{$hashMCGene{$tab}}){
								if($hashMCG{$hashEqui{$g}}{$hashEqui{$hashMCGene{$tab}}}){
									++$hashMCG{$hashEqui{$g}}{$hashEqui{$hashMCGene{$tab}}};
									print FILEG1M "$hashEqui{$g}\t$hashEqui{$hashMCGene{$tab}}\t$hashPROT{$g}\t$hashPROT{$hashMCGene{$tab}}\n";
								}
								elsif($hashMCG{$hashEqui{$hashMCGene{$tab}}}{$hashEqui{$g}}){
									++$hashMCG{$hashEqui{$g}}{$hashEqui{$hashMCGene{$tab}}};
									print FILEG1M "$hashEqui{$g}\t$hashEqui{$hashMCGene{$tab}}\t$hashPROT{$g}\t$hashPROT{$hashMCGene{$tab}}\n";
								}
								else{
									++$hashMCG{$hashEqui{$g}}{$hashEqui{$hashMCGene{$tab}}};
									print FILEG1M "$hashEqui{$g}\t$hashEqui{$hashMCGene{$tab}}\t$hashPROT{$g}\t$hashPROT{$hashMCGene{$tab}}\n";
								}
							}
						}
					}				
				}	
				
				foreach my $g (keys %hashMCDisGene2){
					for(my $i = 0; $i < scalar(@{$hashMCDisGene2{$g}}) ; ++$i){
						my $tab = $i + 1;
						if($g eq $hashMCGene{$tab}){
							if((${$hashMCDisGene2{$g}}[$i] =~ m/0\.0+$/) || (${$hashMCDisGene2{$g}}[$i] eq "0e+00") || (${$hashMCDisGene2{$g}}[$i] eq "0")){
								${$hashMCDisGene2{$g}}[$i] = 1000000;
							}
						}
					}
					my $distMin = min @{$hashMCDisGene2{$g}};

					for(my $i = 0; $i < scalar(@{$hashMCDisGene2{$g}}) ; ++$i){
						#Identification des noeuds et aretes
						my $tab = $i + 1;
						if((${$hashMCDisGene2{$g}}[$i] eq $distMin) && ($hashMCGene{$tab} ne $g)){
							if($hashEqui{$g} ne $hashEqui{$hashMCGene{$tab}}){
								if($hashMCG2{$hashEqui{$g}}{$hashEqui{$hashMCGene{$tab}}}){
									++$hashMCG2{$hashEqui{$g}}{$hashEqui{$hashMCGene{$tab}}};
									print FILEG2M "$hashEqui{$g}\t$hashEqui{$hashMCGene{$tab}}\t$hashPROT{$g}\t$hashPROT{$hashMCGene{$tab}}\n";
								}
								elsif($hashMCG2{$hashEqui{$hashMCGene{$tab}}}{$hashEqui{$g}}){
									++$hashMCG2{$hashEqui{$g}}{$hashEqui{$hashMCGene{$tab}}};
									print FILEG2M "$hashEqui{$g}\t$hashEqui{$hashMCGene{$tab}}\t$hashPROT{$g}\t$hashPROT{$hashMCGene{$tab}}\n";
								}
								else{
									++$hashMCG2{$hashEqui{$g}}{$hashEqui{$hashMCGene{$tab}}};
									print FILEG2M "$hashEqui{$g}\t$hashEqui{$hashMCGene{$tab}}\t$hashPROT{$g}\t$hashPROT{$hashMCGene{$tab}}\n";
								}
							}
						}
					}				
				}
			}
		} 
		#}
		
		
		
		#Creation fichier Network contenant les differentes aretes du reseau
		open(NETGMC,">Network.gene.vehiculeMC.csv") or die ("Creative problem Network.gene.vehiculeMC.csv");
		print NETGMC "Source\tTarget\tWeight\tType\n";

		foreach my $u (sort keys %hashMCG){
			foreach my $v (sort keys %{$hashMCG{$u}}){
				my $weight = 1/$hashMCG{$u}{$v};
				$g4->add_weighted_edge($u, $v, $weight);
				print NETGMC "$u\t$v\t$weight\tUndirected\n";
			}	
		}
		close NETGMC;

		#Creation du fichier Table contenant les informations des différents noeuds du reseau
		open(TABGMC,">Table.gene.vehiculeMC.csv") or die ("Creative problem Table.gene.vehiculeMC.csv");
		print TABGMC "ID\tDescription\n";

		my @Z = $g4->vertices;
		foreach my $z(@Z){
			if($z =~ m/^([0-9]+)_.+/){
				print TABGMC "$z\t$1\n";
			}
		}
		close TABGMC;
		
		
		
		#Creation fichier Network contenant les differentes aretes du reseau
		open(NETGMC2,">Network2.gene.vehiculeMC.csv") or die ("Creative problem Network2.gene.vehiculeMC.csv");
		print NETGMC2 "Source\tTarget\tWeight\tType\n";

		foreach my $u (sort keys %hashMCG2){
			foreach my $v (sort keys %{$hashMCG2{$u}}){
				my $weight = 1/$hashMCG2{$u}{$v};
				$g5->add_weighted_edge($u, $v, $weight);
				print NETGMC2 "$u\t$v\t$weight\tUndirected\n";
			}	
		}
		close NETGMC2;

		#Creation du fichier Table contenant les informations des différents noeuds du reseau
		open(TABGMC2,">Table2.gene.vehiculeMC.csv") or die ("Creative problem Table2.gene.vehiculeMC.csv");
		print TABGMC2 "ID\tDescription\n";

		my @A = $g5->vertices;
		foreach my $a(@A){
			if($a =~ m/^([0-9]+)_.+/){
				print TABGMC2 "$a\t$1\n";
			}
		}
		close TABGMC2;
		
		print("\rProccessing cluster $actualNumCluster of $NbClustBon\n");
	}
	
	@clMatrG =`ls *gene.aln.phylip.rendu`;
	
	my %hashG;
	my %hashG2;
	
	foreach my $clMG (@clMatrG){
		chomp($clMG);
		my $cl = "";
		if($clMG =~ m/(Cluster_[0-9]+_[0-9]+)\.gene\.aln\.phylip\.rendu/){
			$cl = $1;
			
			if($clBon{$cl}){
				
				open(CGM, $clMG) or die ("Opening Problem $clMG");

				my %hashGene;
				my %hashDisGene;
				my %hashDisGene2;
				my $nbgene = 0;
				my $numgene = 0;
				my $matrGene = 0;
				my $idgene; 
				
				while(my $li = <CGM>){
					chomp ($li);
					#Recuperation du nombre de genes dans la matrice
					if($li=~ m/^\s+([0-9]+)$/){
						$nbgene = $1;
					}
					if($li =~ m/^([0-9]+\.\.[0-9]+)\s+([0-9.\s]+)$/){
						++$numgene;
						$idgene = $1;
						my $partd = $2;
						
						#Recuperation ID gene
						if($idgene =~ s/\.\./\|/){
							$hashGene{$numgene} = $idgene;
						}
						#Recuperation distance gene (suite)
						my @liDis = split(/\s/,$partd);
						@{$hashDisGene{$idgene}} = @liDis;
					}
					if(($nbgene > 6) && ($li =~ m/^\s\s([0-9.]+\s.+)/)){
						my @liDis2 = split(/\s/,$1);
						#Recuperation distance gene (suite)
						foreach my $d2 (@liDis2){
							push(@{$hashDisGene{$hashGene{$numgene}}},$d2);
						}
					}
				}
				close CGM;
				
				foreach my $gg (keys %hashDisGene){
					$matrGene = 0;
					foreach my $d (@{$hashDisGene{$gg}}){
						++$matrGene;
						if($hashGene{$matrGene} =~ m/([0-9]+)\|[0-9]+/){
							my $verf = $1;
							if($gg =~ m/([0-9]+)\|[0-9]+/){
								my $verf2 = $1;

								if(($verf ne $verf2) && ($hashGene{$matrGene} ne $gg)){
									push(@{$hashDisGene2{$gg}},$d);
								}
								else{
									my $distance = 1000000;
									push(@{$hashDisGene2{$gg}},$distance);
								}
							}
						}
					}
				}
				
				
				foreach my $g (keys %hashDisGene){
					for(my $i = 0; $i < scalar(@{$hashDisGene{$g}}) ; ++$i){
						my $tab = $i + 1;
						if($g eq $hashGene{$tab}){
							if(${$hashDisGene{$g}}[$i] =~ m/0\.00+$/){
								${$hashDisGene{$g}}[$i] = 1000000;
							}
						}
					}
					my $disMin = min @{$hashDisGene{$g}};

					for(my $i = 0; $i < scalar(@{$hashDisGene{$g}}) ; ++$i){
						my $tab = $i + 1;					
						#Identification des noeuds et aretes
						if((${$hashDisGene{$g}}[$i] eq $disMin) && ($hashGene{$tab} ne $g)){
							if($hashEqui{$g} ne $hashEqui{$hashGene{$tab}}){
								if($hashG{$hashEqui{$g}}{$hashEqui{$hashGene{$tab}}}){
									++$hashG{$hashEqui{$g}}{$hashEqui{$hashGene{$tab}}};
									print FILEG1 "$hashEqui{$g}\t$hashEqui{$hashGene{$tab}}\t$hashPROT{$g}\t$hashPROT{$hashGene{$tab}}\n";
								}
								elsif($hashG{$hashEqui{$hashGene{$tab}}}{$hashEqui{$g}}){
									++$hashG{$hashEqui{$g}}{$hashEqui{$hashGene{$tab}}};
									print FILEG1 "$hashEqui{$g}\t$hashEqui{$hashGene{$tab}}\t$hashPROT{$g}\t$hashPROT{$hashGene{$tab}}\n";
								}
								else{
									++$hashG{$hashEqui{$g}}{$hashEqui{$hashGene{$tab}}};
									print FILEG1 "$hashEqui{$g}\t$hashEqui{$hashGene{$tab}}\t$hashPROT{$g}\t$hashPROT{$hashGene{$tab}}\n";
								}
							}
						}
					}				
				}
				foreach my $g2 (keys %hashDisGene2){
					for(my $i = 0; $i < scalar(@{$hashDisGene2{$g2}}) ; ++$i){
						my $tab = $i + 1;
						if($g2 eq $hashGene{$tab}){
							if(${$hashDisGene2{$g2}}[$i] =~ m/0\.00+$/){
								${$hashDisGene2{$g2}}[$i] = 1000000;
							}
						}
					}
					my $distaMin = min @{$hashDisGene2{$g2}};

					for(my $i = 0; $i < scalar(@{$hashDisGene2{$g2}}) ; ++$i){
						my $tab = $i + 1;					
						#Identification des noeuds et aretes
						if((${$hashDisGene2{$g2}}[$i] eq $distaMin) && ($hashGene{$tab} ne $g2)){
							if($hashEqui{$g2} ne $hashEqui{$hashGene{$tab}}){
								if($hashG2{$hashEqui{$g2}}{$hashEqui{$hashGene{$tab}}}){
									++$hashG2{$hashEqui{$g2}}{$hashEqui{$hashGene{$tab}}};
									print FILEG2 "$hashEqui{$g2}\t$hashEqui{$hashGene{$tab}}\t$hashPROT{$g2}\t$hashPROT{$hashGene{$tab}}\n";
								}
								elsif($hashG2{$hashEqui{$hashGene{$tab}}}{$hashEqui{$g2}}){
									++$hashG2{$hashEqui{$g2}}{$hashEqui{$hashGene{$tab}}};
									print FILEG2 "$hashEqui{$g2}\t$hashEqui{$hashGene{$tab}}\t$hashPROT{$g2}\t$hashPROT{$hashGene{$tab}}\n";
								}
								else{
									++$hashG2{$hashEqui{$g2}}{$hashEqui{$hashGene{$tab}}};
									print FILEG2 "$hashEqui{$g2}\t$hashEqui{$hashGene{$tab}}\t$hashPROT{$g2}\t$hashPROT{$hashGene{$tab}}\n";
								}
							}
						}
					}				
				}
			}
		}
	}
	close FILEG1M;
	close FILEG2M;
	close FILEG1;
	close FILEG2;
	
	print "Output files creation : ";
	
	#Creation fichier Network contenant les differentes aretes du reseau
	open(NETGH,">Network.gene.vehicule.csv") or die ("Creative problem Network.gene.vehicule.csv");
	print NETGH "Source\tTarget\tWeight\tType\n";

	foreach my $u (sort keys %hashG){
		foreach my $v (sort keys %{$hashG{$u}}){
			my $weight = 1/$hashG{$u}{$v};
			$g0->add_weighted_edge($u, $v, $weight);
			print NETGH "$u\t$v\t$weight\tUndirected\n";
			
		}	
	}
	close NETGH;

	#Creation du fichier Table contenant les informations des différents noeuds du reseau
	open(TABGM,">Table.gene.vehicule.csv") or die ("Creative problem Table.gene.vehicule.csv");
	print TABGM "ID\tDescription\n";

	my @V = $g0->vertices;
	foreach my $v(@V){
		if($v =~ m/^([0-9]+)_.+/){
			print TABGM "$v\t$1\n";
		}
	}
	close TABGM;


	#Creation fichier Network contenant les differentes aretes du reseau
	open(NETGH2,">Network2.gene.vehicule.csv") or die ("Creative problem Network2.gene.vehicule.csv");
	print NETGH2 "Source\tTarget\tWeight\tType\n";

	foreach my $u (sort keys %hashG2){
		foreach my $v (sort keys %{$hashG2{$u}}){
			my $weight = 1/$hashG2{$u}{$v};
			$g2->add_weighted_edge($u, $v, $weight);
			print NETGH2 "$u\t$v\t$weight\tUndirected\n";
			
		}	
	}
	close NETGH2;

	#Creation du fichier Table contenant les informations des différents noeuds du reseau
	open(TABGM2,">Table2.gene.vehicule.csv") or die ("Creative problem Table2.gene.vehicule.csv");
	print TABGM2 "ID\tDescription\n";

	my @X = $g2->vertices;
	foreach my $x(@X){
		if($x =~ m/^([0-9]+)_.+/){
			print TABGM2 "$x\t$1\n";
		}
	}
	close TABGM2;
	print "done\n";
	
	print "Vehicle analysis with genes : done\n";
	
}





#########Vehicle analysis with proteins

if($seqTY eq "aa"){
		
	print "Vehicle analysis with proteins : \n";
	
	my %hashG;
	my %hashG2;
	open(FILE1, ">bilan.prot.vehicule.txt") or die ("Error\n");
	print FILE1 "Vehicule1\tVehicule2\tProteinVehicule1\tProteinVehicule2\n";
	
	open(FILE2, ">bilan2.prot.vehicule.txt") or die ("Error\n");
        print FILE2 "Vehicule1\tVehicule2\tProteinVehicule1\tProteinVehicule2\n";
	
	foreach my $clM (keys %clBon){
		my $cl = "";
			
				$cl = "$pathProtM"."$clM".".aln.phylip.rendu";
				open(CPM, $cl) or die ("Opening Problem $cl");
				
				my %hashPro;
				my %hashDisPro;
				my %hashDisPro2;
				my $pID = "";
				my $nbpro = 0;
				my $numpro = 0;
				my $matrPro = 0;
				my $length = "";
				
				while(my $li = <CPM>){
					
					chomp ($li);
					#Recuperation du nombre de proteines dans la matrice
					if($li=~ m/^\s+([0-9]+)$/){
						$nbpro = $1;
					}
					if($li =~ m/^([0-9]+\|[0-9]+)\s+([0-9.\-]+.+)$/){
						++$numpro;
						#Recuperation ID proteine
						$pID = $1;
						$hashPro{$numpro} = $pID;

						#Recuperation distance proteine (suite)
						my @liDis = split(/\s+/,$2);
						@{$hashDisPro{$pID}} = @liDis;
					}
					if(($nbpro > 6) && ($li =~ m/^\s\s([0-9.\-]+\s\s.+)/)){
						my @liDis2 = split(/\s+/,$1);
						
						#Recuperation distance proteine (suite)
						foreach my $d2 (@liDis2){
							push(@{$hashDisPro{$hashPro{$numpro}}},$d2);
						}
					}
					
				}
				close CPM;
				
				foreach my $pp (keys %hashDisPro){
					$matrPro = 0;
					foreach my $d (@{$hashDisPro{$pp}}){
						++$matrPro;
						if($hashPro{$matrPro} =~ m/([0-9]+)\|[0-9]+/){
							my $verf = $1;
							if($pp =~ m/([0-9]+)\|[0-9]+/){
								my $verf2 = $1;

								if(($verf ne $verf2) && ($hashPro{$matrPro} ne $pp)){
									push(@{$hashDisPro2{$pp}},$d);
								}
								else{
									my $distance = 1000000;
									push(@{$hashDisPro2{$pp}},$distance);
								}
							}
						}
					}
				}
				
				
				foreach my $p (keys %hashDisPro){
					for(my $i = 0; $i < scalar(@{$hashDisPro{$p}}) ; ++$i){
						my $tab = $i + 1;
						if($p eq $hashPro{$tab}){
							if(${$hashDisPro{$p}}[$i] =~ m/0\.00+$/){
								${$hashDisPro{$p}}[$i] = 1000000;
							}
						}
					}
					my $disMin = min @{$hashDisPro{$p}};

					for(my $i = 0; $i < scalar(@{$hashDisPro{$p}}) ; ++$i){
						my $tab = $i + 1;
						
						#Identification des noeuds et aretes
						if((${$hashDisPro{$p}}[$i] eq $disMin) && ($hashPro{$tab} ne $p)){
							if($hashEqui{$p} ne $hashEqui{$hashPro{$tab}}){
								if($hashG{$hashEqui{$p}}{$hashEqui{$hashPro{$tab}}}){
									++$hashG{$hashEqui{$p}}{$hashEqui{$hashPro{$tab}}};
									print FILE1 "$hashEqui{$p}\t$hashEqui{$hashPro{$tab}}\t$hashPROT{$p}\t$hashPROT{$hashPro{$tab}}\n";
								}
								elsif($hashG{$hashEqui{$hashPro{$tab}}}{$hashEqui{$p}}){
									++$hashG{$hashEqui{$p}}{$hashEqui{$hashPro{$tab}}};
									print FILE1 "$hashEqui{$p}\t$hashEqui{$hashPro{$tab}}\t$hashPROT{$p}\t$hashPROT{$hashPro{$tab}}\n";
								}
								else{
									++$hashG{$hashEqui{$p}}{$hashEqui{$hashPro{$tab}}};
									print FILE1 "$hashEqui{$p}\t$hashEqui{$hashPro{$tab}}\t$hashPROT{$p}\t$hashPROT{$hashPro{$tab}}\n";
								}
							}
						}
					}				
				}
				
				
				
				foreach my $p2 (keys %hashDisPro2){

					for(my $i = 0; $i < scalar(@{$hashDisPro2{$p2}}) ; ++$i){
						my $tab = $i + 1;
						if($p2 eq $hashPro{$tab}){
							if(${$hashDisPro2{$p2}}[$i] =~ m/0\.00+$/){
								${$hashDisPro2{$p2}}[$i] = 1000000;
							}
						}
					}
					my $disMin2 = min @{$hashDisPro2{$p2}};

					for(my $i = 0; $i < scalar(@{$hashDisPro2{$p2}}) ; ++$i){
						my $tab = $i + 1;
						
						#Identification des noeuds et aretes
						if((${$hashDisPro2{$p2}}[$i] eq $disMin2) && ($hashPro{$tab} ne $p2)){
							if($hashEqui{$p2} ne $hashEqui{$hashPro{$tab}}){
								if($hashG2{$hashEqui{$p2}}{$hashEqui{$hashPro{$tab}}}){
									++$hashG2{$hashEqui{$p2}}{$hashEqui{$hashPro{$tab}}};
									print FILE2 "$hashEqui{$p2}\t$hashEqui{$hashPro{$tab}}\t$hashPROT{$p2}\t$hashPROT{$hashPro{$tab}}\n";
								}
								elsif($hashG2{$hashEqui{$hashPro{$tab}}}{$hashEqui{$p2}}){
									++$hashG2{$hashEqui{$p2}}{$hashEqui{$hashPro{$tab}}};
									print FILE2 "$hashEqui{$p2}\t$hashEqui{$hashPro{$tab}}\t$hashPROT{$p2}\t$hashPROT{$hashPro{$tab}}\n";
								}
								else{
									++$hashG2{$hashEqui{$p2}}{$hashEqui{$hashPro{$tab}}};
									print FILE2 "$hashEqui{$p2}\t$hashEqui{$hashPro{$tab}}\t$hashPROT{$p2}\t$hashPROT{$hashPro{$tab}}\n";
								}
							}
						}
					}				
				}
				
			#}
		#}
	}
	close FILE1;
	close FILE2;
	print "\tOutput files creation : ";
	
	#Creation fichier Network contenant les differentes aretes du reseau
	open(NETH,">Network.prot.vehicule.csv") or die ("Creative problem Network.prot.vehicule.csv");
	print NETH "Source\tTarget\tWeight\tType\n";

	foreach my $u (sort keys %hashG){
		foreach my $v (sort keys %{$hashG{$u}}){
			my $weight = "";
			if($hashG{$u}{$v}=~ m/[0-9]+/){
				$weight = 1/$hashG{$u}{$v};
			}
			else{$weight = 0;}
			
			$g1->add_weighted_edge($u, $v, $weight);
			if($v =~ m/([0-9]+)_.+\[([A-Za-z0-9_\s.\/:-]+)$/){
				$v = $1."_".$2;
			}
			if($u =~ m/([0-9]+)_.+\[([A-Za-z0-9_\s.\/:-]+)$/){
				$u = $1."_".$2;
			}
			print NETH "$u\t$v\t$weight\tUndirected\n";
		}	
	}
	close NETH;

	#Creation du fichier Table contenant les informations des différents noeuds du reseau
	open(TABM,">Table.prot.vehicule.csv") or die ("Creative problem Table.prot.vehicule.csv");
	print TABM "ID\tDescription\n";

	my @W = $g1->vertices;
	foreach my $w(@W){
		if($w =~ m/^([0-9]+)_.+/){
			print TABM "$w\t$1\n";
		}
	}
	close TABM;


	#Creation fichier Network contenant les differentes aretes du reseau
	open(NETH2,">Network2.prot.vehicule.csv") or die ("Creative problem Network2.prot.vehicule.csv");
	print NETH2 "Source\tTarget\tWeight\tType\n";

	foreach my $u (sort keys %hashG2){
		foreach my $v (sort keys %{$hashG2{$u}}){
			my $weight = "";
			if($hashG2{$u}{$v}=~ m/[0-9]+/){
				$weight = 1/$hashG2{$u}{$v};
			}
			else{$weight = 0;}
			$g3->add_weighted_edge($u, $v, $weight);
			if($v =~ m/([0-9]+)_.+\[([A-Za-z0-9_\s.\/:-=]+)$/){
				$v = $1."_".$2;
			}
			if($u =~ m/([0-9]+)_.+\[([A-Za-z0-9_\s.\/:-=]+)$/){
				$u = $1."_".$2;
			}
			print NETH2 "$u\t$v\t$weight\tUndirected\n";
			
		}	
	}
	close NETH2;

	#Creation du fichier Table contenant les informations des différents noeuds du reseau
	open(TABM2,">Table2.prot.vehicule.csv") or die ("Creative problem Table2.prot.vehicule.csv");
	print TABM2 "ID\tDescription\n";

	my @Y = $g3->vertices;
	foreach my $y(@Y){
		if($y =~ m/^([0-9]+)_.+/){
			print TABM2 "$y\t$1\n";
		}
	}
	close TABM2;
	print "done\n";
	
	#print "Vehicle analysis with proteins : done\n";
	
}




if($clean eq "yes" || "Y" || "y" || "Yes" || "YES")
{
	print "Cleaning : ";
	
	if($seqTY eq "nt"){
		if( -e "tmp.aln"){
			system("rm tmp.aln");
		}
		system("mkdir $dir/ $dir/HgC/");
		system("mv *.csv $dir/");
		system("mv *.gene.fasta $dir/HgC/");
		system("rm seq_temp temp");
		system("mkdir $dir/phylip/ $dir/tranalign/ $dir/trimal/ $dir/phyml/ $dir/cophenetic/");
		system("mv *.aln $dir/tranalign/");
		system("mv *.phylip $dir/trimal/");
		system("mv *.rendu $dir/phylip/");
		system("mv *cophenetic.r $dir/cophenetic/");
		system("mv *_phyml_*.txt $dir/phyml/");
	}
	if($seqTY eq "aa"){
		system("mkdir $dir/");
		system("mv *.vehicule.csv $dir/");
	}
	system("mv bilan* $dir/");
	print "done\n";

}




}



