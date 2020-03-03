#!/usr/bin/perl -w
use strict;
use threads;
#	Execution en local de Inparanoid entre chaque paire de genome 
#	Le script Inparanoid a ete modifie afin qu'il n'execute pas les BLAST 
# mais qu'il utilise les resultas BLAST executes au prealable sur le cluster

#perl 2_inparanoid.pl [repertoire de travail] [repertoire contenant les genomes] [repertoire contenant les resultats blasts]
#perl 2_inparanoid.pl /home/cecile/Bureau/miseajourBD/local/Inparanoid/inparanoid_cecile/

#############
# Declaration
#############

# Parametres
my $dirWork = $ARGV[0]; # Repertoire de travail ou se situe les scripts
my $dirInGenome = $ARGV[1]; # INPUT : Repertoire contenant les genomes
my $dirBlast = $ARGV[2]; # INPUT : repertoire contenant les resultats BLAST realises au prealable (etape 1)
my $dirOutG=$ARGV[3];
my $dirOut = $dirOutG.'/sqltable/';  # OUTPUT : repertoire contenant tous les groupes d'orthologues obtenus avec Inparanoid (groupe entre 2 genomes uniquement)
my $dirTmp = $dirOutG.'/tmp/'; # OUTPUT : repertoire contenant tous les fichiers FASTA (car ils doivent etre dans un meme repertoire pour etre utilise par inparanoid)
my $dirInparanoid = $dirWork.'inparanoid_4.1'; # INPUT : repertoire contenant les scripts d'Inparanoid
my $threadsAllowed=$ARGV[4];
# Autres variables
my @OldGenome = ();
my @NewGenome = ();
my @File = ();
my $esp1 = ''; #nom de l'espece1
my $esp2 = ''; #nom de l'espece2
my $file = '';
my @jobs=();

############
# Programme
###########

# Creation des repertoires OUTPUT
mkdir($dirTmp) unless (-e $dirTmp);
mkdir($dirOut) unless (-e $dirOut);

# Recuperation des noms des genomes
if (-e $dirInGenome){
	chdir($dirInGenome);
	@File = glob("*fasta");
	foreach $file (@File) {
		$file =~ s/.fasta$//;
		push(@NewGenome,$file);
	}
}
else{
	die ("Error : $dirInGenome not found\n");
}

# Execution d'Inparanoid entre chaque paire de genomes a moins que le resultat ait deja ete cree
foreach $esp1 (@NewGenome) {
	print "esp1 $esp1 \n";
	foreach $esp2 (@NewGenome) {
	   print "esp2 $esp2 \n";
	   push(@jobs,threads->create(sub{
		unless ($esp1 eq $esp2) {
			print $esp1." ".$esp2."\n"; 
			#launch_ip($esp1,$esp2,$dirBlast,$dirTmp,$dirOut,$dirWork,$dirInparanoid) unless ((-e $dirOut.'sqltable.'.$esp1.'-'.$esp2) || (-e $dirOut.'sqltable.'.$esp2.'-'.$esp1));
			print "perl $dirInparanoid/inparanoid_blastp.pl $esp1.fasta $esp2.fasta $dirBlast $dirOut";
			system("perl $dirInparanoid/inparanoid_blastp.pl $esp1.fasta $esp2.fasta $dirBlast $dirOut");
		}
	   }));
	   sleep(1) while(threads->list(threads::running)>=$threadsAllowed);
	}
}
$_->join for @jobs;
# Suppression du repertoire $dirTmp
chdir($dirTmp);
@File = glob("*");
foreach $file (@File) {
	unlink($dirTmp.$file);
}
rmdir($dirTmp);

