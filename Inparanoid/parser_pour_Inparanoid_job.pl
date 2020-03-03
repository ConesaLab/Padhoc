#!/usr/bin/perl -w
#use threads;
# Description:
# parser les sorties de blast de maniere a avoir des fichiers au bon format pour inparanoid
# format tabulé avec les informations suivantes:
# * Query id
# * Hit id
# * Bit score
# * Query length
# * Hit length
# * Length of longest segment on query. This is defined as the length of the segment from the
#  	first position od the first hsp to the last position of the last hsp. I.e. if the
# 	hsps are 1-10 and 90-100, the length of the longest segment will be 100.
# * Length of longest segment on hit, see explanation above.
# * Total match length on query. This is defined as the total length of all hsps. If the hsps
# 	are 1-10 and 90-100, the length will be 20.
# * Total match length on hit. see explanation above
# * Positions for all segments on the query and on the hit. Positions on the query is written
# 	as p:1-100 and positions on the hit is writteh as h:1-100. Positions for all hsps are
#	  specified separated nt tab. Positions for query and hut for a segment is separated by 
# 	one whitespace. 

#score_cutoff utilisé par default par inparanoid: 40
#format blast d'entree : tabulé
#qseqid0 sseqid1 pident2 length3 mismatch4 gapopen5 qstart6 qend7 sstart8 send9 evalue10 bitscore11 qlen12 slen13

#17/06/2014
#perl parser_pour_Inparanoid.pl [score_cutoff] [dossier contenant les blasts] [dossier de sortie]
#perl parser_pour_Inparanoid.pl 40 /home/cecile/Bureau/miseajourBD/local/result/blastp/ /home/cecile/Bureau/miseajourBD/local/Inparanoid/inparanoid_cecile/result/blast/
#perl parser_pour_Inparanoid.pl 40 /home/cecile/Bureau/miseajourBD/local/Inparanoid/inparanoid_cecile/test/ /home/cecile/Bureau/miseajourBD/local/Inparanoid/inparanoid_cecile/test/

########################################
#Déclaration
############
$score_cutoff=$ARGV[0];
$dirBlast=$ARGV[1];
$dirOut=$ARGV[2];
$threadsAllowed=$ARGV[3];
$f=$ARGV[4];

$sumbs=0;#somme de tous les bits scores des hsp pour une query et un subject
$minqstart=1.e1000;#j'ai initialisé avec un grand nombre
$minsstart=1.e1000;#j'ai initialisé avec un grand nombre
#print "min=".$minsstart."\n";
$maxqend=0;
$maxsend=0;
$tmlq=0;#sum(pour chaque hsp(qend -qstart +1)): Total match length on query. 
$tmls=0;#sum(pour chaque hsp(send -sstart +1)): Total match length on hit.
#$deshsp="";#pour_chaque_hsp_ecrire(q:[qstart]-[qend] s:[sstart]-[send]\t)
$qlen=0;
$slen=0;
$lsq=0;#max([qend])-min([qstart])+1
$lss=0;#max([send])-min([sstart])+1
%vide=();
%hsps=();
@jobs=();

#mkdir($dirOut)unless(-e $dirOut);

#on recupere tous les nom des fichiers blasts
#opendir(BL,$dirBlast);
#@fb=readdir(BL);
#pour chaque fichier blast
#foreach $f(@fb){

unless(($f eq "." )||($f eq "..")){
  # push(@jobs,threads->create(sub{
	print "$f parsage commencé ";
	system("date");
	open(IN,$dirBlast.$f);
	$nomaparser=$f;
	#if($f=~/\d/){#modif par cecile pour les nom de blast avec plus de un _
	#	$nomaparser=~s/_(\d)/\.fa-$1/;
	#}
	#else{
	#	$nomaparser=~s/_/\.fa-/;#modif par cecile pour les nom de blast avec plus de un _
	#	$nomaparser=~s/\.blast$/\.fa/;
	#}
	#$nomaparser=~s/\.blast/\.fa/;
	$nomaparser=~s/\.blast$//;
	$nomaparser=~s/\.fa\.fa/\.fa/g;
	open(OUT,">".$dirOut.$nomaparser)or die($dirOut.$nomaparser);#fichier de sortie
	$queryprecedent="";
	$subjectprecedent="";
	while(<IN>){
		#print $_;
		chomp;
		unless(/^#/){
				@l=split(/\t/,$_);
				#si la query et l'id courant sont identiques à la query et l'id precedent
				if(($l[0] eq $queryprecedent) && ($l[1] eq $subjectprecedent)){
					$overlap=0;#1 si plus de 5% de recouvrement entre deux hsp
					foreach $v(keys(%hsps)){
						$hsps{$v}=~/q:(.+)-(.+) h:(.+)-(.+)/;
						$overlap=check_if_overlap($1,$2,$3,$4,$l[6],$l[7],$l[8],$l[9]);#	my $start_hsp1_query  $end_hsp1_query $start_hsp1_hit $end_hsp1_hit $start_hsp2_query $end_hsp2_query $start_hsp2_hit $end_hsp2_hit 
						if($overlap){
							#print "overlap > a 5%: $l[0] - $l[1] hsp1 $1,$2,$3,$4  hsp2 $l[6],$l[7],$l[8],$l[9]\n";
							last;
						}
					}
					#a moins qu'il y ai plus de 5% de recouvrement entre ce nouvel hsp et au autre deja dans la table de hashage
					unless($overlap){
						#mise a jour du min de qstart
						if($l[6]<$minqstart){
							$minqstart=$l[6];
						}
						#mise a jour du min de sstart
						if($l[8]<$minsstart){
							$minsstart=$l[8];
						}
						#mise a jour du max de qend
						if($l[7]>$maxqend){
							$maxqend=$l[7];
						}
						#mise a jour du max de send
						if($l[9]>$maxsend){
							$maxsend=$l[9];
						}
						$tmlq+=$l[7]-$l[6]+1;#sum(pour chaque hsp(qend -qstart +1)): Total match length on query. 
						$tmls+=$l[9]-$l[8]+1;#sum(pour chaque hsp(send -sstart +1)): Total match length on hit.
						#je place les descriptions des hsps dans une table de hashage afin de pouvoir les triee en fonction de la position du debut du hsp
						if(exists($hsps{$l[6]})){#si un autre hsp commence en cette position
							$hsps{$l[6]}.="q:$l[6]-$l[7] h:$l[8]-$l[9]";
						}
						else{
							$hsps{$l[6]}="q:$l[6]-$l[7] h:$l[8]-$l[9]";
						}
						#$deshsp.="\tq:$l[6]-$l[7] s:$l[8]-$l[9]";#pour_chaque_hsp_ecrire(q:[qstart]-[qend] s:[sstart]-[send]\t)
						$sumbs+=$l[11];
					}
				}
				else{
					#ecrire la ligne correspondant a la query et a la seq dans le fichier de sortie
					unless($queryprecedent eq ""){
						unless($sumbs<$score_cutoff){#on n'ecrit pas les hits avec un score inférieur au score_cutoff
							#print "$maxqend-$minqstart+1\n";
							$lsq=$maxqend-$minqstart+1;#max([qend])-min([qstart])+1
							$lss=$maxsend-$minsstart+1;#max([send])-min([sstart])+1
							#print "$maxsend-$minsstart+1\n";
							
							#print "\n$queryprecedent\t$subjectprecedent\t$sumbs\t$qlen\t$slen\t$lsq\t$lss\t$tmlq\t$tmls";
							print OUT "$queryprecedent\t$subjectprecedent\t$sumbs\t$qlen\t$slen\t$lsq\t$lss\t$tmlq\t$tmls";
							foreach $val(sort{$a<=>$b}keys(%hsps)){
								print OUT "\t$hsps{$val}";
								#print "\t$hsps{$val}";
							}
							print OUT "\n";
						}
					}
					#reinitialisation
					%hsps=%vide;
					$qlen=$l[12];
					$slen=$l[13];
					$sumbs=$l[11];
					$queryprecedent=$l[0];
					$subjectprecedent=$l[1];
					$minqstart=$l[6];
					$minsstart=$l[8];
					$maxqend=$l[7];
					$maxsend=$l[9];
					$tmlq=$l[7]-$l[6]+1;#sum(pour chaque hsp(qend -qstart +1)): Total match length on query. 
					$tmls=$l[9]-$l[8]+1;#sum(pour chaque hsp(send -sstart +1)): Total match length on hit.
					#$deshsp="q:$l[6]-$l[7] s:$l[8]-$l[9]";#pour_chaque_hsp_ecrire(q:[qstart]-[qend] s:[sstart]-[send]\t)
					$hsps{$l[6]}="q:$l[6]-$l[7] h:$l[8]-$l[9]";
				}
		}
	}
	#ecriture de la derniere ligne
	unless($sumbs<$score_cutoff){
		$lsq=$maxqend-$minqstart+1;#max([qend])-min([qstart])+1
		$lss=$maxsend-$minsstart+1;#max([send])-min([sstart])+1
		#print "\n$queryprecedent\t$subjectprecedent\t$sumbs\t$qlen\t$slen\t$lsq\t$lss\t$tmlq\t$tmls";
		print OUT "$queryprecedent\t$subjectprecedent\t$sumbs\t$qlen\t$slen\t$lsq\t$lss\t$tmlq\t$tmls";
		foreach $val(sort{$a<=>$b}keys(%hsps)){
			print OUT "\t$hsps{$val}";
			#print "\t$hsps{$val}";
		}
		print OUT "\n";
	}
	#print "\n";
	close(OUT);
	close(IN);
	print " et terminé à ";
	system("date");
	print "\n";
    #}));
    #sleep(1) while(threads->list(threads::running)>=$threadsAllowed);
 }
#}
#closedir(BL);
#$_->join for @jobs;

# When this subroutine is used, non-linear arrangements of the hsps are not allowed.
sub check_if_overlap {
	my $start_hsp1_query  = shift;	# One hsp
	my $end_hsp1_query=shift;
	my $start_hsp1_hit  = shift;	
	my $end_hsp1_hit=shift;
  my $start_hsp2_query = shift;	# Another hsp
 	my $end_hsp2_query=shift; 
 	my $start_hsp2_hit = shift;	
 	my $end_hsp2_hit=shift; 
 	
# Check which hsp is N-teminal.
  if ($start_hsp1_query eq $start_hsp2_query) { # If the fragments start at the same position, there is an overlap (100% of shorter sequence)
			return 1;
  } elsif ($start_hsp1_query < $start_hsp2_query) { 
			$first_hsp_start_q = $start_hsp1_query;
			$first_hsp_stop_q = $end_hsp1_query;
			$last_hsp_start_q = $start_hsp2_query;
			$last_hsp_stop_q = $end_hsp2_query;
			$first_hsp_start_h = $start_hsp1_hit;
			$first_hsp_stop_h = $end_hsp1_hit;
			$last_hsp_start_h = $start_hsp2_hit;
			$last_hsp_stop_h = $end_hsp2_hit;
	} else {
			$first_hsp_start_q = $start_hsp2_query;
			$first_hsp_stop_q = $end_hsp2_query;
			$last_hsp_start_q = $start_hsp1_query;
			$last_hsp_stop_q = $end_hsp1_query;
			$first_hsp_start_h = $start_hsp2_hit;
			$first_hsp_stop_h = $end_hsp2_hit;
			$last_hsp_start_h = $start_hsp1_hit;
			$last_hsp_stop_h = $end_hsp1_hit;
	}

	# Return whether there is an overlap or not.
	return (_check_overlap_linear($first_hsp_start_q, $first_hsp_stop_q, $last_hsp_start_q, $last_hsp_stop_q) || _check_overlap_linear($first_hsp_start_h, $first_hsp_stop_h, $last_hsp_start_h, $last_hsp_stop_h));
#	return (_check_overlap_linear($first_hsp->[1], $first_hsp->[2], $last_hsp->[1], $last_hsp->[2]) || _check_overlap_linear($first_hsp->[3], $first_hsp->[4], $last_hsp->[3], $last_hsp->[4]));
}

sub _check_overlap_linear {
	my ($start1, $end1, $start2, $end2) = @_;
# Length of segment 1
  my $length1 = $end1 - $start1 + 1;
# Length of segment 2
	my $length2 = $end2 - $start2 + 1;
 # Get the length of the sortest of these segments
  my $shortest_length = ($length1 < $length2)?$length1:$length2;
 # Maxumin of 5% overlap (witg regard to the shorter segment) is allowed
 	#print "($start2 - $end1 - 1) / $shortest_length\n";
  return (($start2 - $end1 - 1) / $shortest_length < - 0.05);  
}



