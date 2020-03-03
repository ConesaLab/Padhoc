#!/usr/bin/perl -w
use threads;
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

#print "min=".$minsstart."\n";
#$deshsp="";#pour_chaque_hsp_ecrire(q:[qstart]-[qend] s:[sstart]-[send]\t)
@jobs=();

mkdir($dirOut)unless(-e $dirOut);

#on recupere tous les nom des fichiers blasts
opendir(BL,$dirBlast);
@fb=readdir(BL);
#pour chaque fichier blast
foreach $f(@fb){
 unless(($f eq "." )||($f eq "..")){
   push(@jobs,threads->create(sub{
   	print "perl parser_pour_Inparanoid_job.pl $score_cutoff $dirBlast $dirOut $threadsAllowed $f";
	system("perl parser_pour_Inparanoid_job.pl $score_cutoff $dirBlast $dirOut $threadsAllowed $f");
   }));
   sleep(1) while(threads->list(threads::running)>=$threadsAllowed);
 }
}
closedir(BL);
$_->join for @jobs;

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



