#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  SCRIPTS_INPARANOID.pl
#
#        USAGE:  ./SCRIPTS_INPARANOID.pl  
#
#  DESCRIPTION:  lancer l'ensemble des scripts pour inparanoid en une seule fois
#
#      OPTIONS:  dossier_res_blastp nb_especes(pr_info)  dos_res nb_cpu
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Cecile Pereira 
#      COMPANY:  
#      VERSION:  1.0
#      CREATED:  26/06/2013 15:24:47
#     REVISION:  ---
#===============================================================================

#use strict;
use warnings;
use Cwd;
$home=cwd."/Inparanoid/";
chdir($home);
$param{'-score_inparanoid'}=40;
#$DOSBLAST="/archives/pereira/ORTHOMCL_PRD_201104/blastp/";
#/home/cpereira/Protein_Reference_Dataset_tmp/blastp
$DOSBLAST=$ARGV[0];
#$DATA="/archives/pereira/BRH_PRD_201104/2011_04_reference_proteomes_fasta/";
#$DATA="/home/cecile/PRD_2011_FUNGI/Fungi_PRD2011/";
$DATA=$ARGV[1];
#$nbf=$ARGV[2];
$nbblasts=$home."/INPARA_blastp/";
$res=$ARGV[2];
$cpu=$ARGV[3];
print "cpu $cpu\n";
mkdir($res);
print "method Inparanoid\n";
system("date");
print "parser blast\n";
system("perl $home/parser_pour_Inparanoid.pl $param{'-score_inparanoid'} $DOSBLAST/ $nbblasts $cpu");
system("date");
print "2_inparanoid\n";
system("perl $home/2_inparanoid.pl $home/ $DATA/ $nbblasts $res/ $cpu");
system("date");


