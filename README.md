# Padhoc

**Dependencies:**
* Neo4j - https://neo4j.com/download-center/
* Vagrant virtual machine - https://app.vagrantup.com/CecilePereira/boxes/AdaptMetrecon​
* TEES - https://github.com/jbjorne/TEES/wiki
* tmChem - https://www.ncbi.nlm.nih.gov/research/bionlp/Tools/tmchem/

Note: If the installation of TEES with the configure.py program do not succeed to install svm_multiclass you may have to correct one line of the TEES file: Utils/DefaultSettings.py.
Try to change:
URL["SVM_MULTICLASS_LINUX"] = "http://download.joachims.org/svm_multiclass/current/svm_multiclass_linux.tar.gz"
for:
URL["SVM_MULTICLASS_LINUX"] = "http://download.joachims.org/svm_multiclass/current/svm_multiclass_linux32.tar.gz"

**How to use the program:**

1) Initialization 

1.1) Initialization of the database

* kill neo4j (top) (research the jobid by doing a restart)
* load the brenda database: neo4j-community-3.2.1$ sudo bin/neo4j-admin load --
* start neo4j again: in a terminal: neo4j-community-3.2.1$ ./bin/neo4j console from=/media/data/StressPath/Scripts_neo4j_2/extract_path_neo4j/neo4j_dbs/athalianaBrenda.dump -database=graph.db --force

1.2) initialization of the VM:

Use the boxe https://app.vagrantup.com/CecilePereira/boxes/AdaptMetrecon

2) Select one of the two fallowing steps:

2.1) for on specie:
The script allowing to extract metabolic pathways from pubmed publications:
* Script_v3.py
The list of the option is available with python2 Script_v3.py -h

Example of command line:
python2 Script_v3.py --brendaDB 0 --rep /home/cecile/data/StressPath/Scripts_neo4j_2/RESULTS/Extract_metabo_abioticstress_CitClem_CitSin_AraTha_PhyPat_16August2017//arabidopsisthaliana --maxm 1000 --repvagrant /home/cecile/data/StressPath/metrecon_cecile/ --part 4 --path "abiotic stress" --organism "Arabidopsis thaliana" --namerelation TMStress

An update of this scripts with other TEES (interactions proteins proteins...) will be available soon.
We are currently working on generate_graph.py Script_v3.py in order to provide this new functionnality.

2.2) For several species at the same time:

The script allowing to extract the metabolic pathway for several species at the time and provide the orthologs of those species (inparanoid):
* Extract_metabo_several_species.py (list of the options with -h)


