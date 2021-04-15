# Padhoc

**Dependencies:**
* Neo4j - https://neo4j.com/download/
* Vagrant virtual machine - https://app.vagrantup.com/scasani/boxes/adaptedMetrecon
* tmChem - https://www.ncbi.nlm.nih.gov/research/bionlp/Tools/tmchem/
* TEES - Provided with Padhoc in TEES-master


**Installing dependencies:**

*Libraries*

Run install_packages.sh file to install all dependencies required by Padhoc.

*Neo4j*

Follow Neo4j installing instructions, download desktop version and set user and password, these are needed for Padhoc to connect to the database.
Create an empty graph within Neo4j App. Graph version must be 3.5 or lower.


*tmChem*

Download and install tmChem following the package instructions, if CRF version from the package is outdated and prevents the installation from working download an updated version from this link: https://taku910.github.io/crfpp/#download
Once downloaded, place it to tmChem folder, change the nameof the folder to CRF. Install CRF following the software instructions and, once installed, repeat tmChem installation.


*TEES*

TEES is distributed by Padhoc and can be found at TEES-master folder. 
Modify line 33 of files Tools/TMCHEM.py and Tools/TMCHEM_noduplicate.py to the path of tmChem folcer location.
Please run configure.py script using the command: python configure.py
Follow the installation instructions, we strongly recommend to use TEES folder to locate .tees and .tees_local_settings.py
Every time Padhoc is used tees local settings file must be exported using command:
export TEES_SETTINGS=~path_to_tees/.tees_local_settings.py

Note: If the installation of TEES with the configure.py program does not succeed to install svm_multiclass you may have to correct one line of the TEES file: Utils/DefaultSettings.py.
Try to change:
URL["SVM_MULTICLASS_LINUX"] = "http://download.joachims.org/svm_multiclass/current/svm_multiclass_linux.tar.gz"
for:
URL["SVM_MULTICLASS_LINUX"] = "http://download.joachims.org/svm_multiclass/current/svm_multiclass_linux32.tar.gz"

*Brenda*

Create BRENDA user and password at https://www.brenda-enzymes.org/register.php, these are required by Padhoc to SOAP access the BRENDA database.


**How to use the program:**

1) Initialize the database

* Enter into Neo4j program and select the database to initialize, click start button
* In case a previously create database aims to be started, click the terminal within neo4j graph, and import the database using neo4j-admin load command. Then change the running database at settings file and reinitialize the database.


2) Run Padhoc either in one species or in multiple species.

2.1) Run Padhoc for one species:
The main script used to run padhoc is: run_padhoc.py
The options can be accessed using command: python run_padhoc.py -h 
Compulsory arguments are: repository, email, password, user graph db, user db password, organism and taxa ID

 - Databases: The available databases include BRENDA (-b), pazar, string, tfData, omnipath and intact database.
 - It is possible to run Padhoc with different input sources:
    - Text search: -k to select the pathway to search. -m to select the maximum number of papers to download.
    - List of PMIDs: --pmid to input a file with pmids, there should be one pmid per line in the file.
    - Input text: -i to select a folder with the available text to start the text mining engine.

Command line example: 
python run_padhoc.py -o output_dir -k path_to_search -m 100 -y -c GE09 GE11 -z path_to_teesFolder --usergdb neo4j_user_db --passgraphdb neo4j_password -b 1 -s 'Homo sapiens' -t taxa_id -e brenda_email -p brenda_password

2.2) Several species:
This option is run using the multiSpecies_newtork.py script. The options to be used are displayed using -h, but they are basically the same as in the case for one species. 
This script runs several times the run_padhoc script, once for each input species, and at the end of the script searches in InParanoid for possible homologies between the proteins of the different species. It also includes paralogies within one species.


3) Graph Compression

All information from databases and text has now been introduced to Neo4j, the database contains thus redundant information that does not necesarily add knowledge to the graph. Padhoc offers the user the possibility to compress the graph database, which clusters the nodes that have a similar name. This step is strongly recommended and is run from padhoc's folder compress_graph as follows:

  - One species: python -m path.to.padhoc.compress_graph.compress_graph
  - More than one species: python -m path.to.padhoc.compress_graph.multiSpecies_compress_graph
  
This step can take long if the subgraph is big, since padhoc compares similarities between all entities in the database and clusters them.


4) Network Visualization

The network is displayed using Neo4j visualization capabilities. Neo4j is accessed using a non-SQL programming language called Cypher, this language allows the user flexibility to retrieve information from the database by using node patterns or properties from nodes and edges, we recommend using the following cypher query to retrieve all compressed nodes, which include nodes that have been found in text that have a database annotation:

(n)-[r:Compressed_relationship]-(y) RETURN n,y

It is possible to retrieve alternative information by using other cypher queries (tutorial for Cypher language: https://neo4j.com/developer/cypher-basics-i/).


**Toy dataset**

We have included a toy dataset to show the principles of Padhoc and for the user to ensure that the installation works properly. The toy example data is included with Padhoc in the folder toy_dataset.

1) Import Homo sapiens database in Neo4j
Create a new grap (Add graph > Create Local Graph), recommended version: 3.5.14. Select Manage and Open Terminal.
Once the terminal is started, write the following command, substituting the "path_to_padhoc" with your local path:

 bin/neo4j-admin load --from="path_to_padhoc"toy_dataset/homoSapiens_db.dump --database=hsapiens_example.db
 
Close the terminal. Change the database that Neo4j is using: Go to Settings > change line 11 where it says "#dbms.active_database=graph.db" to dbms.active_database=hsapiens_example.db (deleting the #). Apply changes.

2) Start Neo4j database
Start the database by clicking the "play" symbol. Once it is running you can navigate the database to visualize the structure.

3) Run Padhoc
If you haven't, export the TEES settings: export TEES_SETTINGS=~path_to_tees/.tees_local_settings.py
Go to your padhoc folder in the terminal and run the following command:

python run_padhoc.py -i toy_dataset/text6U1CR1X4/ -o padhoc_example -y -c GE09 GE11 --usergdb "neo4j_user" -b 0 --passgraphdb "neo4j_password" -z TEES-master -s "Homo sapiens" -t 9606

This process can take around half an hour.

4) Visualization
In Neo4j interface, click on "Open Browser". When the database opens and loads, the user can interact with it through Cypher queries. The database should contain a total of x nodes and x relationships. When writing the following Cypher command:

MATCH (n)-[r]-(y) WHERE n.id = 'O15516' AND y.id = 'CHEBI:15346' RETURN n.uniprotEntryName, r.ECs ,y.compoundName

You should obtain the following result:

"CLOCK_HUMAN"	["2.3.1.48"]	"coenzyme A"

Congratulations,you have succeeded creating your first subgraph :)
