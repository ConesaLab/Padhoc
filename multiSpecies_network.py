#!/usr/bin/env python


from graph_database import graph_database
from Inparanoid import inparanoid
from uniprot_queries import uniprot_queries
from chebi_from_string import chebi_from_string
from activepool import ActivePool

import argparse, os, time, re
import multiprocessing, sys, glob



def run_text_mining(reptmp, s1, taxid, path, brendapass, maxm, email, repvagrant, s, pool, brendaDB, metrecon, tees, corporas, reptees, userdb, passdb, pazar, string, intact, omnipath, tf_data, repp):
    '''
    run text mining scripts.
    '''
    name = multiprocessing.current_process().name
    with s:
        pool.makeActive(name)
        sys.stderr.write("#################################")
        print "Now running: %s" %str(pool)
        stringspecie = "Specie : "+s1
        s2=re.sub(' ','',s1)
        sys.stderr.write(stringspecie + "\n")
        sys.stderr.write(time.strftime("%a, %d %b %Y %H:%M:%S",time.localtime()) + "\n")
        torun = "python2 run_padhoc.py -o " + reptmp + " -e " + email + ' --usergdb ' + userdb + ' --passgraphdb ' + passdb +\
        ' -b ' + str(brendaDB) + ' -s "' + s1 + '" -k "' + path + '" -p ' + brendapass + ' -m ' + str(maxm) + " -t " + taxid
        if repvagrant != '':
        	torun = torun + ' -v ' + repvagrant
        if metrecon == True:
        	torun = torun + ' -w'
        if tees == True:
        	torun = torun + ' -y'
        if len(corporas) > 0:
        	a = ''
        	for corp in corporas:
        		a = a + ' ' + corp
        	torun = torun + ' -c' + a
        if reptees != '':
        	torun = torun + ' -z ' + reptees
        if pazar == True:
        	torun = torun + ' --pazarDB'
        if string == True:
        	torun = torun + ' --stringDB'
        if tf_data == True:
        	torun = torun + ' --tfData'
        if omnipath == True:
        	torun = torun + ' --omnipath'
        if intact == True:
        	torun = torun + ' --intactDB'


        torun = torun + " > " + repp + "/text" + s2 + ".out"
        print torun
        sys.stderr.write(torun + "\n")
        os.system(torun)
        pool.makeInactive(name)


if __name__ == '__main__':


	organisms=""
	taxids = ''

	rep = ''
	email = ''
	repvagrant = ''

	path = ''
	brendapass = ''
	maxa=''
	repfa=""

	reptees=''
	brendaDB = 1

	usergraphdb = ''
	passgraphdb= ''


	parser=argparse.ArgumentParser(description="This script apply the extraction of metabolism reaction from text to several species iteratively and then extract the ortholog relationships from inparanoid.")
	parser.add_argument("-o","--rep",help='Path to the output repository',default=rep)
	parser.add_argument("-e","--email",help="User email (brenda user)", default=email)
	parser.add_argument('-v',"--repvagrant",help="Path to the repository where the vagrant machine is installed",default=repvagrant)
	parser.add_argument('-s',"--organisms",help="List of the organisms of interest (separated by a ',')",default=organisms)
	parser.add_argument('-t',"--taxids",help="List of the taxids linked to the organisms above (separated by a ',')",default=taxids)
	parser.add_argument('-k','--path',help="Keyword or pathway of interest (used for the pubmed request)",default=path)
	parser.add_argument('-p','--brendapass',help="Brenda password",default=brendapass)
	parser.add_argument('-m','--maxm',help='Maximum number of article extracted from the pubmed request for each specie',default=maxa,type=int)
	parser.add_argument('-d','--repb',help='Repository with blast results',default="")
	parser.add_argument('-f','--repfa',help='Repository with species fasta files, the fasta files should come from uniprot',default=repfa)
	parser.add_argument('-i','--repinpara',help='Output repository for inparanoid predictions',default="")
	parser.add_argument('-g','--cpu',help='Number of cpu for blast and inparanoid',default=1)
	parser.add_argument('-a','--allo',help='For the species of interest, keep all the orthologs in the inparanoid database (1) or only the one with the other species of interest (0).',default=1)
	parser.add_argument('--orthoWeb',help='1: Run the research of orthologs on the inparanoid website (done one ID by one ID, can take a long time), 0 squip this step. Default 0',default=0)
	parser.add_argument('-b','--brendaDB',help='Use the brenda database in order to populate the neo4j graph (may take several hours). Default No (0), 1 for yes.',default=brendaDB)
	parser.add_argument('-w', '--metrecon', help='Use metrecon program to obtain the relationships from text',action='store_true')
	parser.add_argument('-y', '--tees',help='Use tees program to obtain the relationships from text, choose corpora with -c, GE11 as default', action='store_true')
	parser.add_argument('-c', '--corpora', nargs='+', help='Choose tees corpora between GE09, GE11, GE13, BB11, BI11, CG13, CO11, DDI13, EPI11, GRN13, GRO13, ID11, PC13, REL11, REN11;'\
	    ' or providing the local path to a local corpora. GE11 as default', default=['GE11'])
	parser.add_argument('-z',"--reptees",help="Path to the repository where tees is installed",default=reptees)
	parser.add_argument("--usergdb",help="Neo4j graph database user ID",default=usergraphdb)
	parser.add_argument("--passgraphdb",help='Neo4j graph database password', default=passgraphdb)

	parser.add_argument("--pazarDB", action='store_true', help="Add the PAZAR DB to the Neo4j database, adding TF-gene regulation info.", default=False)
	parser.add_argument("--stringDB", action='store_true', help="Add the String DB to the Neo4j database.", default=False)
	parser.add_argument("--tfData", action='store_true', help="Add the TF-gene expression relationships to the Neo4j database.", default=False)
	parser.add_argument("--omnipath", action='store_true', help="Add the Omnipath DB to the Neo4j database.", default=False)
	parser.add_argument("--intactDB", action='store_true', help="Add the Intact DB to the Neo4j database.", default=False)
	args=parser.parse_args()


	if not os.path.exists(args.rep):
		os.mkdir(args.rep)

	if args.repfa=="":
		print "WARNING: the fasta repository must be define if at least one specie is not in the inparanoid online database"
	else:
		if not os.path.exists(args.repfa):
			os.mkdir(args.repfa)  

	if args.repb=="":
		args.repb=args.rep+'/blastp/'
		if not os.path.exists(args.repb):
			os.mkdir(args.repb)

	if args.repinpara=='':
		args.repinpara=args.rep+'/inportho/'
		if not os.path.exists(args.repinpara):
			os.mkdir(args.repinpara)

	######################################################
	# List of the arguments
	######################################################
	print "Script: Extract_metabo_several_species.py"
	print time.strftime("%a, %d %b %Y %H:%M:%S",time.localtime())
	print "Arguments:"
	dicargs = vars(args)
	for a in dicargs.keys():
		print a+' : '+str(dicargs[a])

	######################################################
	# 1) run script that extract relations from text and brenda for each of the species
	######################################################
	species = [n.strip() for n in args.organisms.split(',')]
	taxids = [n.strip() for n in args.taxids.split(',')]
	    
	jobs = []
	pool = ActivePool()
	sm = multiprocessing.Semaphore(args.cpu)

	for n in range(0,len(species)): #Foreach input specie
		t = taxids[n]
		s = species[n]
		print s, t
		# Write specie name and time and request to the stdr (in function run_text_mining)
		reptmp = args.rep + '/' + s.replace(' ', '')        #remove spaces, maj etc in the name of the repository
		p = multiprocessing.Process(target = run_text_mining, args = (reptmp, s, t, args.path, args.brendapass, args.maxm, args.email,\
		 args.repvagrant, sm,pool, args.brendaDB, args.metrecon, args.tees, args.corpora, args.reptees, args.usergdb, args.passgraphdb,\
		 args.pazarDB, args.stringDB, args.intactDB, args.omnipath, args.tfData, args.rep))
		jobs.append(p)

	for j in jobs:
		j.start()

	for j in jobs:
		j.join()
		sys.stderr.write('Now running: %s\n' % str(pool))

	######################################################
	# 2) extract the ortholog relationships present in the inparanoid database
	# Cecile: This part can take a lot of time (several days for a test with 4 species)
	######################################################

	uniprot = uniprot_queries(args.organisms.split(',')[0], args.taxids.split(',')[0])
	chebi = chebi_from_string()
	chebi.chebi_connect()

	gd = graph_database(chebi, uniprot, '', '', args.usergdb, args.passgraphdb)
	gd.connect()

	inparanoid = inparanoid(gd, args.repfa, args.repb, args.repinpara) # definition of the inparanoid class
	if args.orthoWeb: #if research of the orthologs on the website of inparanoid
		sys.stderr.write("For all the proteins in the neo4j database: research of the orthologs in inparanoid database\n")
		sys.stderr.write(time.strftime("%a, %d %b %Y %H:%M:%S",time.localtime())+"\n")
		# Now look if orthologs between those species
		for specie in species:
			uniprot_IDs = gd.obtain_uniIDs_from_specie(specie)
			for uniID in uniprot_IDs:
				if uniID != None:
					# how to test if a specie is in inparanoid database?
					inparanoid.research_orthologs(uniID, specie)



	######################################################
	# 3) RUN inparanoid
	######################################################

	sys.stderr.write("Research of the orthologs by applying inparanoid (run for missing)\n")
	sys.stderr.write(time.strftime("%a, %d %b %Y %H:%M:%S",time.localtime())+"\n")

	for specie in species:
		inparanoid.extract_fasta(specie)

	print "Run Blast on the Uniport IDs"
	inparanoid.doblast(species)
	print "Match blast results to extract homologies"
	inparanoid.doinparanoid()

	resfiles = glob.glob(args.repinpara + '/sqltable/sqltable*')
	for fname in resfiles:
		orthologs = inparanoid.parse_inparanoid(fname)
		inparanoid.add_blast_orthology_relationship(orthologs)