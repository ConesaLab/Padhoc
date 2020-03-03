#!/usr/bin/env python

'''
'''


from brenda_annotation import brenda_annotation
from parse_files import parse_file
from graph_database import graph_database
from entity_recognition import entity_recognition
from parse_abrev import parse_abrev
from extract_articles import extract_articles
from run_tees import run_tees
from parse_tmchem_file import parse_tmchemfile
from run_vagrant_metrecon import run_vagrant_metrecon
from uniprot_matcher import uniprot_matcher
from chebi_from_string import chebi_from_string
from uniprot_queries import uniprot_queries
from compare_string_to_db import compare_string_to_db

from pazar_database import pazar_database
from intact_database import intact_database
from string_database import string_database
from TF_data import TF_data
from omnipath_db import omnipath

import argparse
import string
import sys, os, time, glob, re
import random
from time import strftime



def clean(word,specie=''):
    '''
    Clean words to allow them comparison.
    See if I need to add all the cases that I covered with perl before or if the score in the chebi comparison is enouth
    '''
    #lowercases
    cleanw=word.lower() #word in small letters
    specie=specie.lower()
    sp=specie.split(' ')
    specie=specie.strip()
    #no whitecharacters
    cleanw=cleanw.strip()
    #remove parenthesis
    cleanw=re.sub('^(r)-','',cleanw)
    cleanw=re.sub('^(s)-','',cleanw)
    cleanw=re.sub('\)n$','',cleanw) #remove the indications of several time this compound
    cleanw=re.sub('\)m$','',cleanw)
    cleanw=re.sub('^\(\+?\-?\d{0,5}\)','',cleanw)
    cleanw=re.sub('\(','',cleanw)
    cleanw=re.sub('\)','',cleanw)
    cleanw=re.sub('\[','',cleanw)
    cleanw=re.sub('\]','',cleanw)
    cleanw=re.sub('\{','',cleanw)
    cleanw=re.sub('\}','',cleanw)
    cleanw=re.sub('^l-','',cleanw)
    cleanw=re.sub('^d-','',cleanw)
    cleanw=re.sub('^r-','',cleanw)
    cleanw=re.sub('^s-','',cleanw)
    cleanw=re.sub('^ec:','',cleanw)
    cleanw=re.sub('^ec ','',cleanw)
    cleanw=re.sub('^alpha-','',cleanw)
    cleanw=re.sub('^beta-','',cleanw)
    cleanw=re.sub(' genes{0,1}','',cleanw)
    cleanw=re.sub(' proteins{0,1}','',cleanw)
    #cleanw=re.sub('\)','',cleanw)
    #remove words never used in the text
    cleanw=re.sub(' atom','',cleanw)
    cleanw=re.sub('s$','',cleanw) #remove plural
    cleanw=cleanw.replace(' ','') #test 3 july
    cleanw=cleanw.replace('"','')
    cleanw=cleanw.replace('-','')
    cleanw=cleanw.replace("'",'')
    for s in sp:
        cleanw=re.sub(s,'',cleanw)#remove the specie name from the entity name (ex for Arabidopsis thaliana: remove arabidopsis and then remove thaliana)
    if cleanw=='':
        cleanw='empty_string' #to avoid error when we remove all the carracters of the script with this function
    return cleanw


def id_generator(size=8, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))



if __name__ == '__main__':
    

    ############################################
    # Default values:
    email = ''
    listpmidfile = ""
    #repvagrant = '/home/salva/metrecon/'
    repvagrant = ''
    #reptees = '/home/salva/TEES-master'
    reptees = ''
    path = None

    brendapass = 'salvacasani91'
    maxa = 70 # max number of publications (default value)
    brendaDB = 0
	
    rep = ''
    organism = ''
    organism_tax = '' # to avoid problem with several strains (used only with brenda)
    usergraphdb=''
    passgraphdb=''
    
    #################################################
    # PARAMETERS :
    #################################################
    parser=argparse.ArgumentParser(description="For one specie, extract relations from articles (pubmed request) and the reactions from Brenda. The result is stored in a neo4j databaseSalva Casani: \nsalvacasani@gmail.com\nCecile Pereira: pereiracecile@live.fr")
    parser.add_argument("-o","--rep",help='Path to the output repository',default=rep)
    parser.add_argument("-e","--email",help="User email (brenda user)", default=email)
    parser.add_argument('-v',"--repvagrant",help="Path to the repository where the vagrant machine is installed",default=repvagrant)
    parser.add_argument('-s',"--organism",help="Organism of interest",default=organism)
    parser.add_argument('-k','--path',help="Keyword or pathway of interest (used for the pubmed request)",default="")
    parser.add_argument('-p','--brendapass',help="Brenda password",default=brendapass)
    parser.add_argument('-m','--maxm',help='Maximum number of article extracted from the pubmed request',default=500,type=int)
    parser.add_argument('-t','--organism_tax',help='Organism taxID',default=organism_tax)
    parser.add_argument('-b','--brendaDB',help='Use the brenda database in order to populate the neo4j graph (may take several hours). Default No (0), 1 for yes.',default=brendaDB)
    parser.add_argument('-w', '--metrecon', help='Use metrecon program to obtain the relationships from text',action='store_true')
    parser.add_argument('-y', '--tees',help='Use tees program to obtain the relationships from text, choose corpora with -c, GE11 as default', action='store_true')
    parser.add_argument('-c', '--corpora', nargs='+', help='Choose tees corpora between GE09, GE11, GE13, BB11, BI11, CG13, CO11, DDI13, EPI11, GRN13, GRO13, ID11, PC13, REL11, REN11;'\
        ' or providing the local path to a local corpora. GE11 as default', default=['GE11'])
    parser.add_argument('-z',"--reptees",help="Path to the repository where tees is installed",default=reptees)
    parser.add_argument('-i',"--inputtext",help="Folder with a set of articles of interest. Optional to skip article downloading step",default=None)
    parser.add_argument("--usergdb",help="Neo4j graph database user ID",default=usergraphdb)
    parser.add_argument("--passgraphdb",help='Neo4j graph database password', default=passgraphdb)
    parser.add_argument("--pmid",help="File with a series of pmids to download", default=listpmidfile)

    #Databases args
    parser.add_argument("--pazarDB", action='store_true', help="Add the PAZAR DB to the Neo4j database, adding TF-gene regulation info.", default=False)
    parser.add_argument("--stringDB", action='store_true', help="Add the String DB to the Neo4j database.", default=False)
    parser.add_argument("--tfData", action='store_true', help="Add the TF-gene expression relationships to the Neo4j database.", default=False)
    parser.add_argument("--omnipath", action='store_true', help="Add the Omnipath DB to the Neo4j database.", default=False)
    parser.add_argument("--intactDB", action='store_true', help="Add the Intact DB to the Neo4j database.", default=False)

    args=parser.parse_args()
    
    #####################################################
    # PARAMETERS FOR THE TESTS :  
    part = 0 # subpart to execute, 0: run all the program

    ####################################
    # 0) PRE-PROCESSES
    ####################################
    # Create the output repository
    if not os.path.exists(args.rep):
        os.mkdir(args.rep)
    else:
        print "The output repository exist"
        
    randomstring = id_generator() # random string used in order to be able to run metrecon on several repositories at the same time (different names)
    #randomstring = 'G7VY1A6Z' #to remove, only present for the test
    
    #####################################################
    # Write the parameters in a output parameter file
    f = open(args.rep+'/script_parameter_'+randomstring+'.out','w')
    f.write("Script: generate_graph.py\n")
    f.write(strftime("%a, %d %b %Y %H:%M:%S",time.localtime()))
    f.write("\nArguments:\n")
    dicargs=vars(args)
    for a in dicargs.keys():
        f.write(a+' : '+str(dicargs[a])+"\n")
    f.close()

    # Output repositories and requests
    if (args.pmid == "") and (args.inputtext == None):
        if args.organism != None:
            req = args.path + " AND " + args.organism + "[MeSH Terms]" # pubmed request
        else:
            req = args.path


    brendauser = args.email # email for brenda
    namea1a2rep='a1a2_text'+randomstring # output repositories file names for the texmining by TEES and metrecon

    if args.inputtext == None: #no text provided by the user (need to request pubmed and extract the articles)
        textrep= args.rep + '/text'+randomstring 
        nametextrep = '/text'+randomstring
    else: # if the input text files are provided 
        torun = 'cp -r %(textfile)s %(out_rep)s' % {'textfile':args.inputtext, 'out_rep':args.rep}
        os.system(torun)
        if args.inputtext.split('/')[-1] != '':
            text = args.inputtext.split('/')[-1]
        else:
            text = args.inputtext.split('/')[-2]
        textrep = args.rep + '/' + text


    A1A2REP=args.rep+'/'+namea1a2rep+'/' #output repository path

    start = time.time() # start date
    
    if args.pmid == "":
        query=args.organism+"%"+args.path+"%"+str(args.maxm)
    else:
        query = args.pmid
    
    ###############################################################################################
    ## 1) Extract abstracts and complete texts from medline and PMC and recover the abreviations
    ###############################################################################################
    listpmid=[]
    if part == 0:
        if args.inputtext == None: # Unless a repository with publications (text) is already provided
            if(os.path.isfile(args.pmid)):
                print "Extraction of the articles:"
                print "The list of PMID is given as input"
                article=extract_articles('',listpmidfile+'_read',args.rep,args.email,args.maxm,nametextrep)
                listpmid=article.read_pmidlist(args.pmid)
                
            else:
                if args.path != '':
                    firstmessage="1) Extraction of the "+str(args.maxm)+" articles\n"
                    sys.stderr.write(firstmessage)
                    sys.stderr.write(strftime("%a, %d %b %Y %H:%M:%S",time.localtime())+"\n")
                    articles=extract_articles(req,listpmidfile,args.rep,args.email,args.maxm,nametextrep)
                    listpmid=articles.pubmed_request()


    ###############################################################################################
    # 2) Extraction of the relations from the BRENDA database
    ###############################################################################################
    #class needed from salva script: brenda, uniprot (cecile: I only have uniprot_matcher), graph_database (all functions?), chebi (class chebi_from_string?)
    brenda = brenda_annotation(args.email,args.brendapass) # create brenda_annotation object
    brenda.access_protocol() # access to brenda
    chebi = chebi_from_string()
    chebi.chebi_connect()
    uniprot = uniprot_queries(args.organism, args.organism_tax)
    gd=graph_database(chebi,uniprot,brenda,args.organism,args.usergdb,args.passgraphdb) #
    gd.connect()

    #for each long name the corresponding chebiID
    if int(args.brendaDB) == 1:
        # test if brenda already run for the specie only if brendaDB==1
        specie = gd.check_specie() #if specie absent from the database #Cecile: look if OK with the part orthologs when I download ortho from inparanoid website. 
        # OK if I don't download orthologs of species that are not requested
        if specie == False:
            organisms = open('speclist.txt').read()  # this file come from uniprot database (provided in the git repository)
            ec_converter = open('enzyme.dat').read() # (file provided in the git repository) This file may change in several years, TODO see if I can download it from the script. 
            #why do we need the expasy file? Can we do by an other way to avoid deprecacy link to update of this file? 
            #can we download this file automatically?
            
            uniprot_match = uniprot_matcher(organisms)
            uniprot_match.EC2Uniprot(ec_converter, args.organism_tax)

            gd.create_database(uniprot_match.ec2Uniprot_dict)

    ###############################################################################################
    ##  Add the databases that were requested. Multiple possibilities for TF-gene 
    ##  regulation and Protein interaction
    ###############################################################################################


    if int(args.tfData) == 1:
        print 'Adding TF - gene regulation data to the neo4j database\n\n'
        if args.organism == 'Homo sapiens':
            tf = TF_data(gd, uniprot)
            tf.process_files()
        else:
            print 'This information is only available for Homo sapiens\n\n'


    if int(args.pazarDB) == 1:
        print 'Adding TF - gene regulation data fro PAZAR DB to neo4j graph db\n\n'
        pazar = pazar_database(args.organism, gd, uniprot)
        pazar.process_files()
        pazar.fill_database()


    #For the string database, the file for the desired organisms has to be included 
    #within the data folder and added to the following function
    if int(args.stringDB) == 1:
        print 'Adding protein interaction data from String database\n\n'
        #stringFile = 'data/string_data/511145.protein.links.detailed.v10.5.txt'
        string = string_database(gd, args.organism, uniprot)
        if args.organism == 'Homo sapiens':
            string.process_dbFile('data/string_data/9606.protein.links.detailed.v10.5.txt', 500)
        if args.organism == 'Escherichia coli':
            string.process_dbFile('data/string_data/511145.protein.links.detailed.v10.5.txt', 500)
        else:
            print 'Organism is not included within String options, please add it to the data folder as specified in Padhoc guidelines'


    if int(args.intactDB) == 1:
        print 'Adding protein interaction data from IntAct database\n\n'
        intact = intact_database(gd, args.organism_tax, uniprot)
        intact.process_intact()


    if int(args.omnipath) == 1:
        print 'Adding signaling information from Omnipath\n\n'
        if args.organism == 'Homo sapiens':
            omnipath = omnipath(gd, uniprot)
            omnipath.process_omnipath()
        else:
            print 'Information only available for Homo sapiens\n\n'

    # Delete
    #Test this dictionary
    dico_synonyms_prot=gd.extractsynonyms()#dico synonyms proteins
    ###############################################################################################
    ## 2) Run metrecon (adapted version of TEES)
    ## - how to do that without the virtual machine ? 
    ##          * Use virtual machine but put my new scripts into a box: CECILE: I am preparing the box (08/09/2017)
    ## - types of links ? (oriented?) #TEES yes, metrecon ?
    ## - metabolites: BANNER and TMCHEM
    ## - genes: BANNER
    ###############################################################################################
    sys.stderr.write("2) Extraction of the metabolic reactions from the text\n")
    sys.stderr.write(strftime("%a, %d %b %Y %H:%M:%S",time.localtime())+"\n")
    
    if (part == 0 or part == 1) and len(os.listdir(textrep)) > 0:

        if args.metrecon == True: # SEEMS OK
            runmetrecon = run_vagrant_metrecon(textrep, args.rep + '/metrecon/',args.repvagrant)
            runmetrecon.split_files(randomstring)
        # TEST TEES PART
        if args.tees == True: # CECILE: question: Do you run banner and tmchem? if yes just once or for each TEES? 
			      # CECILE: Do you use metabolite or only proteins?
			      # CECILE: Can we use banner / tmchem run of metrecon in order to save time?
            tees = run_tees(args.rep + '/tees/')
            for corpora in args.corpora:
                tees.split_input(randomstring, args.reptees, corpora, textrep) 


    ###############################################################################################
    ## Sentences
    ###############################################################################################
    # read the file placeholder-preprocessed.xml.gz-ner.xml.gz,
    # it contain sentences and position of the compounds/proteins relative to the sentence start
    sys.stderr.write("3) Correspondance sentences / proteins\n")
    sys.stderr.write(strftime("%a, %d %b %Y %H:%M:%S",time.localtime())+"\n")
    
    
    #TODO: Should we include all this in an external function? Maybe within the parse_files class?

    dicoartTsentenceMet={}; dicoartTsentenceTees={}
    dicoPMID_idDoc_met={}; dicoPMID_idDoc_tees={}
    
    if part == 0 or part == 4 :
        if os.path.exists(args.rep + '/metrecon/'):
            repa1a2tmp=glob.glob(args.rep + '/metrecon/a1a2*')
            for r in repa1a2tmp:
                #nerfile = A1A2REP+'/placeholder-preprocessed.xml.gz-ner.xml.gz'
                nerfile=r+'/placeholder-preprocessed.xml.gz-ner.xml.gz' #modif 10 july
                #print "#####################\n nerfile:"
                #print nerfile 
                pner = parse_file(nerfile,dicoartTsentenceMet,dicoPMID_idDoc_met,r)
                dicoartTsentenceMet, dicoPMID_idDoc_met = pner.docTsent() # dictionary with keys [PMID][T0]=sentence and dico PMID_idDOC in TEES
                #print len(dicoartTsentence.keys())
                #print len(dicoPMID_idDoc.keys())
        for corpora in args.corpora:
            if os.path.exists(args.rep + '/tees/' + corpora):
                repa1a2tmp = glob.glob(args.rep + '/tees/' + corpora + '/a1a2*')
                for r in repa1a2tmp:
                    nerfile=r+'/-preprocessed.xml.gz-ner.xml.gz'
                    pner = parse_file(nerfile,dicoartTsentenceTees,dicoPMID_idDoc_tees,r)
                    dicoartTsentenceTees, dicoPMID_idDoc_tees = pner.docTsent()
        

    if part == 0 or part == 4 :
        
        # Initialize the neo4j database
        #gd=graph_database_cecile(args.organism)
        #gd.connect()
        
        ################################################################################################
        # Normalize the entities in the A1 file:
        # 1: READ the file sentences.pubtator.tmChem to obtain the normalization given by tmchem:
        # 2: CONVERT abreviations into long version of the text.
        # 3: Compounds: request to chebi / Proteins: request to uniprot
        
        # Extract normalization from tmChem
        sys.stderr.write("4) Read tmchem result\n")
        if os.path.exists(args.rep + '/metrecon/'):
            #Salva: Why do we separate by articles in this dictionary? Why not the same dictionary for all the files?
            repa1a2tmp_met=glob.glob(args.rep + '/metrecon/a1a2*')
            dictmchem_met={}
            for r in repa1a2tmp_met:
                #print r
                sys.stderr.write(strftime("%a, %d %b %Y %H:%M:%S",time.localtime())+"\n")
                tmchemfile=r+'/sentences.pubtator.tmChem'
                tm=parse_tmchemfile(tmchemfile,dictmchem_met,r)
                dictmchem_met=tm.parse() #in dictmchem: ID doc, name of the compound, ID CHEBI

        
        for corpora in args.corpora:
            if os.path.exists(args.rep + '/tees/' + corpora):
            #Salva: Why do we separate by articles in this dictionary? Why not the same dictionary for all the files?
                repa1a2tmp_tees=glob.glob(args.rep + '/tees/' + corpora + '/a1a2*')
                dictmchem_tees={}
                for r in repa1a2tmp_tees:
                #print r
                    sys.stderr.write(strftime("%a, %d %b %Y %H:%M:%S",time.localtime())+"\n")
                    tmchemfile=r+'/sentences.pubtator.tmChem'
                    tm=parse_tmchemfile(tmchemfile,dictmchem_tees,r)
                    dictmchem_tees=tm.parse()
        print "######################"
        
        sys.stderr.write("5) Access the abreviations\n")
        sys.stderr.write(strftime("%a, %d %b %Y %H:%M:%S",time.localtime())+"\n")
        # Correspondance abrev-long
        # Salva: What is this exactly and where does it come from? It cannot be used if the text is provided instead of downloaded
        abrevrep=args.rep+'/abrev'
        abrev=parse_abrev(abrevrep)
        dicoabrev=abrev.parse() # in dicoabrev: doc PMID (warning, diff from IDdoc), abrev, long version

        ### Metrecon a1 a2 files read, and graph creation
        
        compareDB = compare_string_to_db(gd)
        compareDB.neo4j_to_dictio()

        if os.path.exists(args.rep + '/metrecon/'):

            Entity_recognition = entity_recognition(chebi, uniprot, dicoabrev, gd, dicoartTsentenceMet, dico_synonyms_prot, args.organism)
            #entity_recognition.fill_dictionaries_from_neo4j()

            for r in repa1a2tmp_met:
                listpmid=glob.glob(r + '/*a2')
                for pmid_path in listpmid:
                    pmid = pmid_path.split('/')[-1][:-3]
                    print pmid
                    if os.stat(pmid_path).st_size != 0:
                        Entity_recognition.process_a1_changed(r, pmid, dictmchem_met, dicoPMID_idDoc_met, query, compareDB)
                        Entity_recognition.process_a2(r, pmid, 'metrecon', query)


        Entity_recognition = entity_recognition(chebi, uniprot, dicoabrev, gd, dicoartTsentenceTees, dico_synonyms_prot, args.organism)
        #entity_recognition.fill_dictionaries_from_neo4j()

        
        for corpora in args.corpora:
            print '\n' + corpora + '\n'

            corpora_name = corpora
            if len(corpora.split('/')) > 1:
                corpora_name = corpora.split('/')[-2]

            repa1a2tmp_tees=glob.glob(args.rep + '/tees/' + corpora_name + '/a1a2*')
            if os.path.exists(args.rep + '/tees/' + corpora_name):

                for r in repa1a2tmp_tees:
                    listpmid=glob.glob(r + '/*a2')
                    for pmid_path in listpmid:
                        pmid = pmid_path.split('/')[-1][:-3]
                        print pmid
                        if os.stat(pmid_path).st_size != 0:
                            Entity_recognition.process_a1_changed(r, pmid, dictmchem_tees, dicoPMID_idDoc_tees, query, compareDB)
                            Entity_recognition.process_a2(r, pmid, corpora_name, query)
        
        sys.stderr.write(strftime("%a, %d %b %Y %H:%M:%S",time.localtime())+"\n")
        print strftime("%a, %d %b %Y %H:%M:%S",time.localtime())
        print "---END---"
        sys.stderr.write("---END---\n")

    ################### END ##########################################################################
    
   
