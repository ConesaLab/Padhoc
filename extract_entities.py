#!/usr/bin/env python

'''
Process the entity recognition files, obtain the entities
detected bu the Named Entity recognition algorithm (BANNER)
'''

import os
import re
from sys import argv
import glob
import string
import hashlib
from SOAPpy import WSDL, SOAPProxy
from socket import error as socket_error
from xml.sax import SAXParseException
import mysql.connector
from mysql.connector import errorcode
import time
import itertools
from neo4j.v1 import GraphDatabase, basic_auth
import xml.etree.ElementTree as ET
import zeep
from zeep import Client
  

class entity_recognition():
	'''
	Class to iterate through the a1 files, obtain the 
	entities recognized by the NER algorithm and categorise
	them
	'''


	def __init__(self, path=os.getcwd()):
		
		self.path = path
		self.files_list = []
		self.metabolites = set()
		self.proteins = set()
		self.a1_dictio = {}
		self.index_dictio = {}


	def obtain_a1_files(self):
		'''
		Identify all the a1 files present in the 
		given path
		'''
		for file in glob.glob(self.path + '/*a1'):
			self.files_list.append(file)
		return None


	def create_entity(self):
		'''
		Use a1 files info to create the 
		'''
		self.obtain_a1_files()
		for art in self.files_list:
			file = open(art).read()
			file_name = art.split('/')[-1]
			self.a1_dictio[file_name] = {}
			self.index_dictio[file_name[:-3]] = {}
			for line in file.split('\n'):
				if line != '':
					line = line.split('\t')
					tag = line[0]
					indexes = [line[1].split()[1], line[1].split()[2]]
					entity = line[-1]
					if line[1].split()[0] == 'Protein':
						self.proteins.add(entity)
					elif line[1].split()[0] == 'Metabolite':
						self.metabolites.add(entity)
					else:
						pass
					self.a1_dictio[file_name][tag] = entity
					self.index_dictio[file_name[:-3]][tag] = indexes
		return self.proteins, self.metabolites, self.a1_dictio


class entity_normalization():
	'''
	'''

	def __init__(self, entity):

		self.entity = entity


	def func():
		'''
		'''
		return None


class uniprot_matcher():
	'''
	Match a protein to the UniProt database to check if it
	is a "real" protein, and to obtain the UniProt identifier
	'''

	def __init__(self, organisms_file):

		self.organisms_file = organisms_file
		self.ec2Uniprot_dict = {}

	def uniprot_accession(self, inp, out, query):
		'''
		Protocol to access the uniprot database.
		'''
		url = 'http://www.uniprot.org/uploadlists/'
		params = {
			'from': inp,
			'to': out,
			'format':'tab',
			'query': query
		}
		data = urllib.urlencode(params)
		request = urllib2.Request(url, data)
		response = urllib2.urlopen(request)
		page = response.read(200000)
		return page


	def uniprotOrganism(self, taxa):
		'''
		Convert from organism to UniProt organism using the
		conversion file. For this function we need the 
		organism taxonomic number
		'''
		uni_name = None
		for line in self.organisms_file.split('\n'):
			line = line.split()
			try:
				if line[2][:-1]== taxa:
					uni_name = line[0]
			except IndexError:
				pass
		return uni_name


	def fill_ec_dictio(self, ec, uniprotID, uniprotGene):
		'''
		Fill the dictionary of ECs to uniprot identifiers
		'''
		if ec in self.ec2Uniprot_dict.keys():
			self.ec2Uniprot_dict[ec].append((uniprotID, uniprotGene))
		else:
			self.ec2Uniprot_dict[ec] = [(uniprotID, uniprotGene)]
		return None


	def EC2Uniprot(self, converter, taxa):
		'''
		Use the expasy file to convert from EC to UniProt.
		Filter by organism using the uniprot organisms file
		'''
		organism = self.uniprotOrganism(taxa)
		if organism != None:
			for line in converter.split('//'):
				if line.strip().startswith('ID'):
					line = line.split('\n')
					ec_id = line[1].split()[1]
					for l in line:
						if l.startswith('DR'):
							l = l.replace('DR', '').replace(' ', '').split(';')
							for uniprot in l[:-1]:
								uniprotID = uniprot.split(',')[0]
								uniprotGene = uniprot.split(',')[1]
								if uniprotGene.split('_')[-1] == organism:
									self.fill_ec_dictio(ec_id, uniprotID, uniprotGene)
		else:
			print "Organism not present in the uniprot database, constructing EC-oriented graph"
		return None


class sentence_retriever():

	def __init__(self, path=os.getcwd()):

		self.file = path + '-preprocessed.xml.gz-sentences.xml'
		self.sentence_dict = {}


	def get_sentence(self):
		'''
		Obtain the sentences for each of the articles by parsing
		the sentences xml file. This info will be stored in the 
		dictionary, together with the indexes of the sentence.
		'''
		sentence_file = open(self.file).read()
		root = ET.fromstring(sentence_file)
		tees2id = {}
		#Create a local dictio to pass from TEES id to PMID
		for child in root:
			PMID = child.attrib['origId']
			teesId = child.attrib['id']
			tees2id[teesId] = PMID
			self.sentence_dict[PMID] = {}
		#Fill the index dictio with the index and sentence information
		for sentence in root.iter('sentence'):
			sentence_text = sentence.attrib['text']
			sentence_index = sentence.attrib['charOffset']
			teesId = sentence.attrib['id'].split('.')[0] + '.' + sentence.attrib['id'].split('.')[1]
			PMID = tees2id[teesId]
			self.sentence_dict[PMID][sentence_index] = sentence_text
		return None


class TM_relationship():
	'''
	Recover the relationships obtained by running TEES in the 
	articles.

	Here, we end with a dictionary containing all the relationships 
	(by name), and an extra category specifying which relationships
	are included within other rels, and which contain other rels

	The interesting final result is encoded in the rels_dictio
	'''


	def __init__(self, path=os.getcwd()):

		self.path = path
		self.files_list = []
		self.a2_dictio = {}
		self.a2_rels_dict = {}
		self.a2_tags_dict = {}


	def obtain_a2_files(self):
		'''
		Identify all the a2 files present in the given path
		'''
		for file in glob.glob(self.path + '/*a2'):
			self.files_list.append(file)
		return None


	def fill_rel_dict(self, a1_dict, relationship, file_name):
		'''
		Fill the relationships dictionary using the relationship
		described by the a2 file
		'''
		relationship = relationship.split('\t')
		rel_tag = relationship[0]
		rel = relationship[1]
		self.a2_rels_dict[file_name[:-3]][rel_tag] = {'Theme':[], 'Cause':[], 'Rel':[]}
		for elem in rel.split():
			elem = elem.split(':')
			if elem[0].startswith('Theme'):
				if elem[1].startswith('T'):
					element = a1_dict[file_name[:-1] + '1'][elem[1]]
				elif elem[1].startswith('E') or elem[1].startswith('R'):
					element = elem[1]
				self.a2_rels_dict[file_name[:-3]][rel_tag]['Theme'].append(element)
			elif elem[0].startswith('Cause'):
				if elem[1].startswith('T'):
					entity = a1_dict[file_name[:-1] + '1'][elem[1]]
				elif elem[1].startswith('E') or elem[1].startswith('R'):
					entity = elem[1]
				self.a2_rels_dict[file_name[:-3]][rel_tag]['Cause'].append(entity)
			else:
				relationship = self.a2_dictio[file_name][elem[1]]
				self.a2_rels_dict[file_name[:-3]][rel_tag]['Rel'].append(relationship)
		return None


	def fill_tag_dict(self, relationship, file_name):
		'''
		Use the tags to fill the dictionary
		'''
		relationship = relationship.split('\t')
		tag = relationship[0]
		self.a2_tags_dict[file_name[:-3]][tag] = {'Theme':[], 'Cause':[], 'Rel':[]}
		for elem in relationship[1].split():
			elem = elem.split(':')
			if elem[0].startswith('Theme'):
				element = elem[1]
				self.a2_tags_dict[file_name[:-3]][tag]['Theme'].append(element)
			elif elem[0].startswith('Cause'):
				entity = elem[1]
				self.a2_tags_dict[file_name[:-3]][tag]['Cause'].append(entity)
			else:
				relationship = elem[1]
				self.a2_tags_dict[file_name[:-3]][tag]['Rel'].append(relationship)
		return None


	def recover_relationships(self, a1_dict):
		'''
		Use a2 files to recover the relationships found by TEES

		Here, we end with a dictionary containing all the relationships 
		(by name), and an extra category specifying which relationships
		are included within other rels, and which contain other rels
		'''
		self.obtain_a2_files()
		for art in self.files_list:
			file = open(art).read()
			file_name = art.split('/')[-1]
			self.a2_dictio[file_name] = {}
			self.a2_rels_dict[file_name[:-3]] = {}
			self.a2_tags_dict[file_name[:-3]] = {}
			for line in file.split('\n'):
				if line.startswith('T'):
					tag = line.split('\t')[0]
					entity = line.split('\t')[-1]
					self.a2_dictio[file_name][tag] = entity
				elif line.startswith('E'):
					self.fill_rel_dict(a1_dict, line, file_name)
					self.fill_tag_dict(line, file_name)
		return None


class brenda_annotation():
	'''
	Use Brenda database to annotate the entities found by the 
	NER algorithm
	'''

	def __init__(self):

		self.user = 'salcagal@alumni.uv.es'
		self.pswrd = 'salvacasani91'

	def access_protocol(self):
		'''
		Define the accession protocol information 
		'''
		self.endpointURL = "http://www.brenda-enzymes.org/soap/brenda_server.php"
		self.password = hashlib.sha256(self.pswrd).hexdigest()
		self.client = SOAPProxy(self.endpointURL)
		return None

	def run_function(self, function, ecNumber=None, organism=None, recommendedName=None):
		'''
		Run the brenda function of interest with the required parameters
		'''
		search = ","
		resultString = ''
		if ecNumber != None:
			search += ('ecNumber*' + ecNumber + '#')
		if recommendedName != None:
			search += ('recommendedName*' + recommendedName + '#')
		if organism != None:
			search += ('organism*' + organism + '#')
		parameters = "%s,"%(self.user)+self.password+search
		for n in range(10):
			try:
				resultString = getattr(self.client,function)(parameters)
				break
			except (socket_error, SAXParseException):
				pass
		return resultString


class organism_database():
	'''
	Requirements: Brenda class must have been initialized
	'''

	def __init__(self, organism, db, db_cnx, db_cursor):

		self.organism = organism
		self.synonym2EC = {}
		self.substrate2EC = {}
		self.product2EC = {}
		self.db = db
		self.cnx = db_cnx
		self.cursor = db_cursor

	def obtain_EC(self):
		'''
		Obtain the EC numbers belonging to the desired organism
		'''
		ecNumbers = brenda.run_function('getEcNumbersFromOrganism', organism=self.organism)
		ecNumbers = ecNumbers.split('!')
		return ecNumbers

	def fill_dictionary(self, dictio, entityType, EC, inpu):
		'''
		Fill a dictionary with the information gathered using the Brenda API

		- If used ever again, double check that the 'if entry' clause has the 
		  proper intendation
		'''
		if inpu != '':
			for elem in inpu.split('!'):
				for ent in elem.split('#'):
					if ent.startswith(entityType):
						entry = ent.split('*')[1]
			if entry not in dictio.keys():
				dictio[entry] = set([EC])
			else:
				dictio[entry].add(EC)
		return dictio

	
	def parse_brenda_output(self, input_string, entity_type):
		'''
		Brenda output is a complex string, here the string is parsed
		and a list with all the relevant elements is outputed
		'''
		entity_set = set()
		if input_string != '':
			for elem in input_string.split('!'):
				for entity in elem.split('#'):
					if entity.startswith(entity_type):
						entry = entity.split('*')[1]
						entity_set.add(entry)
		return entity_set

	
	def fill_database(self, EC, substrates, products, synonyms, reactions, recommendedname):
		'''
		Use the brenda information to fill the database. The database is created and
		filled using the manage database class

		- This can be modified by creating a function that can be used to
		  insert a compound to the db.
		'''
		# Add synonym to synonym database
		for synonym in (synonyms|recommendedname):
			self.db.insert_synonym(synonym, EC, self.cnx, self.cursor)
		reac_list = []

		try: 
			recommendedname = recommendedname.pop()
		except KeyError:
			recommendedname = None
		EC_added = False

		# Generate a reactions list
		for reaction in reactions:
			reac_set = set()
			for reac in reaction.split('='):
				for elem in reac.split('+'):
					if elem.strip() != '':
						reac_set.add(elem.strip())
			reac_list.append((reaction, reac_set))
		added_reac = []

		# Add substrate to database
		for substrate in substrates:
			self.db.insert_compound(substrate, EC, self.cnx, self.cursor)
			substrate_added = False
			for reaction in reac_list:
				if substrate in reaction[1]:
					self.db.insert_reaction(reaction[0], EC, substrate, 0, 
						recommendedname, self.cnx, self.cursor)
					added_reac.append(reaction)
					substrate_added = True
					EC_added = True
			if substrate_added == False:
				self.db.insert_reaction(None, EC, substrate, 0, 
					recommendedname, self.cnx, self.cursor)
				EC_added = True
		
		# Add product to database
		for product in products:
			self.db.insert_compound(product, EC, self.cnx, self.cursor)
			product_added = False
			for reaction in reac_list:
				if product in reaction[1]:
					self.db.insert_reaction(reaction[0], EC, product, 1, 
						recommendedname, self.cnx, self.cursor)
					added_reac.append(reaction)
					product_added = True
					EC_added = True
			if product_added == False:
				self.db.insert_reaction(None, EC, product, 1, 
					recommendedname, self.cnx, self.cursor)
				EC_added = True

		# Add the possible reactions that weren't added yet
		for reaction in reac_list:
			if reaction not in added_reac:
				print reaction
				self.db.insert_reaction(reaction[0], EC, None, None, 
					recommendedname, self.cnx, self.cursor)
				EC_added = True

		if EC_added == False:
			self.db.insert_reaction(None, EC, None, None, 
					recommendedname, self.cnx, self.cursor)
		return None

	def generate(self, brenda):
		'''
		Iterate through the different Enzyme Codes to obtain the substrate, 
		product and sysnonyms corresponding to each EC
		'''
		ecNumbers = self.obtain_EC()
		for EC in ecNumbers:
			print EC
			substrates = brenda.run_function('getSubstrate', organism=self.organism, ecNumber=EC)
			substrates = self.parse_brenda_output(substrates, 'substrate')
			products = brenda.run_function('getProduct', organism=self.organism, ecNumber=EC)
			products = self.parse_brenda_output(products, 'product')
			synonyms = brenda.run_function('getSynonyms', organism=self.organism, ecNumber=EC)
			synonyms = self.parse_brenda_output(synonyms, 'synonyms')
			reactions = brenda.run_function('getReaction', organism = self.organism, ecNumber=EC)
			reactions = self.parse_brenda_output(reactions, 'reaction')
			recommendedName = brenda.run_function('getRecommendedName', ecNumber=EC)
			recommendedName = self.parse_brenda_output(recommendedName, 'recommendedName')
			self.fill_database(EC, substrates, products, synonyms, reactions, recommendedName)
		return None


class manage_database():
	'''
	Set of functions to connect, create, and modify a database
	'''

	def __init__(self, username, password, db_name):

	 	self.username = username
	 	self.password = password
	 	self.db_name = db_name
	 	self.tables = {}


	def connect_mysql(self):
		'''
		Connect using username and password
		'''
		cnx = mysql.connector.connect(user=self.username, password=self.password)
		cursor = cnx.cursor(buffered=True)
		return cnx, cursor


	def create_database(self, cursor):
		'''
		Create the database that will be needed to store all the info
		of the BRENDA database
		'''
		try:
			cursor.execute("CREATE DATABASE {} DEFAULT CHARACTER SET 'utf8'".format(self.db_name))
		except mysql.connector.Error as err:
			print("Failed creating database: {}".format(err))
			exit(1)
		return None

	def connect_database(self):
		'''
		Connect to the given database, create it if it does not exist.
		'''
		cnx, cursor = self.connect_mysql() # log in credentials
		try:
			cnx.database = self.db_name
			database_exists = True
		except mysql.connector.Error as err:
			if err.errno == errorcode.ER_BAD_DB_ERROR:
				self.create_database(cursor)
				cnx.database = self.db_name
				database_exists = False
			else:
				print(err)
				exit(1)
		return cnx, cursor, database_exists


	def define_tables(self):
		'''
		Define the tables that will be used to store all the information
		'''
		self.tables['Reaction_table'] = (
		"CREATE TABLE `reaction_table` ("
		"  `EC` varchar(50) NOT NULL,"
		"  `Recommended_name` varchar(300),"
		"  `Compound_ID` varchar(500),"
		"  `S_P` INT,"
		"  `Reaction` VARCHAR(1000)"
		") ENGINE=InnoDB")

		self.tables['name_tbl'] = (
		"CREATE TABLE `name_tbl` ("
		"  `EC_name` varchar(200) NOT NULL,"
		"  `EC` varchar(50) NOT NULL"
		") ENGINE=InnoDB")

		self.tables['compound_tbl'] = (
		"CREATE TABLE `compound_tbl` ("
		"  `Compound_ID` varchar(500) NOT NULL,"
		"  `EC` varchar(50) NOT NULL"
		") ENGINE=InnoDB")
		return None


	def generate_tables(self, cursor):
		'''
		Generate the tables defined, if the table already exists, a
		warning is shown.
		'''
		self.define_tables()
		for name, ddl in self.tables.iteritems():
			try:
				print("Creating table {}: ".format(name))
				cursor.execute(ddl)
			except mysql.connector.Error as err:
				if err.errno == errorcode.ER_TABLE_EXISTS_ERROR:
					print("already exists.")
				else:
					print(err.msg)
			else:
					print("OK")
		return None


	def insert_reaction(self, reaction, ec, compound, s_p, recom_name, cnx, cursor):
		'''
		Insert a reaction into the reaction table
		'''
		new_reaction = ("INSERT INTO reaction_table "
			"(EC, Recommended_name, Compound_ID, S_P, Reaction) "
			"VALUES (%s, %s, %s, %s, %s)")
		data_reaction = (ec, recom_name, compound, s_p, reaction)
		cursor.execute(new_reaction, data_reaction)
		cnx.commit()
		return None


	def insert_compound(self, compound, ec, cnx, cursor):
		'''
		Insert a compound into the compound table
		'''
		new_compound = ("INSERT INTO compound_tbl "
			"(Compound_ID, EC) "
			"VALUES (%s, %s)")
		data_compound = (compound, ec)
		cursor.execute(new_compound, data_compound)
		cnx.commit()
		return None


	def insert_synonym(self, enzyme_name, ec, cnx, cursor):
		'''
		Insert a synonym into the synonym table
		'''
		new_enzyme = ("INSERT INTO name_tbl "
			"(EC_name, EC) "
			"VALUES (%s, %s)")
		data_enzyme = (enzyme_name, ec)
		cursor.execute(new_enzyme, data_enzyme)
		cnx.commit()
		return None


	def extract_reaction(self, EC, cursor):
		'''
		Extract a reaction from the mySQL database using
		the enzyme code as a key
		'''
		obtain_reaction = ("SELECT * FROM reaction_table "
			"WHERE EC = %s")
		EC = (EC)
		cursor.execute(obtain_reaction, (EC,))
		extracted = cursor.fetchall()
		return extracted


	def extract_compound(self, compound, cursor):
		'''
		Extract a compound from the mySQL database, using the 
		compound table
		'''
		obtain_compound = ("SELECT * FROM compound_tbl "
			"WHERE Compound_ID = %s")
		compound = (compound)
		cursor.execute(obtain_compound, (compound,))
		extracted = cursor.fetchall()
		return extracted


	def extract_synonym(self, name, cursor):
		'''
		Extract the enzyme code from the names table using 
		the synonym name as a key
		'''
		obtain_name = ("SELECT * FROM name_tbl "
			"WHERE EC_name = %s")
		synonym = (name)
		cursor.execute(obtain_name, (synonym,))
		extracted = cursor.fetchall()
		return extracted


class annotate_entities():
	'''
	Annotate each of the entities, if they are found as a name
	they will be annotated with the corresponding EC number. 
	If it is a compound, they will be stored in a separate dictionary;
	which will contain the EC number in which reaction it has been 
	found at

	- When we manage to change the annotation of each of the entities,
	  we will need to change this class, distinguishing between 
	  protein and metabolite
	'''

	def __init__(self, db):

		self.names = {}
		self.compounds = {}
		self.non_annotated = []
		self.enzymes = []
		self.db = db
		self.enzymes_dictionary = {}
		self.graph_TM = {}
		self.graph_entities = []


	def fill_dictionary(self, cursor_out, dictio, entity):
		'''
		Use the cursor output to feed the desired dictionary
		'''
		list_of_ECs = []
		for elem in cursor_out:
			list_of_ECs.append(elem[1])
		dictio[entity] = list_of_ECs
		return list_of_ECs


	def annotate(self, entity, cursor):
		'''
		Use the information from the mySQL database to annotate
		each of the entities.
		First we will try to annotate it with the names table,
		and if it does not succeed, we will annotate it using
		the compounds table
		'''
		name = db.extract_synonym(entity, cursor)
		if name == []:
			name = db.extract_compound(entity, cursor)
			if name != []:
				self.fill_dictionary(name, self.compounds, entity)
			else:
				self.non_annotated.append(entity)
		else:
			enz = self.fill_dictionary(name, self.names, entity)
			for e in enz:
				self.enzymes.append(e)
		return self.enzymes


	def fill_enzymes_dict(self, substrates, products, enzyme, enzymeName):
		'''
		Use input info to fill the general enzymes dict, 
		that will be used to generate the SIF file that 
		will give the information to Cytoscape.
		'''
		if enzyme not in self.enzymes_dictionary.keys():
			self.enzymes_dictionary[enzyme] = {'substrates':[], 
			'products':[], 'name':''}
			self.enzymes_dictionary[enzyme]['substrates'] = substrates
			self.enzymes_dictionary[enzyme]['products'] = products
			self.enzymes_dictionary[enzyme]['name'] = enzymeName
		return None

	
	def enzyme_info(self, enzyme, cursor):
		'''
		Obtain the information from the reaction covered by the 
		given enzyme using the mySQL database
		'''
		substrates = []; products = []
		reactions = self.db.extract_reaction(enzyme, cursor)
		for reaction in reactions:
			if reaction[2] != None:
				if reaction[3] == 0:
					substrates.append(reaction[2])
				elif reaction[3] == 1:
					products.append(reaction[2])
		recommended_name = reactions[0][1]
		return substrates, products, recommended_name


	def enzymes_from_compounds(self, cursor):
		'''
		Obtain the enzymes that correspond to each compound
		if they are present in the enzymes found. 
		Use these enzymes as the center of the graph, and 
		use the compounds as a mean to create links between 
		the enzymes
		'''
		for compound, enzymes in self.compounds.iteritems():
			enzyme_found = False
			subs_enzymes = []; prod_enzymes = []
			for enzyme in enzymes:
				if enzyme in self.enzymes:
					enzyme_found = True
					substrates, products, enzyme_name = self.enzyme_info(enzyme, cursor)
					self.fill_enzymes_dict(substrates, products, enzyme, enzyme_name)
			if enzyme_found == True:
				#self.fill_enzymes_dict(subs_enzymes, prod_enzymes, compound)
				#Create compound object?
				pass
		return None

	def compounds_from_enzymes(self, cursor):
		'''
		Don't limit the results to the compounds found, but
		select all the enzymes detected and obtain the 
		Use these enzymes as the center of the graph, and 
		use the compounds as a mean to create links between 
		the enzymes
		'''
		for enzyme in self.enzymes:
			substrates, products, enzyme_name = self.enzyme_info(enzyme, cursor)
			self.fill_enzymes_dict(substrates, products, enzyme, enzyme_name)
		return None


	def graph_annotation(self, graphDB, entity):
		'''
		Use the graph database to annotate the entities
		'''
		result = graphDB.is_enzyme(entity)
		res = [i for i in result]
		if len(res) == 0:
			result = graphDB.is_compound(entity)
			res = [i for i in result]
		if len(res) > 0:
			self.graph_TM[entity] = []
			for elem in res:
				name = elem[0]
				label = elem[1]
				self.graph_TM[entity].append(name)
				self.graph_entities.append(name.encode("utf-8"))
		return None


class create_graph_file():
	'''
	Create the SIF file that will be used as input for 
	Cytoscape. This file will be created from 
	'''

	def __init__(self, enzymes_dictionary):
		
		self.enzymes_dictionary = enzymes_dictionary
		self.enzymes_used = set()
		self.graph_dictionary = {}


	def fill_dictio(self, entity, successor):
		'''
		Fill a dictionary, in which the keys can be either
		compounds or enzymes, and they can be related to 
		either compounds or enzymes
		'''
		if entity in self.graph_dictionary.keys():
			self.graph_dictionary[entity].add(successor)
		else:
			self.graph_dictionary[entity] = set([successor])
		return None


	def create_file(self, compounds_dict):
		'''
		Use the enzymes dictionary to create the file
		'''
		for enzyme in self.enzymes_dictionary.keys():
			for substrate in self.enzymes_dictionary[enzyme]['substrates']:
				has_enzyme = False
				try:
					compound_enzymes = compounds_dict[substrate]
				except KeyError:
					continue
				for enz in compound_enzymes:
					if enz in self.enzymes_dictionary.keys():
						has_enzyme = True
						break
				if has_enzyme == True:
					self.fill_dictio(substrate, enzyme)
					self.enzymes_used.add(enzyme)
			for product in self.enzymes_dictionary[enzyme]['products']:
				try:
					compound_enzymes = compounds_dict[product]
				except KeyError:
					continue
				has_enzyme = False
				for enz in compound_enzymes:
					if enz in self.enzymes_dictionary.keys():
						has_enzyme = True
						break
				if has_enzyme == True:
					self.fill_dictio(enzyme, product)
					self.enzymes_used.add(enzyme)
		return None


	def write_file(self, file):
		'''
		Write the file using the graph dictionary
		'''
		myfile = open(file, 'w')
		myfile.write('in_node'+'\t'+'Relationship'+'\t'+'out_node'+'\n')
		for k, v in self.graph_dictionary.iteritems():
			try:
				k = self.enzymes_dictionary[k]['name']
			except KeyError:
				pass
			k = k.replace(' ', '_')
			for element in v:
				try:
					element = self.enzymes_dictionary[element]['name']
				except KeyError:
					pass
				element = element.replace(' ', '_')
				myfile.write(k + '\t')
				myfile.write('Brenda_relationship' + '\t')
				myfile.write(element + '\n')
		myfile.close()
		return None


class TM_to_dict:
	'''
	Create a dictionary that includes the information from the
	a2 dictionary created from the TEES result.
	'''

	def __init__(self, min_entities=2):

		self.min_entities = min_entities
		self.causal_graph_dictio = {}
		self.no_causal_graph_dictio = {}


	def entities_count(self, art, rel, rels_dictio):
		'''
		Count the number of entities that appear at each 
		of the relationships covered
		'''
		entities_number = 0
		if rel['Theme'] != []:
			for entity in rel['Theme']:
				if re.match(r'^E\d{1}', entity) is not None:
					n = self.entities_count(art, rels_dictio[art][entity], rels_dictio)
					entities_number += n
				else:
					entities_number += 1
		if rel['Cause'] != []:
			for entity in rel['Cause']:
				if re.match(r'^E\d{1}', entity) is not None:
					n = self.entities_count(art, rels_dictio[art][entity], rels_dictio)
					entities_number += n
				else:
					entities_number += 1
		return entities_number


	def obtain_entities(self, art, rel, rels_dictio):
		'''
		Obtain the entities from a 
		'''
		causal_entities = []; theme_entities = []
		if rel['Cause'] != []:
			for elem in rel['Cause']:
				if re.match(r'^E\d{1}', elem) is not None:
					causal, theme = self.obtain_entities(art, rels_dictio[art][elem], rels_dictio)
					for ent in (causal+theme):
						causal_entities.append(ent)
				else:
					causal_entities.append(elem)
		if rel['Theme'] != []:
			for elem in rel['Theme']:
				if re.match(r'^E\d{1}', elem) is not None:
					causal, theme = self.obtain_entities(art, rels_dictio[art][elem], rels_dictio)
					for ent in (causal + theme):
						theme_entities.append(ent)
				else:
					theme_entities.append(elem)
		return causal_entities, theme_entities


	def fill_graph_dictionary(self, art, rel, rels_dictio, sentence):
		'''
		Use the input relationship to fill a dictionary, 
		which will be used to reconstruct the graph
		'''
		causal, theme = self.obtain_entities(art, rel, rels_dictio)
		rel_type = rel['Rel']
		if causal != []:
			for elem in causal:
				if elem not in self.causal_graph_dictio.keys():
					self.causal_graph_dictio[elem] = {}
				for them in theme:
					if them not in self.causal_graph_dictio[elem].keys():
						self.causal_graph_dictio[elem][them] = {'Rel_types':[], 'n_occurrences':0, 'sentence':[]}
					for relat in rel_type:
						self.causal_graph_dictio[elem][them]['Rel_types'].append(relat)
					self.causal_graph_dictio[elem][them]['n_occurrences'] += 1
					self.causal_graph_dictio[elem][them]['sentence'].append(sentence)
		else:
			for elem in theme:
				if elem not in self.no_causal_graph_dictio.keys():
					self.no_causal_graph_dictio[elem] = {}
				for them in theme:
					if them != elem:
						if them not in self.no_causal_graph_dictio[elem].keys():
							self.no_causal_graph_dictio[elem][them] = {'Rel_types':[], 'n_occurrences':0, 'sentence':[]}
						for relat in rel_type:
							self.no_causal_graph_dictio[elem][them]['Rel_types'].append(relat)
						self.no_causal_graph_dictio[elem][them]['n_occurrences'] += 1
						self.no_causal_graph_dictio[elem][them]['sentence'].append(sentence)
		return None


	def obtain_index(self, art, index_dictio, tag, tag_rels):
		'''
		Obtain the index of the relationship. I take a random index from 
		all the entities that are part of the relationship
		'''
		theme, cause = self.obtain_entities(art, tag_rels[art][tag], tag_rels)
		entities = theme + cause
		for entity in entities:
			index = index_dictio[art][entity]
		return index


	def obtain_sentence(self, art, sentences_dictio, index):
		'''
		Obtain the sentences from the sentences dictionary by using the indexes
		'''
		for key, value in sentences_dictio[art].iteritems():
			index1 = int(key.split('-')[0]); index2 = int(key.split('-')[1])
			if index1 <=int(index[0]) <= index2:
				sentence = value
		return sentence


	def process_relationships(self, rels_dictio, index_dict, tag_rels, sentences_dict):
		'''
		Iterate through all the relationships, and reconstruct
		the 'interaction' that will be added to the graph.
		It needs to recover both entities and relationships
		'''
		for art, rels in rels_dictio.iteritems():
			self.rels_used = []
			for tag, rel in rels.iteritems():
				enti_num = self.entities_count(art, rel, rels_dictio)
				if enti_num >= self.min_entities:
					index = self.obtain_index(art, index_dict, tag, tag_rels)
					sentence = self.obtain_sentence(art, sentences_dict, index)
					self.fill_graph_dictionary(art, rel, rels_dictio, sentence)
		return None


class TMdict_to_graph:
	'''
	Use the information gathered from the a2 files to create relationships
	between the different entities. These entities are matched against the
	annotated entities from the a1 file dictionary.
	'''

	def __init__(self, graph_file, ent_an, graph_class):

		self.out_file = graph_file
		self.ent_an = ent_an
		self.graph_class = graph_class


	def fill_graph_file(self, ent1, ent2):
		'''
		Use the info from the TM progress to fill the graph, first
		check if the relationship is already in the dictio
		'''
		ent1 = ent1.replace(' ', '_')
		ent2 = ent2.replace(' ', '_')
		file = open(self.out_file, 'a')
		file.write(ent1 + '\t')
		file.write('TM_relationships' + '\t')
		file.write(ent2 + '\n')
		file.close()
		return None


	def check_relationship(self, ent1, ent2):
		'''
		Check if the relationship was already added to the graph
		First we need to retrieve the enzyme code from the entity
		'''
		try:
			ec = self.ent_an.names[ent1]
			for e in ec:
				if e in self.ent_an.enzymes_dictionary.keys():
					ent1 = e
					break
		except KeyError:
			pass
		try:
			ec = self.ent_an.names[ent2]
			for e in ec:
				if e in self.ent_an.enzymes_dictionary.keys():
					ent2 = e
					break
		except KeyError:
			pass

		rel_graph = False
		if ent1 in self.graph_class.graph_dictionary.keys():
			if ent2 in self.graph_class.graph_dictionary[ent1]:
				rel_graph = True
		return rel_graph


	def check_dictionaries(self, entity):
		'''
		Check if the entity is either within the enzymes or the
		compounds that were annotated in the entities.
		True means that is in the enzymes or compounds dictionary 
		'''
		try:
			ec = self.ent_an.names[entity]
			is_entity = True
			for e in ec:
				if e in self.ent_an.enzymes_dictionary.keys():
					name = self.ent_an.enzymes_dictionary[e]['name']
					break
		except KeyError:
			try:
				ent_an.compounds[entity]
				is_entity = True
				name = entity
			except KeyError:
				is_entity = False
				name = None
		return is_entity, name


	def cause_dictio2graph(self, cause_dictio):
		'''
		Open the graph file and modify it to include the
		relationships summarized in the causal graph. 
		These relationships are unidirectional.
		'''
		for key, value in cause_dictio.iteritems():
			key_entity, kName = self.check_dictionaries(key)
			for v in value.keys():
				v_entity, vName = self.check_dictionaries(v)
				if key_entity and v_entity == True:
					k_v = self.check_relationship(key, v)
					v_k = self.check_relationship(v, key)
					if not k_v or v_k == True:
						self.fill_graph_file(kName, vName)
		return None


	def non_cause_dictio2graph(self, non_cause):
		'''
		Open the graph file and modify it to include the relationships 
		present in the non causal graph. These relationships are 
		bidirectional, because we do not have any information
		about causality.
		'''
		for key, value in non_cause.iteritems():
			key_entity, kName = self.check_dictionaries(key)
			for v in value.keys():
				v_entity, vName = self.check_dictionaries(v)
				if key_entity and v_entity == True:
					k_v = self.check_relationship(key, v)
					v_k = self.check_relationship(v, key)
					if not k_v or v_k == True:
						self.fill_graph_file(kName, vName)
		return None


class graph_database():
	'''
	Create a Neo4j graph database, create the tools to fill it
	'''


	def __init__(self, brenda, organism):

		self.brenda = brenda
		self.organism = organism


	def connect(self):
		'''
		Connect to the graph db using neo4j credentials
		'''
		self.driver = GraphDatabase.driver("bolt://localhost:7687", 
			auth=basic_auth("neo4j", "salva"))
		self.session = self.driver.session()
		return None


	def create_enzyme(self, uniID, uniName, syns_list, ec):
		'''
		Add an enzyme entity to the database
		'''
		self.session.run("CREATE (a:Enzyme:Protein {uniprotID: {uniID}, synonyms: {syns}, "
			"ECs: {ec}, id: {uniID}, uniprotEntryName: {uniEntry}, specie: {specie}})",
			{"uniID": uniID, "syns": syns_list, "ec": [ec], "uniEntry": uniName, "specie":self.organism})
		return None


	def set_enzyme_properties(self, uniID, syns_list, ecs):
		'''
		Use an existing node and change its variable properties
		'''
		self.session.run("MATCH (a) WHERE a.id = {uniID} SET a.synonyms = {synonyms}, a.ECs = {ECs}",
			{"uniID": uniID, "synonyms": syns_list, "ECs":ecs})
		return None


	def add_prop_enzyme(self, uniID, syns_list, ec):
		'''
		Use an existing enzyme and add new properties.
		'''
		result = self.session.run("MATCH (a) WHERE a.id = {uniID} RETURN a.synonyms as "
			"synonyms, a.ECs as EC", {"uniID":uniID})
		for n in result:
			synonyms = set(n["synonyms"])
			EC = set(n["EC"])
			for syn in syns_list:
				synonyms.add(syn)
			EC.add(ec)
			self.set_enzyme_properties(uniID, list(synonyms), list(EC))
		return None


	def create_ec_enzyme(self, EC, syns_list):
		'''
		Create an enzyme node
		Used for the EC-centered case
		'''
		self.session.run("CREATE (a:Enzyme:Protein {synonyms: {syns}, "
			"ECs: {ec}, id: {ec}, specie: {specie}})",
			{"syns": syns_list, "ec": [EC], "specie":self.organism})
		return None


	def add_prop_ec_enzyme(self, EC, syns_list):
		'''
		Update an enzyme node
		Used for the EC-centered case
		'''
		result = self.session.run("MATCH (a) WHERE a.id = {EC} RETURN a.synonyms as "
			"synonyms, a.ECs as EC", {"EC":EC})
		for n in result:
			synonyms = set(n["synonyms"])
			ec = set(n["EC"])
			for syn in syns_list:
				synonyms.add(syn)
			ec.add(EC)
			self.set_enzyme_properties(EC, list(synonyms), list(ec))
		return None


	def create_compound(self, chebiID, chebiName):
		'''
		Add a compound to the database
		'''
		self.session.run("CREATE (a:Compound {chebiID: {chebiID}, "
			"compoundName: {compoundName}, id: {chebiID}})",
			{"chebiID": chebiID, "compoundName": chebiName})
		return None


	def create_brenda_relationship(self, nodeA, nodeB, ec, reaction, specie):
		'''
		Create a relationship between to existing nodes in the database
		'''
		self.session.run('MATCH (a), (b) WHERE a.id={nodeA} AND b.id={nodeB} '
			"CREATE (a)-[:Brenda_relationship {ECs: {ec}, reactionsBrenda: {reaction}, species: {specie}}]"
			"->(b)",{"nodeA":nodeA, "nodeB":nodeB, "ec":[ec], "reaction":[reaction], "specie":[specie]})
		return None


	def check_relationship(self, nodeA, nodeB, relationship):
		'''
		Check if a relationship already exists
		'''
		result = self.session.run('MATCH (a)-[r:%s]->(b) '
			'WHERE a.id="%s" AND b.id="%s" RETURN a,r,b'%(relationship,nodeA, nodeB))
		if len([i for i in result]) > 0:
			is_rel = True
		else:
			is_rel = False
		return is_rel


	def update_brenda_relationship(self, nodeA, nodeB, ec, reaction, specie):
		'''
		Update the information from ane xisting relationship
		ec, reaction and specie must be a list
		'''
		result = self.session.run('MATCH (a)-[r:Brenda_relationship]->(b) '
			'WHERE a.id={nodeA} AND b.id={nodeB} '
			'RETURN r.ECs as ecs, r.reactionsBrenda as reac, r.species as species',
			{"nodeA":nodeA, "nodeB":nodeB})
		for n in result:
			ecs = (set(n["ecs"]+ec))
			reactions = (set(n["reac"])|set(reaction))
			species = (set(n["species"]+specie))
			self.set_reaction_properties(nodeA, nodeB, list(ecs), list(reactions), list(species))
		return None


	def set_reaction_properties(self, nodeA, nodeB, ec, reaction, specie):
		'''
		Create a relationship between to existing nodes in the database
		all properties musy be list type
		'''
		self.session.run('MATCH (a)-[r:Brenda_relationship]->(b) WHERE a.id={nodeA} AND b.id={nodeB} '
			"SET r.ECs = {ec}, r.reactionsBrenda = {reaction}, r.species = {specie}",
			{"nodeA":nodeA, "nodeB":nodeB, "ec":ec, "reaction":reaction, "specie":specie})
		return None


	def create_prop_relationship(self, nodeA, nodeB, relationship, sentence):
		'''
		Create the relationship and add the sentence information
		'''
		self.session.run('MATCH (a), (b) WHERE a.name="%s" AND b.name="%s" '
			"CREATE (a)-[:%s {Sentence: %s}]->(b)"%(nodeA, nodeB, relationship, sentence))
		return None


	def check_protein(self, ID):
		'''
		Check if a node exists. Will be used before adding a node
		'''
		result = self.session.run("MATCH (a:Protein) WHERE a.id={id} "
			"RETURN a.id AS id", 
			{"id": ID})
		if len([i for i in result]) > 0:
			is_node = True
		else:
			is_node = False
		return is_node


	def check_compound(self, ID):
		'''
		Check if a node exists. Will be used before adding a node
		'''
		result = self.session.run("MATCH (a:Compound) WHERE a.id={id} "
			"RETURN a.id AS id", 
			{"id": ID})
		if len([i for i in result]) > 0:
			is_node = True
		else:
			is_node = False
		return is_node


	def is_enzyme(self, syn):
		'''
		Check if the entity is into the synonyms of any enzyme
		'''
		result = self.session.run("MATCH (a) WHERE '%s' IN a.Synonyms "
			"RETURN a.name, labels(a)"%(syn))
		return result


	def is_compound(self, compound):
		'''
		Check if the entity is included within any compound 
		(node) name
		'''
		result = self.session.run('MATCH (a) WHERE a.name = "%s" '
			'RETURN a.name, labels(a)'%(compound))
		return result


	def brenda_obtain_EC(self):
		'''
		Obtain the EC numbers belonging to the desired organism
		'''
		ecNumbers = self.brenda.run_function('getEcNumbersFromOrganism', organism=self.organism)
		ecNumbers = ecNumbers.split('!')
		return ecNumbers


	def parse_brenda_output(self, input_string, entity_type):
		'''
		Brenda output is a complex string, here the string is parsed
		and a list with all the relevant elements is outputed
		'''
		entity_set = set()
		if input_string != '':
			for elem in input_string.split('!'):
				reaction = None
				for entity in elem.split('#'):
					if entity.startswith(entity_type):
						entry = entity.split('*')[1]
					if entity.startswith('reactionPartners'):
						reaction = entity.split('*')[1]
				if entity_type in ('substrate', 'product'):
					entity_set.add((entry,reaction))
				else:
					entity_set.add(entry)
		return entity_set


	def fill_ec_centered_graphdb(self, EC, substrates, products, synonyms):
		'''
		Use Brenda information to recover the relationship between
		ECs and metabolites. In this case the graph will be ec centered
		'''
		is_node = self.check_protein(EC)
		if is_node == True:
			self.add_prop_ec_enzyme(EC, list(synonyms))
		else:
			self.create_ec_enzyme(EC, list(synonyms))
			
		for subs in substrates:
			substrate = subs[0]; reaction = subs[1]
			chebi_name = chebi_from_string(chebi_dictio, substrate)
			chid, name = chebi_name.chebi()
			if chid != None:
				is_node = self.check_compound(chid)
				if is_node == False:
					self.create_compound(chid, name)
				is_rel = self.check_relationship(chid, EC, 'Brenda_relationship')
				if is_rel == False:
					self.create_brenda_relationship(chid, EC, EC, reaction, self.organism)
				else:
					self.update_brenda_relationship(chid, EC, [EC], [reaction], [self.organism])

		for prod in products:
			product = prod[0]; reaction = prod[1]
			chebi_name = chebi_from_string(chebi_dictio, product)
			chid, name = chebi_name.chebi()
			if chid != None:
				is_node = self.check_compound(chid)
				if is_node == False:
					self.create_compound(chid, name)
				is_rel = self.check_relationship(EC, chid, "Brenda_relationship")
				if is_rel == False:
					self.create_brenda_relationship(EC, chid, EC, reaction, self.organism)
				else:
					self.update_brenda_relationship(EC, chid, [EC], [reaction], [self.organism])
		return None
	



	def fill_graph_db(self, ec, substrates, products, synonyms, ec2uniprot, chebi_dictio):
		'''
		Use the Brenda information to fill the graph database, create a node
		for each protein that correspond to the enzyme and establish the relationships 
		according to the connection by compound
		'''
		for prot in ec2uniprot[ec]:
			uniID = prot[0]; uniName = prot[1]
			is_node = self.check_protein(uniID)
			if is_node == True:
				self.add_prop_enzyme(uniID, list(synonyms), ec)
			else:
				self.create_enzyme(uniID, uniName, list(synonyms), ec)

			for subs in substrates:
				substrate = subs[0]; reaction = subs[1]
				chebi_name = chebi_from_string(chebi_dictio, substrate)
				chid, name = chebi_name.chebi()
				if chid != None:
					is_node = self.check_compound(chid)
					if is_node == False:
						self.create_compound(chid, name)
					is_rel = self.check_relationship(chid, uniID, 'Brenda_relationship')
					if is_rel == False:
						self.create_brenda_relationship(chid, uniID, ec, reaction, self.organism)
					else:
						self.update_brenda_relationship(chid, uniID, [ec], [reaction], [self.organism])

			for prod in products:
				product = prod[0]; reaction = prod[1]
				chebi_name = chebi_from_string(chebi_dictio, product)
				chid, name = chebi_name.chebi()
				if chid != None:
					is_node = self.check_compound(chid)
					if is_node == False:
						self.create_compound(chid, name)
					is_rel = self.check_relationship(uniID, chid, "Brenda_relationship")
					if is_rel == False:
						self.create_brenda_relationship(uniID, chid, ec, reaction, self.organism)
					else:
						self.update_brenda_relationship(uniID, chid, [ec], [reaction], [self.organism])
		return None


	def create_database(self, ec2uniprot, chebi_dictio):
		'''
		Just as it was done for the MySQL database, use the BRENDA database
		information to fill the Neo4j database. The main difference with
		the mySQL database is that here we add the relationships as well
		'''
		self.connect()
		ecNumbers = self.brenda_obtain_EC()
		for EC in ecNumbers:
			print EC
			substrates = brenda.run_function('getSubstrate', organism=self.organism, ecNumber=EC)
			substrates = self.parse_brenda_output(substrates, 'substrate')
			products = brenda.run_function('getProduct', organism=self.organism, ecNumber=EC)
			products = self.parse_brenda_output(products, 'product')
			synonyms = brenda.run_function('getSynonyms', organism=self.organism, ecNumber=EC)
			synonyms = self.parse_brenda_output(synonyms, 'synonyms')
			recommendedName = brenda.run_function('getRecommendedName', ecNumber=EC)
			recommendedName = self.parse_brenda_output(recommendedName, 'recommendedName')
			synonyms = synonyms|recommendedName
			try:
				ec2uniprot[EC]
				self.fill_graph_db(EC, substrates, products, synonyms, ec2uniprot, chebi_dictio)
			except KeyError:
				if len(ec2uniprot.keys()) == 0:
					print 'alo'
					if len(substrates) != 0:
						self.fill_ec_centered_graphdb(EC, substrates, products, synonyms)
				pass
		return None


	def check_db_length(self):
		'''
		Check the length of the graph db, this will help to take the
		decission of creating it or using the existing db
		'''
		nodes = None
		driver = GraphDatabase.driver("bolt://localhost:7687", 
			auth=basic_auth("neo4j", "salva"))
		session = driver.session()
		result = session.run('MATCH (n) RETURN count(n)')
		for r in result:
			for n in r:
				nodes = r[0]
		session.close()
		return nodes


	def check_TM_relationship(self, ent1, ent2):
		'''
		Check if a given relationship exists
		'''
		result = self.session.run('MATCH (n)-[r]-(y)' 
			' WHERE n.name = "%s" AND y.name = "%s" RETURN n.name, y.name'%(ent1, ent2))
		res = [i for i in result]
		return res


	def add_TM_check(self, key, value, graph_entity):
		'''
		Iterate over the list of keys and values, check if each of 
		these relationships exist, if the don't add them to the
		Neo4j graph DB
		'''
		if key in graph_entity.keys():
			myKeys = graph_entity[key]
			for v in value.keys():
				sentence = set(value[v]['sentence'])
				if v in graph_entity.keys():
					myVals = graph_entity[v]
					for element in myKeys:
						for target in myVals:
							rel = self.check_TM_relationship(element, target)
							if len(rel) == 0:
								#Modify this statement to include properties to the
								#relationship, this will allow the addition of the text
								#Also need to add the type of relationship (binding, negative, ...)
								self.create_prop_relationship(element, target, 'TM_relationship', sentence)
			#Maybe we need to include a property that states that
			#for the same relationship a lot of different entities 
			#(synonyms) were used
		return None



	def add_TM_relationship(self, cause_dictio, non_cause_dictio, graph_entity):
		'''
		Add the relationships gathered from the TEES text mining
		into the Neo4j graph DB.
		'''
		for key, value in cause_dictio.iteritems():
			self.add_TM_check(key, value, graph_entity)
		for key, value in non_cause_dictio.iteritems():
			self.add_TM_check(key, value, graph_entity)
		return None


class retrieve_graph():

	def __init__(self, graph_database, enzymes):

		self.graph_database = graph_database
		self.graph_entities = enzymes


	def obtain_graph(self):
		'''
		Here I will have to implement a way to access the graph DB 
		directly, to show the results
		'''
		return None


	def accession_query(self):
		'''
		Generate the query that will be used to access the Neo4j 
		graph DB
		'''
		self.graph_entities.remove('ADP')
		self.graph_entities.remove('ATP')
		self.graph_entities.remove('CoA')
		exceptions = ['ATP', 'ADP', 'more', '?', 'H2O', 'CoA']
		#All the nodes
		print '\nRepresent all nodes'
		print 'MATCH (n)'
		print 'WITH %s AS entities, n'%(self.graph_entities)
		print 'WHERE n.name IN entities' 
		print 'RETURN n'%(self.graph_entities)
		#All the connected nodes
		print '\nRepresent all connected nodes'
		print 'MATCH p = (n)-[r]-(y)'
		print 'WITH %s AS entities, p'%(self.graph_entities)
		print 'WHERE n.name IN entities AND y.name IN entities' 
		print 'RETURN p'
		#All the connected nodes allowing for one "space"
		print '\nRepresent all connected nodes, allowing for one intermediate node'
		print 'MATCH p = (n)--(y)--(t)'
		print 'WITH %s AS entities, %s AS exceptions, p'%(self.graph_entities, exceptions)
		print 'WHERE n.name IN entities AND NOT y.name IN exceptions AND t.name IN entities AND n <> t'
		print 'RETURN p'
		return None


class chebi_from_string(): #TODO: test to see what we get in function of the input string. If input string== chebi name then keep chebi. Clean function to compare both?
   
	def __init__(self,chebi_dictio,name=None):
		'''
		search in the chebi database a chebi number corresponding to the chebi name
		'''
		self.chebi_dictio = chebi_dictio
		self.long=name
		self.chebiid=""
    
	def chebi(self):
		'''
		http://docs.python-zeep.org/en/master/
		http://www.ebi.ac.uk/chebi/webServices.do
		'''
		chebiId=None; chebiName=None
		if self.long in self.chebi_dictio.keys():
			chebiId = self.chebi_dictio[self.long][0]
			chebiName = self.chebi_dictio[self.long][1]
		else:
			wsdl='http://www.ebi.ac.uk/webservices/chebi/2.0/webservice?wsdl'
			client = zeep.Client(wsdl=wsdl)
			for n in range(10):
				try:
					listchebiname=client.service.getLiteEntity(self.long,"ALL NAMES",5,'ALL')
					break
				except (zeep.exceptions.Fault, zeep.exceptions.TransportError):
					listchebiname = None
					pass
			if listchebiname!=None:
				for entity in listchebiname:
					chebiId = entity.chebiId
					chebiName = entity.chebiAsciiName
					if self.long == chebiName:
						self.chebi_dictio[self.long] = [chebiId, chebiName]
						return chebiId, chebiName
				if listchebiname[0].searchScore >= 10:
					chebiId = listchebiname[0].chebiId
					chebiName = listchebiname[0].chebiAsciiName
			self.chebi_dictio[self.long] = [chebiId, chebiName]
		return chebiId, chebiName #return also chebi ascii name


def clean(word,specie=''):
    '''
    Clean words to allow them comparison.
    See if I need to add all the cases that I covered with perl before or if the score in the chebi comparison is enouth
    '''
    #lowercases
    cleanw=word.lower()
    specie=specie.lower()
    sp=specie.split(' ')
    specie=specie.strip()
    #no whitecharacters
    cleanw=cleanw.strip()
    #remove parenthesis
    cleanw=re.sub('^\(\+?\-?\d{0,5}\)','',cleanw)
    cleanw=re.sub('\(','',cleanw)
    cleanw=re.sub('\)','',cleanw)
    cleanw=re.sub('\[','',cleanw)
    cleanw=re.sub('\]','',cleanw)
    cleanw=re.sub('\{','',cleanw)
    cleanw=re.sub('\}','',cleanw)
    cleanw=re.sub(' genes{0,1}','',cleanw)
    cleanw=re.sub(' protein{0,1}','',cleanw)
    #cleanw=re.sub('\)','',cleanw)
    #remove words never used in the text
    cleanw=re.sub(' atom','',cleanw)
    cleanw=cleanw.replace(' ','') #test 3 july
    for s in sp:
        cleanw=re.sub(s,'',cleanw)#remove the specie name from the entity name (ex for Arabidopsis thaliana: remove arabidopsis and then remove thaliana)
    if cleanw=='':
        cleanw='empty_string'
    return cleanw


if __name__ == '__main__':



	start = time.time()

	
	########################################################################
	## Before running this script, the user needs to download the desired 
	## papers, and run TEES and BANNER on them. This should be implemented
	## within the script. 
	## The mysql database has to be installed, as well as the neo4j DB
	## in case it needs to be used.
	## Still need to implement the options menu for the algorithm to run.
	########################################################################


	try:
		path = argv[1]
		entity = entity_recognition(path)
		rel = TM_relationship(path)
		sent = sentence_retriever(path)
	except IndexError:
		entity = entity_recognition()
		rel = TM_relationship()
		sent = sentence_retriever()

	organism = 'Homo sapiens'; organism_tax = '9606'

	proteins, metabolites, a1_dictio = entity.create_entity()
	rel.recover_relationships(a1_dictio)
	sent.get_sentence()



	# for art, elems in rel.a2_rels_dict.iteritems():
	# 	print art
	# 	for k,v in elems.iteritems():
	# 		print k, v
	


	######################################################################
	## Create a mySQL database to store the relationships through all the 
	## Enzyme Codes. The keys to look through the db would be the substrates, 
	## compounds and synonyms.
	##
	## Brenda annotation has to be initialized, otherwise the SOAP 
	## credentials are not loaded
    ##
    ## This module requires mySQL database to be install, with the give 
    ## user and password
	######################################################################



	username = 'salvacasani'; password = 'salvacasani'
	db_name = 'homo_sapiens'

	cytoscape_graph = 'graph.SIF'
	
	#Describe the db management class

	db = manage_database(username, password, db_name)
	cnx, cursor, db_exists = db.connect_database()

	organisms = open('speclist.txt').read()
	ec_converter = open('enzyme.dat').read()
	uniprot = uniprot_matcher(organisms)
	uniprot.EC2Uniprot(ec_converter, organism_tax)

	if db_exists == False:

		print 'Generating %s database, this step may take a couple of hours\n'%(db_name)

		db.generate_tables(cursor)         # Generate database if it is not created, generate tables if they are not created

		# Obtain the desired brenda annotation information
		brenda = brenda_annotation()
		brenda.access_protocol()

		hSapiens = organism_database(organism, db, cnx, cursor)
		hSapiens.generate(brenda)
	
	else:

		print '\nDatabase exists, advancing to analysis step'


	########################################################################
	## Now the database is completed, and we have a class to interact with
	## it, now we can use it to annotate the entities found.
	########################################################################


	ent_an = annotate_entities(db)

	for protein in proteins:
		enz = ent_an.annotate(protein, cursor)

	ent_an.enzymes_from_compounds(cursor)
	#ent_an.compounds_from_enzymes(cursor)


	########################################################################
	## We can use the annotated entities to generate the graph. We can do it
	## using the annotated compounds or using the enzymes.
	########################################################################


	graph = create_graph_file(ent_an.enzymes_dictionary)
	graph.create_file(ent_an.compounds)
	graph.write_file(cytoscape_graph)


	########################################################################
	## Once the 'preliminary' graph is created, we can use the relationships
	## recovered by TEES using the Text Mining algorithms in order to add
	## extra information to the reactions described via BRENDA
	########################################################################


	graph_rels = TM_to_dict()
	graph_rels.process_relationships(rel.a2_rels_dict, entity.index_dictio, 
		rel.a2_tags_dict, sent.sentence_dict)


	
	########################################################################
	## Use the dictionaries created from the Text Mining steps, and 
	## include this information in the graph file that will be introduced
	## in cytoscape. We need to take into account that this dictionary
	## will contain the info for all the corpora that were used.
	########################################################################


	TM2graph = TMdict_to_graph(cytoscape_graph, ent_an, graph)
	TM2graph.cause_dictio2graph(graph_rels.causal_graph_dictio)
	TM2graph.non_cause_dictio2graph(graph_rels.no_causal_graph_dictio)



	# print ent_an.names
	# print ent_an.compounds
	# print ent_an.non_annotated
	# print ent_an.enzymes
	# print ent_an.enzymes_dictionary


	########################################################################
	## Use the brenda API to recover all information from BRENDA DB and 
	## introduce it into the graph database (Neo4j)
	########################################################################


	brenda = brenda_annotation()
	brenda.access_protocol()


	
	graph_db = graph_database(brenda, organism)
	nodes_number = graph_db.check_db_length()
	
	chebi_dictio = {}

	if nodes_number == 0:
		print '\nCreating Neo4j graph database'
		graph_db.create_database(uniprot.ec2Uniprot_dict, chebi_dictio)
	else:
		print '\nExisting Neo4j Graph Database, annotating entities'
		graph_db.connect()

	# Will need to remove this line!!
	#graph_db.create_database(uniprot.ec2Uniprot_dict, chebi_dictio)

	for protein in proteins:
		ent_an.graph_annotation(graph_db, protein)


	graph_db.add_TM_relationship(graph_rels.causal_graph_dictio, 
		graph_rels.no_causal_graph_dictio, ent_an.graph_TM)


	retrieve_graph = retrieve_graph(graph_db, ent_an.graph_entities)
	retrieve_graph.accession_query()


	cursor.close()
	cnx.close()
	graph_db.session.close()
	
	end = time.time()
	print end - start