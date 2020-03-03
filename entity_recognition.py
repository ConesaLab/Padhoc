#!/usr/bin/env python

import time

from chebi_from_string import chebi_from_string
from parse_files import parse_file
from parse_tmchem_file import parse_tmchemfile
from libchebipy import ChebiEntity
import re, sys, itertools, neo4j
from graph_database import graph_database
from brenda_annotation import brenda_annotation
from uniprot_queries import uniprot_queries
from compare_string_to_db import compare_string_to_db
from neo4j import exceptions


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
    cleanw=re.sub('<[^>]+>','',cleanw)
    cleanw=re.sub('^(r)-','',cleanw)
    cleanw=re.sub('^(s)-','',cleanw)
    cleanw=re.sub('\)n$','',cleanw)
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
        cleanw='empty_string'
    return cleanw

class entity_recognition():

	def __init__(self, chebi, uniprot, abreviations, db, sentence_dict, protein_synonyms, species):

		self.chebi = chebi
		self.uniprot = uniprot
		self.db = db
		self.a1_dictio = {}
		self.a2_dictio = {}
		self.abreviations = abreviations
		self.chebiname = {} # To check if the metaboite was already tested, but relies on the info from the requests
		self.sentenceDict = sentence_dict
		self.protein_synonyms = protein_synonyms
		self.a1_prot_uniprot = {}
		self.species = species
		self.protein_dictio = {}
		self.fill_dictionaries_from_neo4j()
		self.false_entities = []


	def process_a1(self, path, pmid, dictmchem, dicoPMID_idDoc, query):
		'''
		Process the a1 file, where the entities are called.
		'''
		self.a1_dictio[pmid] = {}
		path = path + '/' + pmid + '.a1'
		with open(path) as a1_file:
			for line in a1_file:
				line = line.rstrip().split('\t')
				self.a1_dictio[pmid][line[0]]={}
				self.a1_dictio[pmid][line[0]]['textname'] = line[-1].replace('"','')
				typecmp = line[1].split(' ')
				if typecmp[0]=='Metabolite':
					self.a1_dictio[pmid][line[0]]['chebi'] = None
					self.get_metabolite(line, dictmchem, dicoPMID_idDoc, pmid)
					self.add_metabolite(line, pmid, query)
				else:
					self.a1_dictio[pmid][line[0]]['id']=''
					self.a1_dictio[pmid][line[0]]['uniprotID'] = ''
					self.get_protein(line, pmid)
					self.add_protein(query, pmid, line)
		return None


	def process_a1_changed(self, path, pmid, dictmchem, dicoPMID_idDoc, query, compareDB):
		'''
		Process the a1 file, where the entities are called.
		'''
		self.a1_dictio[pmid] = {}
		path = path + '/' + pmid + '.a1'
		start = time.time()
		with open(path) as a1_file:
			for line in a1_file:
				end1 = time.time()
				start = time.time()
				line = line.rstrip().split('\t'); is_protein = False
				self.a1_dictio[pmid][line[0]]={}
				self.a1_dictio[pmid][line[0]]['textname'] = line[-1].replace('"','')
				typecmp = line[1].split(' ')
				# Need to search first the database to doublecheck if the entity is already present there
				if clean(line[-1]) not in self.false_entities:			
					present_in_DB, is_metabolite, is_protein = self.search_entity(line, pmid)
					
					if present_in_DB == False:
						is_metabolite = self.get_metabolite(line, dictmchem, dicoPMID_idDoc, pmid)
						if is_metabolite == False:
							uniID = self.get_protein(line, pmid)
							end_uniprot = time.time() - start
							if uniID != None:
								is_protein = True

						metab = self.node_in_DB(line[0], self.a1_dictio[pmid][line[0]]['textname'], compareDB, pmid, is_metabolite, is_protein)
						if metab == True:
							is_metabolite = True
						elif metab == False:
							is_protein = True
					if is_metabolite == False and is_protein == False:
						self.false_entities.append(clean(line[-1]))
				else:
					is_metabolite = False
					is_protein = False

				if is_metabolite == False and is_protein == False:
					if typecmp[0] == 'Metabolite':
						is_metabolite = True
					else:
						is_protein = True

				if is_metabolite == True:
					self.add_metabolite(line, pmid, query)
				if is_protein == True:
					self.add_protein(query, pmid, line)
		return None

	
	def node_in_DB(self, element, entity, compareDB, pmid, metabolite, protein):
		'''
		Find if the name is present in the DB using the similarity function
		Fill in the nodes dictionary with the information that will be used

		If the metabolite has already been found, it adds the info into the
		existing dictionary
		'''
		start = time.time()
		nodes, is_metabolite = compareDB.similar_nodes(entity)
		if len(nodes) > 0:
			if is_metabolite == True:
				if metabolite == True:
					for node in nodes:
						if node not in self.a1_dictio[pmid][element]['chebi']:
							self.a1_dictio[pmid][element]['chebi'].append(node)
						if node not in self.chebi.nameid[clean(entity)]:
							self.chebi.nameid[clean(entity)].append(node)
				# The DB_comparison is more likely to fail than the API
				elif protein == True:
					is_metabolite = False
#					print entity
#					print nodes
				else:
					self.a1_dictio[pmid][element]['chebi'] = nodes
					self.chebi.nameid[clean(entity)] = nodes
			elif is_metabolite == False:
				if protein == True:
					for node in nodes:
						self.a1_dictio[pmid][element]['uniprotID'].append(node)
						if clean(entity) not in self.a1_prot_uniprot.keys():
							self.a1_prot_uniprot[clean(entity)] = {}
							self.a1_prot_uniprot[clean(entity)]['uniprotID'] = [node]
						else:
							if node not in self.a1_prot_uniprot[clean(entity)]['uniprotID']:
								self.a1_prot_uniprot[clean(entity)]['uniprotID'].append(node)
				# The DB_comparison is more likely to fail than the API
				elif metabolite == True:
#					print entity
#					print nodes
					is_metabolite = True
				else:
					self.a1_dictio[pmid][element]['uniprotID'] = nodes
					self.a1_prot_uniprot[clean(entity)] = {}
					self.a1_prot_uniprot[clean(entity)]['uniprotID'] = nodes

		# If it wasn't found in the DB, add the info searching through the API
		if is_metabolite == None:
			if metabolite == True:
				is_metabolite = True
			elif protein == True:
				is_metabolite = False
		return is_metabolite


	def obtain_neo4j(self):
		'''
		Obtain neo4j db using a cypher request. Here we
		will only ask for the properties that we are interested,

		'''
		request = 'MATCH (n:Compound) RETURN n.chebiID as id, n.compoundName as name'
		compoundResult = self.db.session.run(request)
		return compoundResult


	def fill_dictionaries_from_neo4j(self):
		'''
		Use the neo4j database to recover all the information from
		the nodes and fill the dictionaries from this module with them,
		thus extracting plenty of information without the need of 
		requesting to the APIs
		'''
		compResult = self.obtain_neo4j()
		for r in compResult:
			chebiID = r['id']
			chebiName = r['name']
			if chebiID not in self.chebiname.keys():
				self.chebiname[chebiID] = chebiName
		return None
	

	def check_chebi(self, chebi):
		'''
		Check if a chebi entity is already in the chebi dictio,
		otherwise add it
		'''
		if chebi not in self.chebiname.keys():
			chebi_name = ChebiEntity(chebi).get_name()
			self.chebiname[chebi] = chebi_name
		return None


	def check_chebiClass(self, pmid, metabolite):
		'''
		Check the chebi class dictionaries for other occurrences
		of these metabolites 
		'''
		tmchemadded = False
		try:
			chebi_ids = self.chebi.nameid[clean(metabolite[-1])]
			for chebi in chebi_ids:
				self.check_chebi(chebi)
				try:
					self.a1_dictio[pmid][metabolite[0]]['chebi'].append(chebi)
				except KeyError:	
					self.a1_dictio[pmid][metabolite[0]]['chebi'] = [chebi]
			tmchemadded = True
		except KeyError:
			#Check if the name is found with an abreviation
			if pmid in self.abreviations.keys():
				if metabolite[-1] in self.abreviations[pmid].keys():
					metab_name = self.abreviations[pmid][metabolite[-1]]
					try:
						chebi_ids = self.chebi.nameid[clean(metab_name)]
						for chebi in chebi_ids:
							self.check_chebi(chebi)
							try:
								self.a1_dictio[pmid][metabolite[0]]['chebi'].append(chebi)
							except KeyError:
								self.a1_dictio[pmid][metabolite[0]]['chebi'] = [chebi]
							try:
								self.chebi.nameid[clean(metabolite[-1])].append(chebi)
							except KeyError:
								self.chebi.nameid[clean(metabolite[-1])] = [chebi]
						tmchemadded = True
					except KeyError:
						pass
		return tmchemadded


	def check_tmchem_result(self, metabolite, dictmchem, pmidtmchem, pmid):
		'''
		Check the result from tmchem when the TM algorithm was run
		'''
		tmchemAdded = False
		# Check if there is a tmchem recognition in the article
		if pmidtmchem[pmid][0] in dictmchem.keys():
		#Check if the metabolite in question is recognized by tmchem
			if metabolite[-1] in dictmchem[pmidtmchem[pmid][0]]:
				chebi = dictmchem[pmidtmchem[pmid][0]][metabolite[-1]]
				self.check_chebi(chebi)
				self.a1_dictio[pmid][metabolite[0]]['chebi'] = [chebi]
				self.a1_dictio[pmid][metabolite[0]]['chebiName'] = self.chebiname[chebi]
				try:
					self.chebi.nameid[clean(metabolite[-1])].append(chebi)
				except KeyError:
					self.chebi.nameid[clean(metabolite[-1])] = [chebi]
				tmchemAdded = True
		return tmchemAdded


	def request_to_chebi(self, metabolite, pmid):
		'''
		Request the name to the chebi database
		'''
		foundChebi = False
		# Include a way to check if the entity was already tested
		chebiId, chebiName = self.chebi.chebi(metabolite[-1])
		#print chebiId, chebiName
		if chebiId != None:
			self.a1_dictio[pmid][metabolite[0]]['chebi'] = [chebiId]
			self.a1_dictio[pmid][metabolite[0]]['chebiName'] = [chebiName]
			self.chebiname[chebiId] = chebiName
			try:
				self.chebi.nameid[clean(metabolite[-1])].append(chebiId)
			except KeyError:
				self.chebi.nameid[clean(metabolite[-1])] = [chebiId]
			foundChebi = True
		#Check also the abreviation
		else:
			if pmid in self.abreviations.keys():
				if metabolite[-1] in self.abreviations[pmid].keys():
						metab_name = self.abreviations[pmid][metabolite[-1]]
						metab_abreviation = [metabolite[0], metabolite[1], metab_name]
						foundChebi = self.request_to_chebi(metab_abreviation, pmid)
		#print chebiId, chebiName
		return foundChebi


	def get_metabolite(self, metabolite, dictmchem, pmidtmchem, pmid):
		'''
		Metabolite normalization function, done using ChEBI
		'''
		#Check if the ChEBI id was found with tmChem
		tmchemadded = self.check_tmchem_result(metabolite, dictmchem, pmidtmchem, pmid)
		if tmchemadded == False:
			#Check if it is already in the chebi dictionary
			tmchemadded = self.check_chebiClass(pmid, metabolite)
		# If none of the previous methods worked, make a request
		# to the chebi database
		if tmchemadded == False:
			tmchemadded = self.request_to_chebi(metabolite, pmid)
		return tmchemadded


	def create_new_metabolite(self, pmId, sentence, textName, query):
		'''
		Create a new metabolite into the graph database
		'''
		request = 'MERGE (n:Compound {id:"%s",textname: ["%s"],sentences: ["%s"],PMID_Tnb: ["%s"]})' %(pmId,textName,sentence,query)
		self.db.session.run(request)
		return None


	def search_entity(self, entity, pmid):
		'''
		Search the entity in the proteins and metabolite dictionary.
		In case it exists add it to the a1 dictionary.
		'''
		entity_name = clean(entity[-1])
		search_result = False
		try:
			chebi_ids = self.chebi.nameid[entity_name]
			self.a1_dictio[pmid][entity[0]]['chebi'] = chebi_ids
			search_result = True; is_metabolite = True
		except KeyError:
			is_metabolite = False
		try:
			uniIDs = self.a1_prot_uniprot[entity_name]
			self.a1_dictio[pmid][entity[0]]['uniprotID'] = uniIDs
			is_protein = True; search_result = True
		except KeyError:
			is_protein = False
		return search_result, is_protein, is_metabolite


	def add_sentence(self, dbId, sentence, r, entityType):
		'''
		Add a sentence to the entity of interest
		'''
		if r['s'] != None:
			if sentence not in r['s']:
				request = 'MATCH (n:%s) WHERE n.id="%s" SET n += {sentences: n.sentences + "%s"}'%(entityType, dbId, sentence)
				self.db.session.run(request)
		else:
			request = 'MATCH (n:%s) WHERE n.id="%s" SET n.sentences = ["%s"]'%(entityType, dbId, sentence)
			self.db.session.run(request)
		return None


	def add_textName(self, dbId, textName, r, entityType):
		'''
		Add the text name to the entity of interest
		'''
		if r['tn'] != None:
			if textName not in r['tn']:
				request = 'MATCH (n:%s) WHERE n.id="%s" SET n += {textname: n.textname + "%s"}'%(entityType, dbId, textName)
				self.db.session.run(request)
		else:
			request = 'MATCH (n:%s) WHERE n.id="%s" SET n.textname = ["%s"]'%(entityType, dbId, textName)
			self.db.session.run(request)
		return None


	def add_query(self, dbId, query, r, entityType):
		'''
		Add the PmID to the entity of interest
		'''
		if r['pmid'] != None:
			if query not in r['pmid']:
				request = 'MATCH (n:%s) WHERE n.id="%s" SET n += {PMID_Tnb: n.PMID_Tnb + "%s"}'%(entityType, dbId, query)
				self.db.session.run(request)
		else:
			request = 'MATCH (n:%s) WHERE n.id="%s" SET n.PMID_Tnb = ["%s"]'%(entityType, dbId, query)
			self.db.session.run(request)
		return None


	def add_existing_metabolite(self, dbId, query, sentence, textName):
		'''
		Match an existing metabolite, and add the information that
		is lacking into the entity in the graph database
		'''
		request = 'MATCH (n:Compound) WHERE n.id="%s" RETURN n.sentences as s,'\
		' n.textname as tn, n.PMID_Tnb as pmid'%(dbId)
		result = self.db.session.run(request)
		for r in result:
			self.add_sentence(dbId, sentence, r, 'Compound')
			self.add_textName(dbId, textName, r, 'Compound')
			self.add_query(dbId, query, r, 'Compound')
		return None


	def add_metabolite(self, metabolite, pmid, query):
		'''
		Add the metabolite to the database
		'''
		graphQuery = query + '_' + pmid + ':' + metabolite[0]
		sentence = pmid + ':' + self.sentenceDict[pmid][metabolite[0]]
		#print self.a1_dictio[pmid][metabolite[0]].keys(), self.a1_dictio[pmid][metabolite[0]]
		if 'chebi' not in self.a1_dictio[pmid][metabolite[0]].keys():
			textname = self.a1_dictio[pmid][metabolite[0]]['textname']
			dbId = clean(textname)
			self.a1_dictio[pmid][metabolite[0]]['id'] = [dbId]
			#textName = self.a1_dictio[pmid][metabolite[0]]['textname']
			result = self.db.session.run('MATCH (n:Compound) WHERE n.id="%s" RETURN n.id'%(dbId))
			if len([i for i in result]) > 0:
				self.add_existing_metabolite(dbId, graphQuery, sentence, textname)
			else:
				self.create_new_metabolite(dbId, sentence, textname, graphQuery)
		else:
			dbIds = self.a1_dictio[pmid][metabolite[0]]['chebi']
			#chebiName = self.a1_dictio[pmid][metabolite[0]]['chebiName']
			self.a1_dictio[pmid][metabolite[0]]['id'] = dbIds
			for dbId in dbIds:
				chebiName = self.chebiname[dbId]
				textName = self.a1_dictio[pmid][metabolite[0]]['textname']
				result = self.db.session.run('MATCH (n:Compound) WHERE n.id="%s" RETURN n.id'%(dbId))
				if len([i for i in result]) > 0:
					self.add_existing_metabolite(dbId, graphQuery, sentence, textName)
				else:
					self.create_new_metabolite(dbId, sentence, textName, graphQuery)
					self.db.session.run('MATCH (n:Compound) WHERE n.id="%s" SET n.compoundName = "%s", n.chebiID="%s"'%(dbId, chebiName, dbId))
		return None


	def get_protein(self, protein, pmid):
		'''
		Get the protein name and obtain the Uniprot ID for the protein of interest.
		'''
		uniID = None; uniEntry = None
		geneNames = None; protNames = None
		if clean(protein[-1]) in self.protein_synonyms.keys():
			uniID = self.protein_synonyms[clean(protein[-1])]
			self.a1_dictio[pmid][protein[0]]['uniprotID'] = [uniID]
			if clean(protein[-1]) not in self.a1_prot_uniprot.keys():
				self.a1_prot_uniprot[clean(protein[-1])] = {}
				self.a1_prot_uniprot[clean(protein[-1])]['uniprotID'] = [uniID]
			if uniID not in self.protein_dictio.keys():
				self.protein_dictio[uniID] = {}
				# Need to extract the infor from the uniprot ID
				self.protein_dictio[uniID]['uniprotEntry'] = ''
				self.protein_dictio[uniID]['uniprotgenenames'] = ''
				self.protein_dictio[uniID]['uniprotprotnames'] = ''
		else:
			if clean(protein[-1]) in self.a1_prot_uniprot.keys():
				self.a1_dictio[pmid][protein[0]]['uniprotID'] = self.a1_prot_uniprot[clean(protein[-1])]['uniprotID']
				uniID = self.a1_dictio[pmid][protein[0]]['uniprotID']
			else:
				uniID, uniEntry, geneNames, protNames = self.uniprot.query_id(protein[-1])
				#print uniID, uniEntry, geneNames, protNames
				if uniID != '':
					self.a1_dictio[pmid][protein[0]]['uniprotID'] = [uniID]
					if uniID not in self.protein_dictio.keys():
						self.protein_dictio[uniID] = {}
						self.protein_dictio[uniID]['uniprotEntry'] = uniEntry
						self.protein_dictio[uniID]['uniprotgenenames'] = geneNames
						self.protein_dictio[uniID]['uniprotprotnames'] = protNames
					self.a1_prot_uniprot[clean(protein[-1])] = {}
					self.a1_prot_uniprot[clean(protein[-1])]['uniprotID'] = [uniID]
		if pmid in self.abreviations.keys():
			if uniID == '':
				#print self.abreviations[pmid].keys(), '\n\n\n\n\n\n'
				if protein[-1] in self.abreviations[pmid].keys():
					prot_name = self.abreviations[pmid][protein[-1]]
					protein_abreviation = [protein[0], protein[1], prot_name]
					uniID = self.get_protein(protein_abreviation, pmid)
		if uniID == '':
			uniID = None
		return uniID


	def create_new_protein(self, dbId, sentence, textName, query, graphQuery, protein, uniprot, pmid):
		'''
		Add a new protein node to the graph database. It depends wether the protein
		has been found on uniprot or not
		'''
		if uniprot == False:
			request = 'MERGE (a:Protein {id:"%s",name: "%s",textname: \
			["%s"], sentences: ["%s"], PMID_Tnb: ["%s"], specie: "%s", query: ["%s"]})'\
			%(dbId, textName, textName, sentence, graphQuery, self.species, query)
			try:
				self.db.session.run(request)
			except neo4j.exceptions.CypherSyntaxError:
				request = 'MERGE (a:Protein {id:"%s",name: "%s",textname: \
				["%s"], sentences: ["%s"], PMID_Tnb: ["%s"], specie: "%s", query: ["%s"]})'\
				%(dbId, textName, textName, "No sentence to show", graphQuery, self.species, query)
				self.db.session.run(request)
		else:
			uniEntry = self.protein_dictio[dbId]['uniprotEntry']
			geneNames = self.protein_dictio[dbId]['uniprotgenenames']
			protNames =  self.protein_dictio[dbId]['uniprotprotnames']
			geneList = []; protList = []
			geneNames = geneNames.split(' ')
			for gene in geneNames:
				geneList.append(re.sub(':', '', gene).encode('utf8'))
			for prot in protNames.split('('):
				protList.append(re.sub('\)', '', prot).encode('utf8'))
			request = 'MERGE (a:Protein {id:"%s",name: "%s",textname: \
			["%s"], sentences: ["%s"], PMID_Tnb: ["%s"], specie: "%s",\
			query: ["%s"], uniprotID: "%s", uniProtEntryName: "%s", uniprotGenesNames: %s,\
			 uniprotProteinNames: %s})'\
			%(dbId, textName, textName, sentence, graphQuery, self.species, \
			 query, dbId, uniEntry, geneList, protList)
			self.db.session.run(request)
		return None


	def add_existing_protein(self, dbId, query, sentence, textName, uniprot):
		'''
		Match an existing metabolite, and add the information that
		is lacking into the entity in the graph database
		'''
		request = 'MATCH (n:Protein) WHERE n.id="%s" RETURN n.sentences as s,'\
		' n.textname as tn, n.PMID_Tnb as pmid'%(dbId)
		result = self.db.session.run(request)
		for r in result:
			self.add_sentence(dbId, sentence, r, 'Protein')
			self.add_textName(dbId, textName, r, 'Protein')
			self.add_query(dbId, query, r, 'Protein')
		return None


	def add_protein_node(self, dbId, sentence, query, graphQuery, protein, uniprot, pmid):
		'''
		Check if the node is already in the db, if it is not add it, 
		otherwise update the existing node
		'''
		result = self.db.session.run('MATCH (n:Protein) WHERE n.id="%s" RETURN n.id' %(dbId))
		textName = self.a1_dictio[pmid][protein[0]]['textname']
		if len([i for i in result]) > 0:
			self.add_existing_protein(dbId, query, sentence, textName, uniprot)
		else:
			self.create_new_protein(dbId, sentence, textName, query, graphQuery, protein, uniprot, pmid)
		return None


	def add_protein(self, query, pmid, protein):
		'''
		Add the protein to the Neo4j database using all the information recovered
		'''
		graphQuery = query + '_' + pmid + ':' + protein[0]
		sentence = pmid + ':' + self.sentenceDict[pmid][protein[0]]
		sentence.replace('\\', '')
		if 'uniprotID' in self.a1_dictio[pmid][protein[0]].keys():
			#If there is self.sentenceDictno uniprot ID to use as ID in the db
			self.a1_dictio[pmid][protein[0]]['id'] = self.a1_dictio[pmid][protein[0]]['uniprotID']
			for uniId in self.a1_dictio[pmid][protein[0]]['uniprotID']:
				uniprot = True
				dbId = uniId
				self.add_protein_node(dbId, sentence, query, graphQuery, protein, uniprot, pmid)
		else:
			#If there is uniprot ID to use as ID in the db
			uniprot = False
			textname = self.a1_dictio[pmid][protein[0]]['textname']
			dbId = clean(textname)
			self.a1_dictio[pmid][protein[0]]['id'] = [dbId]
			self.add_protein_node(dbId, sentence, query, graphQuery, protein, uniprot, pmid)
		return None


	def get_E_nodes(self, pmid, entity, final_nodes):
		'''
		Obtain the nodes that are enclosed by the E entity
		'''
		nodes = self.a2_dictio[pmid][entity]['nodes']
		for node in nodes:
			if node.startswith('T'):
				final_nodes.add(node)
			elif node.startswith('E'):
				final_nodes = self.get_E_nodes(pmid, node, final_nodes)
		return final_nodes


	def add_a2dictio(self, theme, pmid, E_entities):
		'''
		Add the relationship to the a2 dictionary, which 
		contains the entities that form the relationship
		'''
		nodes = []; relationships = []; nodes_number = 0
		for t in theme:
			entity = t.split(':')[1]
			if entity.startswith('T'):
				try:
					self.a1_dictio[pmid][entity]
					nodes.append(entity)
					nodes_number += 1
				except KeyError:
					ent_name = self.a2_dictio[pmid][entity]
					relationships.append(ent_name)
			elif entity.startswith('E'):
				E_nodes = self.get_E_nodes(pmid, entity, set())
				nodes.append(entity)
				nodes_number += len(E_nodes)
				if entity in E_entities:
					E_entities.remove(entity)
		return nodes, relationships, nodes_number, E_entities


	def create_a2_ditio(self, fileName, pmid):
		'''
		Use the file to create the a2 dictionary
		'''
		E_entities = [] # E entities that should be used to create the graph
		with open(fileName) as file:
			for line in file:
				entities = []
				line = line.rstrip()
				if line.startswith('T'):
					line = line.split('\t')
					self.a2_dictio[pmid][line[0]] = line[-1]
				elif line.startswith('E'):
					line = line.split('\t')
					theme = line[1].split(' ')
					nodes, relationships, node_number, E_entities = self.add_a2dictio(theme, pmid, E_entities)
					self.a2_dictio[pmid][line[0]] = {'nodes':nodes, 'relationships':relationships}
					if node_number > 1:
						E_entities.append(line[0])
		return E_entities

	
	def new_relationship(self, nodeA, nodeB, search_Type, relationships, query, sentence):
		'''
		Target two existing nodes and create a new relationship between them
		'''
		request = 'MATCH (n), (y) WHERE n.id="%s" AND y.id = "%s" CREATE (n)-[r:TM_relationship {species:["%s"],\
		 reactionTypes: %s, sentences:%s, query:["%s"], nbs: 1, corpora:["%s"]}]->(y)'\
		 %(nodeA, nodeB, self.species, relationships, sentence, query, search_Type)
		result = self.db.session.run(request)
		return None


	def update_relationship(self, nodeA, nodeB, prop, value):
		'''
		Update the required information from the relationship in neo4j
		'''
		request = 'MATCH (n)-[r]->(y) WHERE n.id="%s" AND y.id="%s"\
		 SET r.%s=%s'%(nodeA, nodeB, prop, value)
		self.db.session.run(request)
		return None

	def modify_relationship_properties(self, nodeA, nodeB, prop, value):
		'''
		Modify a value in a relationship
		'''
		request = 'MATCH (n)-[r]->(y) WHERE n.id="%s" AND y.id="%s" RETURN r.%s as prop'\
		 %(nodeA, nodeB, prop)
		result = self.db.session.run(request)
		for r in result:
			if prop != 'nbs':
				if not isinstance(value, list):
					value = [value]
				if r['prop'] != None:
					properties = [elem.encode('utf-8') for elem in r['prop']]
					for elem in value:
						if elem not in properties:
							properties.append(elem)
				else:
					properties = value
				self.update_relationship(nodeA, nodeB, prop, properties)
			else:
				properties = r['prop']
				if r['prop'] != None:
					number = r['prop'] + 1
				else:
					number = 1
				self.update_relationship(nodeA, nodeB, prop, number)
		return None


	def modify_relationship(self, nodeA, nodeB, search_Type, relationships, query, sentence):
		'''
		Target an existing relationship and modify its properties
		'''
		self.modify_relationship_properties(nodeA, nodeB, 'reactionTypes', relationships)
		self.modify_relationship_properties(nodeA, nodeB, 'sentences', sentence)
		self.modify_relationship_properties(nodeA, nodeB, 'query', query)
		self.modify_relationship_properties(nodeA, nodeB, 'corpora', search_Type)
		self.modify_relationship_properties(nodeA, nodeB, 'species', self.species)
		self.modify_relationship_properties(nodeA, nodeB, 'nbs', None)
		return None


	def add_relationship_gd(self, nodes, relationships, pmid, search_Type, query):
		'''
		Add the relationship to the database, first the
		entity has to be mathced with an existing node
		'''
		node_groups = [r for r in itertools.combinations(nodes,2)]
		for nods in node_groups:
			nodeA_id = self.a1_dictio[pmid][nods[0]]['id']
			nodeB_id = self.a1_dictio[pmid][nods[1]]['id']
			for nodeA in nodeA_id:
				for nodeB in nodeB_id:
					request = 'MATCH (n)-[r:TM_relationship]-(y) WHERE n.id="%s" AND\
					 y.id="%s" RETURN r as r'%(nodeA, nodeB)
					result = self.db.session.run(request)
					sentence = list(set([self.sentenceDict[pmid][nods[0]], self.sentenceDict[pmid][nods[1]]]))
					if len([r for r in result]) > 0:
						self.modify_relationship(nodeA, nodeB, search_Type, relationships, query, sentence)
					else:
						self.new_relationship(nodeA, nodeB, search_Type, relationships, query, sentence)
		return None


	def create_relationship(self, pmid, nodes, relationships, search_Type, query):
		'''
		Create set the nodes that will be linked together
		'''
		if any(node.startswith('E') for node in nodes):
			new_nodes = []; new_rels = relationships
			for node in nodes:
				if node.startswith('E'):
					for nod in self.a2_dictio[pmid][node]['nodes']:
						new_nodes.append(nod)
					for rel in self.a2_dictio[pmid][node]['relationships']:
						new_rels.append(rel)
				else:
					new_nodes.append(node)
			self.create_relationship(pmid, new_nodes, new_rels, search_Type, query)
		else:
			self.add_relationship_gd(nodes, relationships, pmid, search_Type, query)
		return None

	def process_a2(self, path, pmid, search_Type, query):
		'''
		Extract the relationships from the a2 file and add 
		them to the graph database
		'''
		fileName = path + '/' + pmid + '.a2'
		self.a2_dictio[pmid] = {}
		E_entities = self.create_a2_ditio(fileName, pmid)
		for entity in E_entities:
			nodes = self.a2_dictio[pmid][entity]['nodes']
			relationships = self.a2_dictio[pmid][entity]['relationships']
			self.create_relationship(pmid, nodes, relationships, search_Type, query)
		return None


if __name__ == '__main__':

	
	path = '/home/salva/human_signaling/tees/GE11/a1a2_TTTAQPFT_P1/'
	pmid = '23892279'
	pmid = '26473953'
	nerFile = '/home/salva/human_signaling/tees/GE11/a1a2_TTTAQPFT_P1/-preprocessed.xml.gz-ner.xml.gz'
	tmChem_file = '/home/salva/human_signaling/tees/GE11/a1a2_TTTAQPFT_P1/sentences.pubtator.tmChem'
	repository = '/home/salva/human_signaling/tees/GE11/a1a2_TTTAQPFT_P1/'

	#nerFile = '/home/salva/human_signaling/metrecon/a1a2_1HSYGRAR_P2/placeholder-preprocessed.xml.gz-ner.xml.gz'
	#tmChem_file = '/home/salva/human_signaling/metrecon/a1a2_1HSYGRAR_P2/sentences.pubtator.tmChem'
	#repository = '/home/salva/human_signaling/metrecon/a1a2_1HSYGRAR_P2/'

	path = '/home/salva/ecoli_salva_pantothenate/metrecon/a1a2_TTTAQPFT_P1/'
	pmid = '23892279'
	pmid = '26366567'
	nerFile = '/home/salva/ecoli_salva_pantothenate/metrecon/a1a2_TTTAQPFT_P1/placeholder-preprocessed.xml.gz-ner.xml.gz'
	tmChem_file = '/home/salva/ecoli_salva_pantothenate/metrecon/a1a2_TTTAQPFT_P1/sentences.pubtator.tmChem'
	repository = '/home/salva/ecoli_salva_pantothenate/metrecon/a1a2_TTTAQPFT_P1/'

	email = 'salcagal@alumni.uv.es'; brendapass = 'salvacasani91'

	chebi = chebi_from_string()
	chebi.chebi_connect()
	uniprot = uniprot_queries('Homo sapiens', '9606')

	brenda = brenda_annotation(email,brendapass) # create brenda_annotation object
	brenda.access_protocol() # access to brenda
	gd=graph_database(chebi, brenda,'Homo sapiens','neo4j','salva') #
	gd.connect()

	pner = parse_file(nerFile, {}, {}, repository)
	dicoartTsentence, dicoPMID_idDoc = pner.docTsent()

	tm = parse_tmchemfile(tmChem_file, {}, repository)
	dictmchem = tm.parse()

	
	compareDB = compare_string_to_db(gd)

	entity_recognition = entity_recognition_salva(chebi, uniprot, {}, gd, dicoartTsentence, {}, 'Homo sapiens')
	#entity_recognition.a1_dictio = {'26607036': {'T398': {'uniprotID': '', 'textname': 'H2A', 'id': 'h2a'}, 'T399': {'uniprotID': '', 'textname': 'D90', 'id': 'd90'}, 'T390': {'uniprotID': u'P08779', 'textname': 'K16', 'uniprotprotnames': u'Keratin, type I cytoskeletal 16 (Cytokeratin-16) (CK-16) (Keratin-16) (K16)', 'uniprotEntry': u'K1C16_HUMAN', 'id': u'p08779', 'uniprotgenenames': u'KRT16 KRT16A'}, 'T391': {'chebi': None, 'textname': 'H2B'}, 'T392': {'chebi': None, 'textname': 'K'}, 'T393': {'uniprotID': '', 'textname': 'histone', 'id': 'histone'}, 'T394': {'uniprotID': '', 'textname': 'H2A', 'id': 'h2a'}, 'T395': {'uniprotID': '', 'textname': 'H2B', 'id': 'h2b'}, 'T396': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T397': {'uniprotID': '', 'textname': 'H18', 'id': 'h18'}, 'T72': {'uniprotID': '', 'textname': 'H3-K56ac', 'id': 'h3k56ac'}, 'T73': {'uniprotID': '', 'textname': 'histone H3', 'id': 'histoneh3'}, 'T70': {'chebi': 'CHEBI:15356', 'textname': 'cysteine'}, 'T71': {'uniprotID': '', 'textname': 'histone', 'id': 'histone'}, 'T76': {'chebi': None, 'textname': 'acetyllysine'}, 'T77': {'chebi': 'CHEBI:33709', 'textname': 'amino acids'}, 'T74': {'uniprotID': '', 'textname': 'H3-K56', 'id': 'h3k56'}, 'T75': {'uniprotID': '', 'textname': 'chromatin compaction37', 'id': 'chromatincompaction37'}, 'T78': {'uniprotID': '', 'textname': 'protein synthesis42', 'id': 'proteinsynthesis42'}, 'T79': {'uniprotID': u'P62805', 'textname': 'histone H4 proteins', 'uniprotprotnames': u'Histone H4', 'uniprotEntry': u'H4_HUMAN', 'id': u'p62805', 'uniprotgenenames': u'HIST1H4A H4/A H4FA; HIST1H4B H4/I H4FI; HIST1H4C H4/G H4FG; HIST1H4D H4/B H4FB; HIST1H4E H4/J H4FJ; HIST1H4F H4/C H4FC; HIST1H4H H4/H H4FH; HIST1H4I H4/M H4FM; HIST1H4J H4/E H4FE; HIST1H4K H4/D H4FD; HIST1H4L H4/K H4FK; HIST2H4A H4/N H4F2 H4FN HIST2H4; HIST2H4B H4/O H4FO; HIST4H4'}, 'T365': {'uniprotID': '', 'textname': 'H2B-K34', 'id': 'h2bk34'}, 'T364': {'uniprotID': '', 'textname': 'R31', 'id': 'r31'}, 'T367': {'uniprotID': '', 'textname': 'S36', 'id': 's36'}, 'T366': {'uniprotID': '', 'textname': 'H2B', 'id': 'h2b'}, 'T361': {'chebi': None, 'textname': 'N'}, 'T360': {'uniprotID': '', 'textname': 'molecule-D H2B', 'id': 'moleculedh2b'}, 'T248': {'uniprotID': '', 'textname': 'histones H3', 'id': 'histonesh3'}, 'T249': {'uniprotID': '', 'textname': 'H2A', 'id': 'h2a'}, 'T246': {'uniprotID': '', 'textname': 'S3a-h', 'id': 's3ah'}, 'T247': {'chebi': 'CHEBI:33704', 'textname': 'amino acid'}, 'T244': {'uniprotID': '', 'textname': 'NCPs', 'id': 'ncp'}, 'T245': {'uniprotID': u'Q9HD26', 'textname': 'Fig', 'uniprotprotnames': u'Golgi-associated PDZ and coiled-coil motif-containing protein (CFTR-associated ligand) (Fused in glioblastoma) (PDZ protein interacting specifically with TC10) (PIST)', 'uniprotEntry': u'GOPC_HUMAN', 'id': u'q9hd26', 'uniprotgenenames': u'GOPC CAL FIG'}, 'T242': {'chebi': 'CHEBI:36976', 'textname': 'nucleotides'}, 'T243': {'uniprotID': '', 'textname': 'NCPs', 'id': 'ncp'}, 'T240': {'uniprotID': '', 'textname': 'Fig. 4c', 'id': 'fig.4c'}, 'T241': {'chebi': None, 'textname': 'N'}, 'T98': {'uniprotID': '', 'textname': 'PylRS', 'id': 'pylr'}, 'T99': {'uniprotID': '', 'textname': 'KacRS_6mt', 'id': 'kacrs_6mt'}, 'T149': {'chebi': None, 'textname': 'Mg++'}, 'T148': {'uniprotID': '', 'textname': 'H2B', 'id': 'h2b'}, 'T419': {'uniprotID': u'Q8WVN6', 'textname': 'K12', 'uniprotprotnames': u'Secreted and transmembrane protein 1 (Protein K-12)', 'uniprotEntry': u'SCTM1_HUMAN', 'id': u'q8wvn6', 'uniprotgenenames': u'SECTM1 K12'}, 'T418': {'uniprotID': u'P05787', 'textname': 'K8', 'uniprotprotnames': u'Keratin, type II cytoskeletal 8 (Cytokeratin-8) (CK-8) (Keratin-8) (K8) (Type-II keratin Kb8)', 'uniprotEntry': u'K2C8_HUMAN', 'id': u'p05787', 'uniprotgenenames': u'KRT8 CYK8'}, 'T145': {'uniprotID': '', 'textname': 'H2A', 'id': 'h2a'}, 'T144': {'uniprotID': '', 'textname': 'Fig. 1a', 'id': 'fig.1a'}, 'T147': {'uniprotID': '', 'textname': 'H2A', 'id': 'h2a'}, 'T146': {'uniprotID': '', 'textname': 'H2B', 'id': 'h2b'}, 'T141': {'uniprotID': '', 'textname': 'H4-tetra-acetylated NCP', 'id': 'h4tetraacetylatedncp'}, 'T140': {'chebi': None, 'textname': 'di-acetylated H4'}, 'T143': {'uniprotID': '', 'textname': 'histones', 'id': 'histone'}, 'T142': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T94': {'uniprotID': '', 'textname': 'human histone H4 ORF', 'id': 'humanhistoneh4orf'}, 'T48': {'chebi': None, 'textname': 'methyllysine'}, 'T95': {'chebi': 'CHEBI:16443', 'textname': 'TAG'}, 'T96': {'uniprotID': '', 'textname': 'T7 RNA polymerase', 'id': 't7rnapolymerase'}, 'T97': {'uniprotID': '', 'textname': 'pyrrolysyl-tRNA synthetase', 'id': 'pyrrolysyltrnasynthetase'}, 'T213': {'uniprotID': '', 'textname': 'NCPs', 'id': 'ncp'}, 'T90': {'uniprotID': u'P62805', 'textname': 'histone H4 proteins', 'uniprotprotnames': u'Histone H4', 'uniprotEntry': u'H4_HUMAN', 'id': u'p62805', 'uniprotgenenames': u'HIST1H4A H4/A H4FA; HIST1H4B H4/I H4FI; HIST1H4C H4/G H4FG; HIST1H4D H4/B H4FB; HIST1H4E H4/J H4FJ; HIST1H4F H4/C H4FC; HIST1H4H H4/H H4FH; HIST1H4I H4/M H4FM; HIST1H4J H4/E H4FE; HIST1H4K H4/D H4FD; HIST1H4L H4/K H4FK; HIST2H4A H4/N H4F2 H4FN HIST2H4; HIST2H4B H4/O H4FO; HIST4H4'}, 'T199': {'chebi': None, 'textname': 'N'}, 'T91': {'chebi': None, 'textname': 'acetyllysine'}, 'T563': {'uniprotID': '', 'textname': 'H4-tetra-acetylated NCP DNA', 'id': 'h4tetraacetylatedncpdna'}, 'T43': {'chebi': 'CHEBI:25094', 'textname': 'lysine'}, 'T561': {'uniprotID': '', 'textname': 'R17', 'id': 'r17'}, 'T560': {'uniprotID': '', 'textname': 'residues K16ac-V21', 'id': 'residuesk16acv21'}, 'T567': {'chebi': 'CHEBI:36976', 'textname': 'nucleotide'}, 'T92': {'uniprotID': '', 'textname': 'Supplementary Fig', 'id': 'supplementaryfig'}, 'T89': {'chebi': None, 'textname': 'acetyllysine'}, 'T88': {'uniprotID': u'P08779', 'textname': 'K16', 'uniprotprotnames': u'Keratin, type I cytoskeletal 16 (Cytokeratin-16) (CK-16) (Keratin-16) (K16)', 'uniprotEntry': u'K1C16_HUMAN', 'id': u'p08779', 'uniprotgenenames': u'KRT16 KRT16A'}, 'T87': {'uniprotID': u'Q8WVN6', 'textname': 'K12', 'uniprotprotnames': u'Secreted and transmembrane protein 1 (Protein K-12)', 'uniprotEntry': u'SCTM1_HUMAN', 'id': u'q8wvn6', 'uniprotgenenames': u'SECTM1 K12'}, 'T86': {'uniprotID': u'P05787', 'textname': 'K8', 'uniprotprotnames': u'Keratin, type II cytoskeletal 8 (Cytokeratin-8) (CK-8) (Keratin-8) (K8) (Type-II keratin Kb8)', 'uniprotEntry': u'K2C8_HUMAN', 'id': u'p05787', 'uniprotgenenames': u'KRT8 CYK8'}, 'T85': {'uniprotID': '', 'textname': 'H4-K5', 'id': 'h4k5'}, 'T84': {'chebi': None, 'textname': 'N'}, 'T83': {'chebi': None, 'textname': 'acetyllysine'}, 'T82': {'uniprotID': u'P62805', 'textname': 'histone H4', 'uniprotprotnames': u'Histone H4', 'uniprotEntry': u'H4_HUMAN', 'id': u'p62805', 'uniprotgenenames': u'HIST1H4A H4/A H4FA; HIST1H4B H4/I H4FI; HIST1H4C H4/G H4FG; HIST1H4D H4/B H4FB; HIST1H4E H4/J H4FJ; HIST1H4F H4/C H4FC; HIST1H4H H4/H H4FH; HIST1H4I H4/M H4FM; HIST1H4J H4/E H4FE; HIST1H4K H4/D H4FD; HIST1H4L H4/K H4FK; HIST2H4A H4/N H4F2 H4FN HIST2H4; HIST2H4B H4/O H4FO; HIST4H4'}, 'T81': {'uniprotID': '', 'textname': 'histone N-terminal tails', 'id': 'histonenterminaltail'}, 'T80': {'chebi': None, 'textname': 'acetyllysine'}, 'T40': {'chebi': 'CHEBI:25094', 'textname': 'lysine'}, 'T192': {'chebi': None, 'textname': 'sugar'}, 'T483': {'uniprotID': '', 'textname': 'E. coli JM109', 'id': 'e.colijm109'}, 'T218': {'uniprotID': u'Q9NQF3', 'textname': 'SHL', 'uniprotprotnames': u'Serine hydrolase-like protein (SHL) (EC 3.1.-.-)', 'uniprotEntry': u'SERHL_HUMAN', 'id': u'q9nqf3', 'uniprotgenenames': u'SERHL SERHL2'}, 'T45': {'chebi': 'CHEBI:29016', 'textname': 'arginine'}, 'T44': {'chebi': 'CHEBI:25094', 'textname': 'lysine'}, 'T363': {'uniprotID': '', 'textname': 'H2B', 'id': 'h2b'}, 'T362': {'uniprotID': '', 'textname': 'molecule-D H2B', 'id': 'moleculedh2b'}, 'T369': {'chebi': None, 'textname': 'dA51'}, 'T368': {'chebi': None, 'textname': 'N'}, 'T38': {'uniprotID': '', 'textname': 'H4 N-terminal tail', 'id': 'h4nterminaltail'}, 'T39': {'uniprotID': '', 'textname': 'chromatin31', 'id': 'chromatin31'}, 'T36': {'chebi': None, 'textname': 'N'}, 'T37': {'chebi': None, 'textname': 'mono-NCP121282930'}, 'T34': {'uniprotID': '', 'textname': 'chromatin fibers2627', 'id': 'chromatinfibers2627'}, 'T35': {'chebi': None, 'textname': 'N'}, 'T32': {'uniprotID': '', 'textname': 'chromatin2122232425', 'id': 'chromatin2122232425'}, 'T33': {'chebi': None, 'textname': 'N'}, 'T30': {'uniprotID': '', 'textname': 'histone DNA', 'id': 'histonedna'}, 'T31': {'chebi': None, 'textname': 'N'}, 'T329': {'uniprotID': u'P08779', 'textname': 'K16', 'uniprotprotnames': u'Keratin, type I cytoskeletal 16 (Cytokeratin-16) (CK-16) (Keratin-16) (K16)', 'uniprotEntry': u'K1C16_HUMAN', 'id': u'p08779', 'uniprotgenenames': u'KRT16 KRT16A'}, 'T328': {'uniprotID': u'Q8WVN6', 'textname': 'K12', 'uniprotprotnames': u'Secreted and transmembrane protein 1 (Protein K-12)', 'uniprotEntry': u'SCTM1_HUMAN', 'id': u'q8wvn6', 'uniprotgenenames': u'SECTM1 K12'}, 'T537': {'uniprotID': '', 'textname': 'TAF1', 'id': 'taf1'}, 'T321': {'uniprotID': '', 'textname': 'TAF1-BrD', 'id': 'taf1brd'}, 'T320': {'uniprotID': '', 'textname': 'acetylated tail peptides52', 'id': 'acetylatedtailpeptides52'}, 'T323': {'uniprotID': '', 'textname': 'S1', 'id': 's1'}, 'T322': {'uniprotID': '', 'textname': 'human NCPs', 'id': 'humanncp'}, 'T325': {'chebi': None, 'textname': 'N'}, 'T324': {'chebi': None, 'textname': 'acetyllysine'}, 'T327': {'uniprotID': u'P05787', 'textname': 'K8', 'uniprotprotnames': u'Keratin, type II cytoskeletal 8 (Cytokeratin-8) (CK-8) (Keratin-8) (K8) (Type-II keratin Kb8)', 'uniprotEntry': u'K2C8_HUMAN', 'id': u'p05787', 'uniprotgenenames': u'KRT8 CYK8'}, 'T326': {'uniprotID': '', 'textname': 'K5', 'id': 'k5'}, 'T518': {'chebi': 'CHEBI:32588', 'textname': 'potassium chloride'}, 'T519': {'chebi': None, 'textname': 'trehalose'}, 'T512': {'chebi': 'CHEBI:6636', 'textname': 'MgCl2'}, 'T513': {'chebi': None, 'textname': 'potassium cacodylate'}, 'T510': {'uniprotID': '', 'textname': 'NCPs', 'id': 'ncp'}, 'T511': {'chebi': 'CHEBI:6636', 'textname': 'MgCl2'}, 'T516': {'chebi': None, 'textname': 'potassium cacodylate'}, 'T517': {'chebi': 'CHEBI:63041', 'textname': 'manganese(II) chloride'}, 'T514': {'chebi': 'CHEBI:63041', 'textname': 'manganese(II) chloride'}, 'T515': {'chebi': 'CHEBI:32588', 'textname': 'potassium chloride'}, 'T215': {'uniprotID': '', 'textname': 'NCPs', 'id': 'ncp'}, 'T214': {'chebi': 'CHEBI:49637', 'textname': 'hydrogen'}, 'T217': {'chebi': 'CHEBI:36976', 'textname': 'nucleotide'}, 'T216': {'uniprotID': '', 'textname': 'NCPs', 'id': 'ncp'}, 'T211': {'uniprotID': '', 'textname': 'E110', 'id': 'e110'}, 'T210': {'uniprotID': '', 'textname': 'H2A', 'id': 'h2a'}, 'T198': {'uniprotID': u'Q9HD26', 'textname': 'Fig', 'uniprotprotnames': u'Golgi-associated PDZ and coiled-coil motif-containing protein (CFTR-associated ligand) (Fused in glioblastoma) (PDZ protein interacting specifically with TC10) (PIST)', 'uniprotEntry': u'GOPC_HUMAN', 'id': u'q9hd26', 'uniprotgenenames': u'GOPC CAL FIG'}, 'T212': {'uniprotID': '', 'textname': 'H2B', 'id': 'h2b'}, 'T196': {'chebi': None, 'textname': 'phosphate'}, 'T197': {'uniprotID': '', 'textname': 'tetra-acetylated H4 tail', 'id': 'tetraacetylatedh4tail'}, 'T194': {'uniprotID': '', 'textname': 'Rb+-coordinating NCP147', 'id': 'rb+coordinatingncp147'}, 'T195': {'chebi': None, 'textname': 'sugar'}, 'T219': {'uniprotID': '', 'textname': 'NCPs1', 'id': 'ncps1'}, 'T193': {'chebi': None, 'textname': 'phosphate'}, 'T190': {'uniprotID': '', 'textname': 'Fig. 2a', 'id': 'fig.2a'}, 'T191': {'uniprotID': '', 'textname': 'Fig. 2c', 'id': 'fig.2c'}, 'T358': {'uniprotID': '', 'textname': 'DNA1', 'id': 'dna1'}, 'T359': {'chebi': None, 'textname': 'N'}, 'T354': {'uniprotID': u'P08779', 'textname': 'K16', 'uniprotprotnames': u'Keratin, type I cytoskeletal 16 (Cytokeratin-16) (CK-16) (Keratin-16) (K16)', 'uniprotEntry': u'K1C16_HUMAN', 'id': u'p08779', 'uniprotgenenames': u'KRT16 KRT16A'}, 'T355': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T356': {'chebi': None, 'textname': 'N'}, 'T357': {'uniprotID': '', 'textname': 'H4 N-terminal tail', 'id': 'h4nterminaltail'}, 'T350': {'uniprotID': u'P05787', 'textname': 'K8', 'uniprotprotnames': u'Keratin, type II cytoskeletal 8 (Cytokeratin-8) (CK-8) (Keratin-8) (K8) (Type-II keratin Kb8)', 'uniprotEntry': u'K2C8_HUMAN', 'id': u'p05787', 'uniprotgenenames': u'KRT8 CYK8'}, 'T351': {'uniprotID': u'Q8WVN6', 'textname': 'K12', 'uniprotprotnames': u'Secreted and transmembrane protein 1 (Protein K-12)', 'uniprotEntry': u'SCTM1_HUMAN', 'id': u'q8wvn6', 'uniprotgenenames': u'SECTM1 K12'}, 'T352': {'uniprotID': u'P08779', 'textname': 'K16', 'uniprotprotnames': u'Keratin, type I cytoskeletal 16 (Cytokeratin-16) (CK-16) (Keratin-16) (K16)', 'uniprotEntry': u'K1C16_HUMAN', 'id': u'p08779', 'uniprotgenenames': u'KRT16 KRT16A'}, 'T353': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T174': {'uniprotID': '', 'textname': 'H4 protein', 'id': 'h4'}, 'T175': {'uniprotID': '', 'textname': 'H4 proteins', 'id': 'h4'}, 'T176': {'chebi': None, 'textname': 'N'}, 'T177': {'uniprotID': '', 'textname': 'histone proteins', 'id': 'histone'}, 'T170': {'uniprotID': u'P05787', 'textname': 'K8', 'uniprotprotnames': u'Keratin, type II cytoskeletal 8 (Cytokeratin-8) (CK-8) (Keratin-8) (K8) (Type-II keratin Kb8)', 'uniprotEntry': u'K2C8_HUMAN', 'id': u'p05787', 'uniprotgenenames': u'KRT8 CYK8'}, 'T171': {'uniprotID': u'Q8WVN6', 'textname': 'K12', 'uniprotprotnames': u'Secreted and transmembrane protein 1 (Protein K-12)', 'uniprotEntry': u'SCTM1_HUMAN', 'id': u'q8wvn6', 'uniprotgenenames': u'SECTM1 K12'}, 'T172': {'uniprotID': u'P08779', 'textname': 'K16', 'uniprotprotnames': u'Keratin, type I cytoskeletal 16 (Cytokeratin-16) (CK-16) (Keratin-16) (K16)', 'uniprotEntry': u'K1C16_HUMAN', 'id': u'p08779', 'uniprotgenenames': u'KRT16 KRT16A'}, 'T173': {'uniprotID': '', 'textname': 'NCPs', 'id': 'ncp'}, 'T178': {'uniprotID': '', 'textname': 'Supplementary Fig', 'id': 'supplementaryfig'}, 'T179': {'chebi': 'CHEBI:8984', 'textname': 'SDS'}, 'T589': {'chebi': None, 'textname': 'Mn++'}, 'T588': {'chebi': None, 'textname': 'D77'}, 'T585': {'chebi': 'CHEBI:36976', 'textname': 'nucleotide'}, 'T584': {'uniprotID': '', 'textname': 'BrD', 'id': 'brd'}, 'T587': {'chebi': None, 'textname': 'V48'}, 'T586': {'chebi': None, 'textname': 'carbonyl O'}, 'T581': {'uniprotID': '', 'textname': 'TAF1', 'id': 'taf1'}, 'T580': {'uniprotID': '', 'textname': 'TAF1-BrD', 'id': 'taf1brd'}, 'T583': {'uniprotID': '', 'textname': 'TAF1', 'id': 'taf1'}, 'T582': {'uniprotID': '', 'textname': 'BrD', 'id': 'brd'}, 'T101': {'chebi': None, 'textname': 'Kac'}, 'T100': {'chebi': 'CHEBI:21860', 'textname': 'pyrrolysine'}, 'T103': {'uniprotID': u'P62805', 'textname': 'histone H4 protein', 'uniprotprotnames': u'Histone H4', 'uniprotEntry': u'H4_HUMAN', 'id': u'p62805', 'uniprotgenenames': u'HIST1H4A H4/A H4FA; HIST1H4B H4/I H4FI; HIST1H4C H4/G H4FG; HIST1H4D H4/B H4FB; HIST1H4E H4/J H4FJ; HIST1H4F H4/C H4FC; HIST1H4H H4/H H4FH; HIST1H4I H4/M H4FM; HIST1H4J H4/E H4FE; HIST1H4K H4/D H4FD; HIST1H4L H4/K H4FK; HIST2H4A H4/N H4F2 H4FN HIST2H4; HIST2H4B H4/O H4FO; HIST4H4'}, 'T102': {'uniprotID': '', 'textname': 'resultant Kac-tRNAPyl', 'id': 'resultantkactrnapyl'}, 'T105': {'uniprotID': u'P62805', 'textname': 'histone H4 proteins', 'uniprotprotnames': u'Histone H4', 'uniprotEntry': u'H4_HUMAN', 'id': u'p62805', 'uniprotgenenames': u'HIST1H4A H4/A H4FA; HIST1H4B H4/I H4FI; HIST1H4C H4/G H4FG; HIST1H4D H4/B H4FB; HIST1H4E H4/J H4FJ; HIST1H4F H4/C H4FC; HIST1H4H H4/H H4FH; HIST1H4I H4/M H4FM; HIST1H4J H4/E H4FE; HIST1H4K H4/D H4FD; HIST1H4L H4/K H4FK; HIST2H4A H4/N H4F2 H4FN HIST2H4; HIST2H4B H4/O H4FO; HIST4H4'}, 'T104': {'chebi': None, 'textname': 'acetyllysine'}, 'T107': {'uniprotID': '', 'textname': 'S1b', 'id': 's1b'}, 'T106': {'uniprotID': u'Q9HD26', 'textname': 'Fig', 'uniprotprotnames': u'Golgi-associated PDZ and coiled-coil motif-containing protein (CFTR-associated ligand) (Fused in glioblastoma) (PDZ protein interacting specifically with TC10) (PIST)', 'uniprotEntry': u'GOPC_HUMAN', 'id': u'q9hd26', 'uniprotgenenames': u'GOPC CAL FIG'}, 'T109': {'uniprotID': '', 'textname': 'system lacked release factor 1', 'id': 'systemlackedreleasefactor1'}, 'T108': {'chebi': 'CHEBI:8984', 'textname': 'SDS'}, 'T459': {'uniprotID': '', 'textname': 'pyrrolysyl aminoacyl tRNA synthetase', 'id': 'pyrrolysylaminoacyltrnasynthetase'}, 'T458': {'uniprotID': '', 'textname': 'H4 proteins', 'id': 'h4'}, 'T415': {'uniprotID': u'P62805', 'textname': 'histone H4', 'uniprotprotnames': u'Histone H4', 'uniprotEntry': u'H4_HUMAN', 'id': u'p62805', 'uniprotgenenames': u'HIST1H4A H4/A H4FA; HIST1H4B H4/I H4FI; HIST1H4C H4/G H4FG; HIST1H4D H4/B H4FB; HIST1H4E H4/J H4FJ; HIST1H4F H4/C H4FC; HIST1H4H H4/H H4FH; HIST1H4I H4/M H4FM; HIST1H4J H4/E H4FE; HIST1H4K H4/D H4FD; HIST1H4L H4/K H4FK; HIST2H4A H4/N H4F2 H4FN HIST2H4; HIST2H4B H4/O H4FO; HIST4H4'}, 'T414': {'uniprotID': '', 'textname': 'TAF1-BrD', 'id': 'taf1brd'}, 'T413': {'uniprotID': '', 'textname': 'TAF1-BrD', 'id': 'taf1brd'}, 'T412': {'uniprotID': '', 'textname': 'TAF1', 'id': 'taf1'}, 'T181': {'uniprotID': '', 'textname': 'H4-tetra-acetylated H4 protein', 'id': 'h4tetraacetylatedh4'}, 'T260': {'uniprotID': u'Q9HD26', 'textname': 'Fig', 'uniprotprotnames': u'Golgi-associated PDZ and coiled-coil motif-containing protein (CFTR-associated ligand) (Fused in glioblastoma) (PDZ protein interacting specifically with TC10) (PIST)', 'uniprotEntry': u'GOPC_HUMAN', 'id': u'q9hd26', 'uniprotgenenames': u'GOPC CAL FIG'}, 'T411': {'uniprotID': '', 'textname': 'acetyllysine-binding double bromodomain', 'id': 'acetyllysinebindingdoublebromodomain'}, 'T410': {'uniprotID': '', 'textname': 'histone histone', 'id': 'histonehistone'}, 'T556': {'uniprotID': '', 'textname': 'molecule-B H4', 'id': 'moleculebh4'}, 'T180': {'uniprotID': '', 'textname': 'S2a', 'id': 's2a'}, 'T502': {'uniprotID': '', 'textname': 'NCPs', 'id': 'ncp'}, 'T557': {'uniprotID': '', 'textname': 'SHLs', 'id': 'shl'}, 'T453': {'uniprotID': '', 'textname': 'MKDHLIHNHHKHEHAHALVPRGSHM', 'id': 'mkdhlihnhhkhehahalvprgshm'}, 'T452': {'uniprotID': '', 'textname': 'TEV-protease recognition sequence', 'id': 'tevproteaserecognitionsequence'}, 'T451': {'uniprotID': '', 'textname': 'thrombin', 'id': 'thrombin'}, 'T527': {'uniprotID': '', 'textname': 'NCPs', 'id': 'ncp'}, 'T450': {'uniprotID': u'P62805', 'textname': 'histone H4 protein', 'uniprotprotnames': u'Histone H4', 'uniprotEntry': u'H4_HUMAN', 'id': u'p62805', 'uniprotgenenames': u'HIST1H4A H4/A H4FA; HIST1H4B H4/I H4FI; HIST1H4C H4/G H4FG; HIST1H4D H4/B H4FB; HIST1H4E H4/J H4FJ; HIST1H4F H4/C H4FC; HIST1H4H H4/H H4FH; HIST1H4I H4/M H4FM; HIST1H4J H4/E H4FE; HIST1H4K H4/D H4FD; HIST1H4L H4/K H4FK; HIST2H4A H4/N H4F2 H4FN HIST2H4; HIST2H4B H4/O H4FO; HIST4H4'}, 'T526': {'uniprotID': '', 'textname': 'Coot59', 'id': 'coot59'}, 'T457': {'chebi': 'CHEBI:27570', 'textname': 'histidine'}, 'T525': {'uniprotID': '', 'textname': 'PHENIX suite48', 'id': 'phenixsuite48'}, 'T456': {'chebi': None, 'textname': 'N'}, 'T389': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T388': {'chebi': None, 'textname': 'H2B'}, 'T455': {'uniprotID': '', 'textname': 'thrombin', 'id': 'thrombin'}, 'T383': {'chebi': None, 'textname': 'N'}, 'T382': {'uniprotID': '', 'textname': 'S3', 'id': 's3'}, 'T381': {'uniprotID': u'Q9HD26', 'textname': 'Fig', 'uniprotprotnames': u'Golgi-associated PDZ and coiled-coil motif-containing protein (CFTR-associated ligand) (Fused in glioblastoma) (PDZ protein interacting specifically with TC10) (PIST)', 'uniprotEntry': u'GOPC_HUMAN', 'id': u'q9hd26', 'uniprotgenenames': u'GOPC CAL FIG'}, 'T380': {'uniprotID': '', 'textname': 'histone histone', 'id': 'histonehistone'}, 'T387': {'chebi': 'CHEBI:23019', 'textname': 'carbonyl'}, 'T386': {'chebi': 'CHEBI:49637', 'textname': 'hydrogen'}, 'T385': {'uniprotID': u'P08779', 'textname': 'K16', 'uniprotprotnames': u'Keratin, type I cytoskeletal 16 (Cytokeratin-16) (CK-16) (Keratin-16) (K16)', 'uniprotEntry': u'K1C16_HUMAN', 'id': u'p08779', 'uniprotgenenames': u'KRT16 KRT16A'}, 'T384': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T521': {'chebi': 'CHEBI:25555', 'textname': 'nitrogen'}, 'T520': {'chebi': None, 'textname': '2-methyl-2,4-pentanediol1'}, 'T310': {'uniprotID': u'Q9HD26', 'textname': 'Fig', 'uniprotprotnames': u'Golgi-associated PDZ and coiled-coil motif-containing protein (CFTR-associated ligand) (Fused in glioblastoma) (PDZ protein interacting specifically with TC10) (PIST)', 'uniprotEntry': u'GOPC_HUMAN', 'id': u'q9hd26', 'uniprotgenenames': u'GOPC CAL FIG'}, 'T311': {'uniprotID': '', 'textname': 'TAF1', 'id': 'taf1'}, 'T312': {'uniprotID': '', 'textname': 'BrD', 'id': 'brd'}, 'T313': {'uniprotID': '', 'textname': 'K5', 'id': 'k5'}, 'T314': {'uniprotID': u'P05787', 'textname': 'K8', 'uniprotprotnames': u'Keratin, type II cytoskeletal 8 (Cytokeratin-8) (CK-8) (Keratin-8) (K8) (Type-II keratin Kb8)', 'uniprotEntry': u'K2C8_HUMAN', 'id': u'p05787', 'uniprotgenenames': u'KRT8 CYK8'}, 'T315': {'uniprotID': u'Q8WVN6', 'textname': 'K12', 'uniprotprotnames': u'Secreted and transmembrane protein 1 (Protein K-12)', 'uniprotEntry': u'SCTM1_HUMAN', 'id': u'q8wvn6', 'uniprotgenenames': u'SECTM1 K12'}, 'T316': {'uniprotID': u'P08779', 'textname': 'K16', 'uniprotprotnames': u'Keratin, type I cytoskeletal 16 (Cytokeratin-16) (CK-16) (Keratin-16) (K16)', 'uniprotEntry': u'K1C16_HUMAN', 'id': u'p08779', 'uniprotgenenames': u'KRT16 KRT16A'}, 'T317': {'chebi': None, 'textname': 'K8-di-acetylated NCP'}, 'T251': {'uniprotID': u'Q9HD26', 'textname': 'Fig', 'uniprotprotnames': u'Golgi-associated PDZ and coiled-coil motif-containing protein (CFTR-associated ligand) (Fused in glioblastoma) (PDZ protein interacting specifically with TC10) (PIST)', 'uniprotEntry': u'GOPC_HUMAN', 'id': u'q9hd26', 'uniprotgenenames': u'GOPC CAL FIG'}, 'T250': {'uniprotID': '', 'textname': 'H2B', 'id': 'h2b'}, 'T253': {'uniprotID': u'Q9HD26', 'textname': 'Fig', 'uniprotprotnames': u'Golgi-associated PDZ and coiled-coil motif-containing protein (CFTR-associated ligand) (Fused in glioblastoma) (PDZ protein interacting specifically with TC10) (PIST)', 'uniprotEntry': u'GOPC_HUMAN', 'id': u'q9hd26', 'uniprotgenenames': u'GOPC CAL FIG'}, 'T252': {'uniprotID': '', 'textname': 'S3b', 'id': 's3b'}, 'T255': {'uniprotID': '', 'textname': 'NCPs', 'id': 'ncp'}, 'T254': {'uniprotID': '', 'textname': 'S3c', 'id': 's3c'}, 'T257': {'uniprotID': '', 'textname': 'S3c', 'id': 's3c'}, 'T256': {'uniprotID': '', 'textname': 'molecule-A H3', 'id': 'moleculeah3'}, 'T138': {'uniprotID': '', 'textname': 'S1f', 'id': 's1f'}, 'T139': {'chebi': 'CHEBI:6636', 'textname': 'MgCl2'}, 'T428': {'uniprotID': '', 'textname': 'TAF1', 'id': 'taf1'}, 'T429': {'uniprotID': '', 'textname': 'BrD', 'id': 'brd'}, 'T130': {'uniprotID': '', 'textname': 'NCPs', 'id': 'ncp'}, 'T131': {'uniprotID': '', 'textname': 'S1e', 'id': 's1e'}, 'T132': {'uniprotID': '', 'textname': 'K5', 'id': 'k5'}, 'T133': {'uniprotID': u'P05787', 'textname': 'K8', 'uniprotprotnames': u'Keratin, type II cytoskeletal 8 (Cytokeratin-8) (CK-8) (Keratin-8) (K8) (Type-II keratin Kb8)', 'uniprotEntry': u'K2C8_HUMAN', 'id': u'p05787', 'uniprotgenenames': u'KRT8 CYK8'}, 'T134': {'uniprotID': u'Q8WVN6', 'textname': 'K12', 'uniprotprotnames': u'Secreted and transmembrane protein 1 (Protein K-12)', 'uniprotEntry': u'SCTM1_HUMAN', 'id': u'q8wvn6', 'uniprotgenenames': u'SECTM1 K12'}, 'T135': {'uniprotID': u'P08779', 'textname': 'K16', 'uniprotprotnames': u'Keratin, type I cytoskeletal 16 (Cytokeratin-16) (CK-16) (Keratin-16) (K16)', 'uniprotEntry': u'K1C16_HUMAN', 'id': u'p08779', 'uniprotgenenames': u'KRT16 KRT16A'}, 'T136': {'uniprotID': '', 'textname': 'histone H4 (H4-tetra-acetylated NCP)', 'id': 'histoneh4h4tetraacetylatedncp'}, 'T137': {'uniprotID': u'Q9HD26', 'textname': 'Fig', 'uniprotprotnames': u'Golgi-associated PDZ and coiled-coil motif-containing protein (CFTR-associated ligand) (Fused in glioblastoma) (PDZ protein interacting specifically with TC10) (PIST)', 'uniprotEntry': u'GOPC_HUMAN', 'id': u'q9hd26', 'uniprotgenenames': u'GOPC CAL FIG'}, 'T550': {'chebi': None, 'textname': 'Mg++'}, 'T562': {'uniprotID': '', 'textname': 'V21', 'id': 'v21'}, 'T228': {'chebi': None, 'textname': 'N'}, 'T229': {'chebi': None, 'textname': 'N'}, 'T558': {'chebi': 'CHEBI:25555', 'textname': 'N'}, 'T559': {'chebi': 'CHEBI:25555', 'textname': 'N'}, 'T224': {'uniprotID': '', 'textname': 'Fig. 3a', 'id': 'fig.3a'}, 'T225': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T226': {'uniprotID': u'Q9HD26', 'textname': 'Fig', 'uniprotprotnames': u'Golgi-associated PDZ and coiled-coil motif-containing protein (CFTR-associated ligand) (Fused in glioblastoma) (PDZ protein interacting specifically with TC10) (PIST)', 'uniprotEntry': u'GOPC_HUMAN', 'id': u'q9hd26', 'uniprotgenenames': u'GOPC CAL FIG'}, 'T227': {'chebi': None, 'textname': 'N'}, 'T220': {'uniprotID': '', 'textname': 'superhelical location', 'id': 'superhelicallocation'}, 'T221': {'uniprotID': '', 'textname': 'NCPs', 'id': 'ncp'}, 'T222': {'chebi': 'CHEBI:36976', 'textname': 'nucleotides'}, 'T93': {'uniprotID': '', 'textname': 'S1a', 'id': 's1a'}, 'T565': {'uniprotID': '', 'textname': 'H4-tetra-acetylated NCPs', 'id': 'h4tetraacetylatedncp'}, 'T554': {'uniprotID': u'P08779', 'textname': 'K16', 'uniprotprotnames': u'Keratin, type I cytoskeletal 16 (Cytokeratin-16) (CK-16) (Keratin-16) (K16)', 'uniprotEntry': u'K1C16_HUMAN', 'id': u'p08779', 'uniprotgenenames': u'KRT16 KRT16A'}, 'T564': {'chebi': 'CHEBI:36976', 'textname': 'nucleotide'}, 'T49': {'uniprotID': '', 'textname': 'H4 N-terminal tail', 'id': 'h4nterminaltail'}, 'T440': {'chebi': None, 'textname': 'N'}, 'T569': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T568': {'chebi': 'CHEBI:25555', 'textname': 'N'}, 'T454': {'uniprotID': '', 'textname': 'N11', 'id': 'n11'}, 'T542': {'uniprotID': u'P62805', 'textname': 'histone H4', 'uniprotprotnames': u'Histone H4', 'uniprotEntry': u'H4_HUMAN', 'id': u'p62805', 'uniprotgenenames': u'HIST1H4A H4/A H4FA; HIST1H4B H4/I H4FI; HIST1H4C H4/G H4FG; HIST1H4D H4/B H4FB; HIST1H4E H4/J H4FJ; HIST1H4F H4/C H4FC; HIST1H4H H4/H H4FH; HIST1H4I H4/M H4FM; HIST1H4J H4/E H4FE; HIST1H4K H4/D H4FD; HIST1H4L H4/K H4FK; HIST2H4A H4/N H4F2 H4FN HIST2H4; HIST2H4B H4/O H4FO; HIST4H4'}, 'T566': {'chebi': 'CHEBI:36976', 'textname': 'nucleotide'}, 'T524': {'uniprotID': '', 'textname': 'histone octamer', 'id': 'histoneoctamer'}, 'T69': {'chebi': 'CHEBI:25094', 'textname': 'lysine'}, 'T555': {'uniprotID': '', 'textname': 'NCPs', 'id': 'ncp'}, 'T68': {'chebi': None, 'textname': 'Kcme3'}, 'T67': {'chebi': None, 'textname': 'N-trimethyl-aminoethylcysteine'}, 'T66': {'chebi': None, 'textname': 'Kcme2'}, 'T536': {'chebi': 'CHEBI:42106', 'textname': 'DTT'}, 'T42': {'chebi': None, 'textname': 'N'}, 'T189': {'uniprotID': '', 'textname': 'NCPs', 'id': 'ncp'}, 'T188': {'uniprotID': '', 'textname': 'H4-tetra-acetylated NCP', 'id': 'h4tetraacetylatedncp'}, 'T262': {'chebi': None, 'textname': 'N'}, 'T263': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T264': {'uniprotID': '', 'textname': 'H2B heterodimer1', 'id': 'h2bheterodimer1'}, 'T265': {'uniprotID': '', 'textname': 'NCPs', 'id': 'ncp'}, 'T266': {'uniprotID': '', 'textname': 'H2A', 'id': 'h2a'}, 'T267': {'uniprotID': '', 'textname': 'D90', 'id': 'd90'}, 'T268': {'uniprotID': u'Q9HD26', 'textname': 'Fig', 'uniprotprotnames': u'Golgi-associated PDZ and coiled-coil motif-containing protein (CFTR-associated ligand) (Fused in glioblastoma) (PDZ protein interacting specifically with TC10) (PIST)', 'uniprotEntry': u'GOPC_HUMAN', 'id': u'q9hd26', 'uniprotgenenames': u'GOPC CAL FIG'}, 'T51': {'uniprotID': u'P08779', 'textname': 'K16', 'uniprotprotnames': u'Keratin, type I cytoskeletal 16 (Cytokeratin-16) (CK-16) (Keratin-16) (K16)', 'uniprotEntry': u'K1C16_HUMAN', 'id': u'p08779', 'uniprotgenenames': u'KRT16 KRT16A'}, 'T52': {'uniprotID': '', 'textname': 'H4-K16Q', 'id': 'h4k16q'}, 'T182': {'uniprotID': '', 'textname': 'H4 protein', 'id': 'h4'}, 'T185': {'uniprotID': '', 'textname': 'H4 protein', 'id': 'h4'}, 'T184': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T187': {'uniprotID': '', 'textname': 'acetyllysine recognition antibodies', 'id': 'acetyllysinerecognitionantibodie'}, 'T186': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T349': {'uniprotID': '', 'textname': 'H4-K5', 'id': 'h4k5'}, 'T348': {'uniprotID': u'P08779', 'textname': 'K16', 'uniprotprotnames': u'Keratin, type I cytoskeletal 16 (Cytokeratin-16) (CK-16) (Keratin-16) (K16)', 'uniprotEntry': u'K1C16_HUMAN', 'id': u'p08779', 'uniprotgenenames': u'KRT16 KRT16A'}, 'T347': {'uniprotID': u'Q8WVN6', 'textname': 'K12', 'uniprotprotnames': u'Secreted and transmembrane protein 1 (Protein K-12)', 'uniprotEntry': u'SCTM1_HUMAN', 'id': u'q8wvn6', 'uniprotgenenames': u'SECTM1 K12'}, 'T346': {'uniprotID': u'P05787', 'textname': 'K8', 'uniprotprotnames': u'Keratin, type II cytoskeletal 8 (Cytokeratin-8) (CK-8) (Keratin-8) (K8) (Type-II keratin Kb8)', 'uniprotEntry': u'K2C8_HUMAN', 'id': u'p05787', 'uniprotgenenames': u'KRT8 CYK8'}, 'T345': {'uniprotID': '', 'textname': 'H4-K5', 'id': 'h4k5'}, 'T344': {'uniprotID': u'P08779', 'textname': 'K16', 'uniprotprotnames': u'Keratin, type I cytoskeletal 16 (Cytokeratin-16) (CK-16) (Keratin-16) (K16)', 'uniprotEntry': u'K1C16_HUMAN', 'id': u'p08779', 'uniprotgenenames': u'KRT16 KRT16A'}, 'T343': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T41': {'chebi': 'CHEBI:29016', 'textname': 'arginine'}, 'T341': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T340': {'uniprotID': u'P08779', 'textname': 'K16', 'uniprotprotnames': u'Keratin, type I cytoskeletal 16 (Cytokeratin-16) (CK-16) (Keratin-16) (K16)', 'uniprotEntry': u'K1C16_HUMAN', 'id': u'p08779', 'uniprotgenenames': u'KRT16 KRT16A'}, 'T167': {'uniprotID': '', 'textname': 'thrombin', 'id': 'thrombin'}, 'T166': {'chebi': 'CHEBI:27570', 'textname': 'histidine'}, 'T165': {'chebi': None, 'textname': 'N'}, 'T164': {'uniprotID': '', 'textname': 'thrombin', 'id': 'thrombin'}, 'T163': {'uniprotID': '', 'textname': 'H4 protein', 'id': 'h4'}, 'T162': {'uniprotID': '', 'textname': 'H4-tetra-acetylated NCP', 'id': 'h4tetraacetylatedncp'}, 'T161': {'chebi': None, 'textname': 'N'}, 'T160': {'uniprotID': u'P08779', 'textname': 'K16', 'uniprotprotnames': u'Keratin, type I cytoskeletal 16 (Cytokeratin-16) (CK-16) (Keratin-16) (K16)', 'uniprotEntry': u'K1C16_HUMAN', 'id': u'p08779', 'uniprotgenenames': u'KRT16 KRT16A'}, 'T169': {'uniprotID': '', 'textname': 'H4 protein tetra-acetylated at K5', 'id': 'h4tetraacetylatedatk5'}, 'T168': {'uniprotID': '', 'textname': 'TEV protease', 'id': 'tevprotease'}, 'T530': {'uniprotID': '', 'textname': 'NCPs', 'id': 'ncp'}, 'T523': {'uniprotID': '', 'textname': 'HKL200058', 'id': 'hkl200058'}, 'T377': {'uniprotID': u'P05787', 'textname': 'K8', 'uniprotprotnames': u'Keratin, type II cytoskeletal 8 (Cytokeratin-8) (CK-8) (Keratin-8) (K8) (Type-II keratin Kb8)', 'uniprotEntry': u'K2C8_HUMAN', 'id': u'p05787', 'uniprotgenenames': u'KRT8 CYK8'}, 'T505': {'chebi': 'CHEBI:64755', 'textname': 'EDTA'}, 'T504': {'chebi': None, 'textname': 'Tris-HCl'}, 'T507': {'chebi': 'CHEBI:6636', 'textname': 'MgCl2'}, 'T506': {'chebi': None, 'textname': 'Tris-HCl'}, 'T501': {'chebi': None, 'textname': 'Mg++'}, 'T500': {'chebi': 'CHEBI:42106', 'textname': 'dithiothreitol'}, 'T468': {'chebi': 'CHEBI:17883', 'textname': 'HCl'}, 'T469': {'chebi': 'CHEBI:26710', 'textname': 'NaCl'}, 'T466': {'chebi': None, 'textname': 'Tris-HCl'}, 'T47': {'chebi': None, 'textname': 'acetyllysine'}, 'T464': {'uniprotID': '', 'textname': 'S30 extract', 'id': 's30extract'}, 'T465': {'uniprotID': '', 'textname': 'H4 proteins', 'id': 'h4'}, 'T462': {'uniprotID': '', 'textname': 'KacRS_6mt', 'id': 'kacrs_6mt'}, 'T463': {'uniprotID': '', 'textname': 'UAG-recognizing tRNAPyl', 'id': 'uagrecognizingtrnapyl'}, 'T460': {'uniprotID': '', 'textname': 'PylRS', 'id': 'pylr'}, 'T461': {'uniprotID': '', 'textname': 'C348F', 'id': 'c348f'}, 'T484': {'uniprotID': '', 'textname': 'DE3', 'id': 'de3'}, 'T485': {'uniprotID': '', 'textname': 'histones H2A type 1-B/E', 'id': 'histonesh2atype1b/e'}, 'T486': {'uniprotID': '', 'textname': 'H2B type 1-J', 'id': 'h2btype1j'}, 'T487': {'uniprotID': '', 'textname': 'H2A', 'id': 'h2a'}, 'T480': {'uniprotID': '', 'textname': 'TEV protease', 'id': 'tevprotease'}, 'T481': {'uniprotID': '', 'textname': 'Histone proteins', 'id': 'histone'}, 'T482': {'uniprotID': u'P62805', 'textname': 'histone H4 protein', 'uniprotprotnames': u'Histone H4', 'uniprotEntry': u'H4_HUMAN', 'id': u'p62805', 'uniprotgenenames': u'HIST1H4A H4/A H4FA; HIST1H4B H4/I H4FI; HIST1H4C H4/G H4FG; HIST1H4D H4/B H4FB; HIST1H4E H4/J H4FJ; HIST1H4F H4/C H4FC; HIST1H4H H4/H H4FH; HIST1H4I H4/M H4FM; HIST1H4J H4/E H4FE; HIST1H4K H4/D H4FD; HIST1H4L H4/K H4FK; HIST2H4A H4/N H4F2 H4FN HIST2H4; HIST2H4B H4/O H4FO; HIST4H4'}, 'T46': {'chebi': None, 'textname': 'N'}, 'T408': {'chebi': 'CHEBI:25094', 'textname': 'lysine'}, 'T488': {'uniprotID': '', 'textname': 'H2B', 'id': 'h2b'}, 'T489': {'uniprotID': '', 'textname': 'NCPs', 'id': 'ncp'}, 'T409': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T544': {'chebi': None, 'textname': 'acetyllysine'}, 'T405': {'uniprotID': '', 'textname': 'S1', 'id': 's1'}, 'T406': {'uniprotID': u'P08779', 'textname': 'K16', 'uniprotprotnames': u'Keratin, type I cytoskeletal 16 (Cytokeratin-16) (CK-16) (Keratin-16) (K16)', 'uniprotEntry': u'K1C16_HUMAN', 'id': u'p08779', 'uniprotgenenames': u'KRT16 KRT16A'}, 'T407': {'chebi': None, 'textname': 'N'}, 'T403': {'uniprotID': '', 'textname': 'H2B', 'id': 'h2b'}, 'T261': {'uniprotID': '', 'textname': 'S3c', 'id': 's3c'}, 'T552': {'uniprotID': u'P05787', 'textname': 'K8', 'uniprotprotnames': u'Keratin, type II cytoskeletal 8 (Cytokeratin-8) (CK-8) (Keratin-8) (K8) (Type-II keratin Kb8)', 'uniprotEntry': u'K2C8_HUMAN', 'id': u'p05787', 'uniprotgenenames': u'KRT8 CYK8'}, 'T522': {'uniprotID': '', 'textname': 'XDS57', 'id': 'xds57'}, 'T553': {'uniprotID': u'Q8WVN6', 'textname': 'K12', 'uniprotprotnames': u'Secreted and transmembrane protein 1 (Protein K-12)', 'uniprotEntry': u'SCTM1_HUMAN', 'id': u'q8wvn6', 'uniprotgenenames': u'SECTM1 K12'}, 'T129': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T14': {'uniprotID': '', 'textname': 'H2B', 'id': 'h2b'}, 'T15': {'uniprotID': '', 'textname': 'H3', 'id': 'h3'}, 'T16': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T17': {'uniprotID': '', 'textname': 'core histone', 'id': 'corehistone'}, 'T10': {'uniprotID': u'P08779', 'textname': 'K16', 'uniprotprotnames': u'Keratin, type I cytoskeletal 16 (Cytokeratin-16) (CK-16) (Keratin-16) (K16)', 'uniprotEntry': u'K1C16_HUMAN', 'id': u'p08779', 'uniprotgenenames': u'KRT16 KRT16A'}, 'T11': {'uniprotID': '', 'textname': 'histone octamer', 'id': 'histoneoctamer'}, 'T12': {'chebi': None, 'textname': 'N'}, 'T13': {'uniprotID': '', 'textname': 'core histones H2A', 'id': 'corehistonesh2a'}, 'T538': {'chebi': 'CHEBI:28619', 'textname': 'acrylamide'}, 'T18': {'chebi': None, 'textname': 'N'}, 'T19': {'uniprotID': '', 'textname': 'histone', 'id': 'histone'}, 'T539': {'chebi': 'CHEBI:28619', 'textname': 'acrylamide'}, 'T303': {'uniprotID': u'P21675', 'textname': 'TAF1 double bromodomain', 'uniprotprotnames': u'Transcription initiation factor TFIID subunit 1 (EC 2.3.1.48) (EC 2.7.11.1) (Cell cycle gene 1 protein) (TBP-associated factor 250 kDa) (p250) (Transcription initiation factor TFIID 250 kDa subunit) (TAF(II)250) (TAFII-250) (TAFII250)', 'uniprotEntry': u'TAF1_HUMAN', 'id': u'p21675', 'uniprotgenenames': u'TAF1 BA2R CCG1 CCGS TAF2A'}, 'T302': {'uniprotID': '', 'textname': 'acetylated NCPs', 'id': 'acetylatedncp'}, 'T301': {'chebi': 'CHEBI:33704', 'textname': 'amino acid'}, 'T300': {'uniprotID': '', 'textname': 'histone molecules', 'id': 'histonemolecule'}, 'T307': {'uniprotID': '', 'textname': 'TAF1-BrD', 'id': 'taf1brd'}, 'T306': {'uniprotID': '', 'textname': 'TAF1', 'id': 'taf1'}, 'T305': {'uniprotID': '', 'textname': 'Kac', 'id': 'kac'}, 'T304': {'uniprotID': '', 'textname': 'histone', 'id': 'histone'}, 'T309': {'uniprotID': '', 'textname': 'TAF1-BrD', 'id': 'taf1brd'}, 'T308': {'uniprotID': '', 'textname': 'NCP TAF1', 'id': 'ncptaf1'}, 'T442': {'uniprotID': u'P62805', 'textname': 'histone H4 proteins', 'uniprotprotnames': u'Histone H4', 'uniprotEntry': u'H4_HUMAN', 'id': u'p62805', 'uniprotgenenames': u'HIST1H4A H4/A H4FA; HIST1H4B H4/I H4FI; HIST1H4C H4/G H4FG; HIST1H4D H4/B H4FB; HIST1H4E H4/J H4FJ; HIST1H4F H4/C H4FC; HIST1H4H H4/H H4FH; HIST1H4I H4/M H4FM; HIST1H4J H4/E H4FE; HIST1H4K H4/D H4FD; HIST1H4L H4/K H4FK; HIST2H4A H4/N H4F2 H4FN HIST2H4; HIST2H4B H4/O H4FO; HIST4H4'}, 'T534': {'chebi': 'CHEBI:26710', 'textname': 'NaCl'}, 'T443': {'chebi': None, 'textname': 'acetyllysine'}, 'T535': {'chebi': 'CHEBI:64755', 'textname': 'EDTA'}, 'T444': {'uniprotID': u'P62805', 'textname': 'histone H4 proteins', 'uniprotprotnames': u'Histone H4', 'uniprotEntry': u'H4_HUMAN', 'id': u'p62805', 'uniprotgenenames': u'HIST1H4A H4/A H4FA; HIST1H4B H4/I H4FI; HIST1H4C H4/G H4FG; HIST1H4D H4/B H4FB; HIST1H4E H4/J H4FJ; HIST1H4F H4/C H4FC; HIST1H4H H4/H H4FH; HIST1H4I H4/M H4FM; HIST1H4J H4/E H4FE; HIST1H4K H4/D H4FD; HIST1H4L H4/K H4FK; HIST2H4A H4/N H4F2 H4FN HIST2H4; HIST2H4B H4/O H4FO; HIST4H4'}, 'T288': {'uniprotID': '', 'textname': 'H2A', 'id': 'h2a'}, 'T441': {'chebi': None, 'textname': 'acetyllysine'}, 'T282': {'uniprotID': '', 'textname': 'H2B', 'id': 'h2b'}, 'T283': {'uniprotID': '', 'textname': 'molecules C and D', 'id': 'moleculescandd'}, 'T280': {'chebi': None, 'textname': 'phosphate'}, 'T128': {'uniprotID': '', 'textname': 'NCPs', 'id': 'ncp'}, 'T286': {'uniprotID': '', 'textname': 'H4 tail of molecule F', 'id': 'h4tailofmoleculef'}, 'T287': {'uniprotID': '', 'textname': 'molecules G and H', 'id': 'moleculesgandh'}, 'T284': {'uniprotID': '', 'textname': 'H2A', 'id': 'h2a'}, 'T285': {'uniprotID': '', 'textname': 'H2B', 'id': 'h2b'}, 'T123': {'uniprotID': '', 'textname': 'H4 proteins', 'id': 'h4'}, 'T122': {'uniprotID': u'P62805', 'textname': 'histone H4', 'uniprotprotnames': u'Histone H4', 'uniprotEntry': u'H4_HUMAN', 'id': u'p62805', 'uniprotgenenames': u'HIST1H4A H4/A H4FA; HIST1H4B H4/I H4FI; HIST1H4C H4/G H4FG; HIST1H4D H4/B H4FB; HIST1H4E H4/J H4FJ; HIST1H4F H4/C H4FC; HIST1H4H H4/H H4FH; HIST1H4I H4/M H4FM; HIST1H4J H4/E H4FE; HIST1H4K H4/D H4FD; HIST1H4L H4/K H4FK; HIST2H4A H4/N H4F2 H4FN HIST2H4; HIST2H4B H4/O H4FO; HIST4H4'}, 'T121': {'chebi': None, 'textname': 'acetyllysine'}, 'T120': {'uniprotID': '', 'textname': 'S1c and 1d', 'id': 's1cand1d'}, 'T127': {'chebi': 'CHEBI:42820', 'textname': 'guanidine'}, 'T126': {'uniprotID': '', 'textname': 'H4 proteins', 'id': 'h4'}, 'T125': {'uniprotID': '', 'textname': 'S1b', 'id': 's1b'}, 'T124': {'uniprotID': u'Q9HD26', 'textname': 'Fig', 'uniprotprotnames': u'Golgi-associated PDZ and coiled-coil motif-containing protein (CFTR-associated ligand) (Fused in glioblastoma) (PDZ protein interacting specifically with TC10) (PIST)', 'uniprotEntry': u'GOPC_HUMAN', 'id': u'q9hd26', 'uniprotgenenames': u'GOPC CAL FIG'}, 'T532': {'uniprotID': '', 'textname': 'TAF1-BrD', 'id': 'taf1brd'}, 'T531': {'uniprotID': u'P62805', 'textname': 'histone H4 proteins', 'uniprotprotnames': u'Histone H4', 'uniprotEntry': u'H4_HUMAN', 'id': u'p62805', 'uniprotgenenames': u'HIST1H4A H4/A H4FA; HIST1H4B H4/I H4FI; HIST1H4C H4/G H4FG; HIST1H4D H4/B H4FB; HIST1H4E H4/J H4FJ; HIST1H4F H4/C H4FC; HIST1H4H H4/H H4FH; HIST1H4I H4/M H4FM; HIST1H4J H4/E H4FE; HIST1H4K H4/D H4FD; HIST1H4L H4/K H4FK; HIST2H4A H4/N H4F2 H4FN HIST2H4; HIST2H4B H4/O H4FO; HIST4H4'}, 'T533': {'chebi': None, 'textname': 'Tris-HCl'}, 'T549': {'uniprotID': '', 'textname': 'NCPs', 'id': 'ncp'}, 'T548': {'chebi': None, 'textname': 'Mg++'}, 'T239': {'chebi': None, 'textname': 'acetyllysine'}, 'T238': {'chebi': None, 'textname': 'N'}, 'T237': {'chebi': None, 'textname': 'F'}, 'T236': {'uniprotID': '', 'textname': 'K12ac', 'id': 'k12ac'}, 'T235': {'uniprotID': u'Q9HD26', 'textname': 'Fig', 'uniprotprotnames': u'Golgi-associated PDZ and coiled-coil motif-containing protein (CFTR-associated ligand) (Fused in glioblastoma) (PDZ protein interacting specifically with TC10) (PIST)', 'uniprotEntry': u'GOPC_HUMAN', 'id': u'q9hd26', 'uniprotgenenames': u'GOPC CAL FIG'}, 'T234': {'chebi': None, 'textname': 'N'}, 'T233': {'chebi': None, 'textname': 'N'}, 'T232': {'uniprotID': '', 'textname': 'SHL -2', 'id': 'shl2'}, 'T231': {'uniprotID': '', 'textname': 'molecule-D H2B', 'id': 'moleculedh2b'}, 'T230': {'uniprotID': '', 'textname': 'H4 molecules', 'id': 'h4molecule'}, 'T570': {'chebi': 'CHEBI:25555', 'textname': 'N'}, 'T571': {'uniprotID': '', 'textname': 'NCPs', 'id': 'ncp'}, 'T435': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T434': {'chebi': None, 'textname': 'acetyllysine'}, 'T437': {'chebi': None, 'textname': 'NCP'}, 'T436': {'uniprotID': u'P62805', 'textname': 'histone H4', 'uniprotprotnames': u'Histone H4', 'uniprotEntry': u'H4_HUMAN', 'id': u'p62805', 'uniprotgenenames': u'HIST1H4A H4/A H4FA; HIST1H4B H4/I H4FI; HIST1H4C H4/G H4FG; HIST1H4D H4/B H4FB; HIST1H4E H4/J H4FJ; HIST1H4F H4/C H4FC; HIST1H4H H4/H H4FH; HIST1H4I H4/M H4FM; HIST1H4J H4/E H4FE; HIST1H4K H4/D H4FD; HIST1H4L H4/K H4FK; HIST2H4A H4/N H4F2 H4FN HIST2H4; HIST2H4B H4/O H4FO; HIST4H4'}, 'T431': {'uniprotID': u'Q9HD26', 'textname': 'Fig', 'uniprotprotnames': u'Golgi-associated PDZ and coiled-coil motif-containing protein (CFTR-associated ligand) (Fused in glioblastoma) (PDZ protein interacting specifically with TC10) (PIST)', 'uniprotEntry': u'GOPC_HUMAN', 'id': u'q9hd26', 'uniprotgenenames': u'GOPC CAL FIG'}, 'T430': {'uniprotID': '', 'textname': 'TAF1-BrD', 'id': 'taf1brd'}, 'T433': {'uniprotID': '', 'textname': 'NCPs', 'id': 'ncp'}, 'T432': {'uniprotID': '', 'textname': 'S1b', 'id': 's1b'}, 'T439': {'uniprotID': '', 'textname': 'NCPs', 'id': 'ncp'}, 'T438': {'uniprotID': '', 'textname': 'histone DNA', 'id': 'histonedna'}, 'T578': {'uniprotID': u'P21675', 'textname': 'TAF1 double bromodomain', 'uniprotprotnames': u'Transcription initiation factor TFIID subunit 1 (EC 2.3.1.48) (EC 2.7.11.1) (Cell cycle gene 1 protein) (TBP-associated factor 250 kDa) (p250) (Transcription initiation factor TFIID 250 kDa subunit) (TAF(II)250) (TAFII-250) (TAFII250)', 'uniprotEntry': u'TAF1_HUMAN', 'id': u'p21675', 'uniprotgenenames': u'TAF1 BA2R CCG1 CCGS TAF2A'}, 'T446': {'chebi': None, 'textname': 'Kac'}, 'T579': {'uniprotID': '', 'textname': 'TAF1-BrD', 'id': 'taf1brd'}, 'T546': {'chebi': None, 'textname': 'acetyllysine'}, 'T447': {'uniprotID': u'P62805', 'textname': 'human histone H4 cDNA', 'uniprotprotnames': u'Histone H4', 'uniprotEntry': u'H4_HUMAN', 'id': u'p62805', 'uniprotgenenames': u'HIST1H4A H4/A H4FA; HIST1H4B H4/I H4FI; HIST1H4C H4/G H4FG; HIST1H4D H4/B H4FB; HIST1H4E H4/J H4FJ; HIST1H4F H4/C H4FC; HIST1H4H H4/H H4FH; HIST1H4I H4/M H4FM; HIST1H4J H4/E H4FE; HIST1H4K H4/D H4FD; HIST1H4L H4/K H4FK; HIST2H4A H4/N H4F2 H4FN HIST2H4; HIST2H4B H4/O H4FO; HIST4H4'}, 'T420': {'uniprotID': u'P08779', 'textname': 'K16', 'uniprotprotnames': u'Keratin, type I cytoskeletal 16 (Cytokeratin-16) (CK-16) (Keratin-16) (K16)', 'uniprotEntry': u'K1C16_HUMAN', 'id': u'p08779', 'uniprotgenenames': u'KRT16 KRT16A'}, 'T551': {'uniprotID': '', 'textname': 'H4-K5', 'id': 'h4k5'}, 'T445': {'chebi': None, 'textname': 'acetyllysine'}, 'T58': {'chebi': 'CHEBI:25094', 'textname': 'lysine'}, 'T59': {'chebi': None, 'textname': 'N'}, 'T541': {'chebi': 'CHEBI:4883', 'textname': 'ethidium bromide'}, 'T540': {'chebi': None, 'textname': 'bisacrylamide'}, 'T543': {'uniprotID': '', 'textname': 'histone proteins', 'id': 'histone'}, 'T273': {'chebi': None, 'textname': 'N'}, 'T272': {'uniprotID': '', 'textname': 'S3g', 'id': 's3g'}, 'T271': {'uniprotID': '', 'textname': 'molecule-D H2B', 'id': 'moleculedh2b'}, 'T270': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T277': {'chebi': 'CHEBI:49637', 'textname': 'hydrogen'}, 'T276': {'uniprotID': u'Q9NQF3', 'textname': 'SHL +3', 'uniprotprotnames': u'Serine hydrolase-like protein (SHL) (EC 3.1.-.-)', 'uniprotEntry': u'SERHL_HUMAN', 'id': u'q9nqf3', 'uniprotgenenames': u'SERHL SERHL2'}, 'T275': {'uniprotID': '', 'textname': 'SHL -5', 'id': 'shl5'}, 'T274': {'uniprotID': '', 'textname': 'H2B', 'id': 'h2b'}, 'T65': {'chebi': None, 'textname': 'N-dimethyl-aminoethylcysteine'}, 'T64': {'chebi': 'CHEBI:497734', 'textname': 'aminoethylcysteine'}, 'T279': {'chebi': 'CHEBI:49637', 'textname': 'hydrogen'}, 'T278': {'chebi': None, 'textname': 'phosphate'}, 'T61': {'uniprotID': '', 'textname': 'H4-K20', 'id': 'h4k20'}, 'T60': {'uniprotID': '', 'textname': 'H3-K79', 'id': 'h3k79'}, 'T63': {'uniprotID': '', 'textname': 'H4-K20me3', 'id': 'h4k20me3'}, 'T62': {'uniprotID': '', 'textname': 'H3-K79me2', 'id': 'h3k79me2'}, 'T50': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T378': {'uniprotID': u'Q8WVN6', 'textname': 'K12', 'uniprotprotnames': u'Secreted and transmembrane protein 1 (Protein K-12)', 'uniprotEntry': u'SCTM1_HUMAN', 'id': u'q8wvn6', 'uniprotgenenames': u'SECTM1 K12'}, 'T379': {'uniprotID': u'P08779', 'textname': 'K16', 'uniprotprotnames': u'Keratin, type I cytoskeletal 16 (Cytokeratin-16) (CK-16) (Keratin-16) (K16)', 'uniprotEntry': u'K1C16_HUMAN', 'id': u'p08779', 'uniprotgenenames': u'KRT16 KRT16A'}, 'T547': {'uniprotID': '', 'textname': 'NCPs', 'id': 'ncp'}, 'T269': {'uniprotID': '', 'textname': 'S3d', 'id': 's3d'}, 'T372': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T373': {'uniprotID': u'Q9NQF3', 'textname': 'SHL', 'uniprotprotnames': u'Serine hydrolase-like protein (SHL) (EC 3.1.-.-)', 'uniprotEntry': u'SERHL_HUMAN', 'id': u'q9nqf3', 'uniprotgenenames': u'SERHL SERHL2'}, 'T370': {'uniprotID': '', 'textname': 'SHL -2', 'id': 'shl2'}, 'T371': {'chebi': 'CHEBI:36976', 'textname': 'nucleotides'}, 'T376': {'uniprotID': '', 'textname': 'H4-K5', 'id': 'h4k5'}, 'T183': {'uniprotID': '', 'textname': 'S2b', 'id': 's2b'}, 'T374': {'uniprotID': '', 'textname': 'H2A1', 'id': 'h2a1'}, 'T375': {'chebi': None, 'textname': 'N'}, 'T152': {'uniprotID': '', 'textname': 'H4 proteins', 'id': 'h4'}, 'T153': {'uniprotID': '', 'textname': 'Fig. 1b', 'id': 'fig.1b'}, 'T150': {'uniprotID': '', 'textname': 'NCPs454647', 'id': 'ncps454647'}, 'T151': {'chebi': None, 'textname': 'Mg++'}, 'T156': {'chebi': None, 'textname': 'Mg++'}, 'T157': {'uniprotID': '', 'textname': 'K5', 'id': 'k5'}, 'T154': {'chebi': None, 'textname': 'Mg++'}, 'T155': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T404': {'chebi': None, 'textname': 'N'}, 'T54': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T158': {'uniprotID': u'P05787', 'textname': 'K8', 'uniprotprotnames': u'Keratin, type II cytoskeletal 8 (Cytokeratin-8) (CK-8) (Keratin-8) (K8) (Type-II keratin Kb8)', 'uniprotEntry': u'K2C8_HUMAN', 'id': u'p05787', 'uniprotgenenames': u'KRT8 CYK8'}, 'T159': {'uniprotID': u'Q8WVN6', 'textname': 'K12', 'uniprotprotnames': u'Secreted and transmembrane protein 1 (Protein K-12)', 'uniprotEntry': u'SCTM1_HUMAN', 'id': u'q8wvn6', 'uniprotgenenames': u'SECTM1 K12'}, 'T400': {'chebi': None, 'textname': 'N'}, 'T401': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T402': {'uniprotID': '', 'textname': 'H2A', 'id': 'h2a'}, 'T55': {'uniprotID': '', 'textname': 'core histones', 'id': 'corehistone'}, 'T56': {'uniprotID': '', 'textname': 'histone histone', 'id': 'histonehistone'}, 'T57': {'uniprotID': '', 'textname': 'histone DNA', 'id': 'histonedna'}, 'T479': {'uniprotID': '', 'textname': 'H4 protein', 'id': 'h4'}, 'T478': {'uniprotID': u'P08779', 'textname': 'K16', 'uniprotprotnames': u'Keratin, type I cytoskeletal 16 (Cytokeratin-16) (CK-16) (Keratin-16) (K16)', 'uniprotEntry': u'K1C16_HUMAN', 'id': u'p08779', 'uniprotgenenames': u'KRT16 KRT16A'}, 'T572': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T573': {'uniprotID': '', 'textname': 'Nucleosomal histone', 'id': 'nucleosomalhistone'}, 'T574': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T575': {'uniprotID': '', 'textname': 'K16ac-D24', 'id': 'k16acd24'}, 'T576': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T577': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T471': {'chebi': 'CHEBI:14434', 'textname': 'imidazole'}, 'T470': {'chebi': 'CHEBI:14434', 'textname': 'imidazole'}, 'T473': {'chebi': 'CHEBI:41218', 'textname': '2-mercaptoethanol'}, 'T472': {'uniprotID': '', 'textname': 'H4 proteins', 'id': 'h4'}, 'T475': {'uniprotID': '', 'textname': 'K5', 'id': 'k5'}, 'T474': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T477': {'uniprotID': u'Q8WVN6', 'textname': 'K12', 'uniprotprotnames': u'Secreted and transmembrane protein 1 (Protein K-12)', 'uniprotEntry': u'SCTM1_HUMAN', 'id': u'q8wvn6', 'uniprotgenenames': u'SECTM1 K12'}, 'T476': {'uniprotID': u'P05787', 'textname': 'K8', 'uniprotprotnames': u'Keratin, type II cytoskeletal 8 (Cytokeratin-8) (CK-8) (Keratin-8) (K8) (Type-II keratin Kb8)', 'uniprotEntry': u'K2C8_HUMAN', 'id': u'p05787', 'uniprotgenenames': u'KRT8 CYK8'}, 'T259': {'chebi': None, 'textname': 'N'}, 'T258': {'uniprotID': '', 'textname': 'S3g', 'id': 's3g'}, 'T497': {'uniprotID': '', 'textname': 'NCPs', 'id': 'ncp'}, 'T496': {'uniprotID': '', 'textname': 'NCPs', 'id': 'ncp'}, 'T495': {'chebi': 'CHEBI:6636', 'textname': 'MgCl2'}, 'T494': {'chebi': None, 'textname': 'acetyllysine'}, 'T493': {'uniprotID': '', 'textname': 'human alpha-satellite DNA region', 'id': 'humanalphasatellitednaregion'}, 'T492': {'uniprotID': '', 'textname': 'H2B', 'id': 'h2b'}, 'T491': {'uniprotID': '', 'textname': 'histones H2A', 'id': 'histonesh2a'}, 'T490': {'uniprotID': '', 'textname': 'H4 proteins', 'id': 'h4'}, 'T499': {'chebi': 'CHEBI:64755', 'textname': 'EDTA'}, 'T498': {'chebi': None, 'textname': 'Tris-HCl'}, 'T342': {'uniprotID': u'P08779', 'textname': 'K16', 'uniprotprotnames': u'Keratin, type I cytoskeletal 16 (Cytokeratin-16) (CK-16) (Keratin-16) (K16)', 'uniprotEntry': u'K1C16_HUMAN', 'id': u'p08779', 'uniprotgenenames': u'KRT16 KRT16A'}, 'T318': {'uniprotID': '', 'textname': 'TAF1-BrD', 'id': 'taf1brd'}, 'T319': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T223': {'uniprotID': '', 'textname': 'NCPs', 'id': 'ncp'}, 'T281': {'uniprotID': '', 'textname': 'histones H2A', 'id': 'histonesh2a'}, 'T29': {'uniprotID': '', 'textname': 'H2B', 'id': 'h2b'}, 'T28': {'chebi': None, 'textname': 'N'}, 'T21': {'chebi': None, 'textname': 'N'}, 'T20': {'uniprotID': '', 'textname': 'histone octamer structure', 'id': 'histoneoctamerstructure'}, 'T23': {'chebi': 'CHEBI:29016', 'textname': 'arginine'}, 'T22': {'chebi': 'CHEBI:25094', 'textname': 'lysine'}, 'T25': {'uniprotID': '', 'textname': 'core histones', 'id': 'corehistone'}, 'T24': {'uniprotID': '', 'textname': 'human4, NCPs containing histone variants567', 'id': 'human4,ncpscontaininghistonevariants567'}, 'T27': {'uniprotID': '', 'textname': 'NCP19', 'id': 'ncp19'}, 'T26': {'chebi': None, 'textname': 'N'}, 'T338': {'uniprotID': u'P05787', 'textname': 'K8', 'uniprotprotnames': u'Keratin, type II cytoskeletal 8 (Cytokeratin-8) (CK-8) (Keratin-8) (K8) (Type-II keratin Kb8)', 'uniprotEntry': u'K2C8_HUMAN', 'id': u'p05787', 'uniprotgenenames': u'KRT8 CYK8'}, 'T339': {'uniprotID': u'Q8WVN6', 'textname': 'K12', 'uniprotprotnames': u'Secreted and transmembrane protein 1 (Protein K-12)', 'uniprotEntry': u'SCTM1_HUMAN', 'id': u'q8wvn6', 'uniprotgenenames': u'SECTM1 K12'}, 'T336': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T337': {'uniprotID': '', 'textname': 'K5', 'id': 'k5'}, 'T334': {'uniprotID': '', 'textname': 'Fig. 1a', 'id': 'fig.1a'}, 'T335': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T332': {'uniprotID': '', 'textname': 'H4 tetra-acetylated and unmodified NCPs', 'id': 'h4tetraacetylatedandunmodifiedncp'}, 'T333': {'uniprotID': '', 'textname': 'H4 tetra-acetylation', 'id': 'h4tetraacetylation'}, 'T330': {'uniprotID': '', 'textname': 'histone', 'id': 'histone'}, 'T331': {'chebi': None, 'textname': 'N'}, 'T53': {'uniprotID': '', 'textname': 'K16ac', 'id': 'k16ac'}, 'T299': {'uniprotID': '', 'textname': 'histone histone', 'id': 'histonehistone'}, 'T298': {'uniprotID': '', 'textname': 'S3e', 'id': 's3e'}, 'T529': {'uniprotID': '', 'textname': 'TAF1-BrD', 'id': 'taf1brd'}, 'T528': {'uniprotID': u'P21675', 'textname': 'TAF1 double bromodomain', 'uniprotprotnames': u'Transcription initiation factor TFIID subunit 1 (EC 2.3.1.48) (EC 2.7.11.1) (Cell cycle gene 1 protein) (TBP-associated factor 250 kDa) (p250) (Transcription initiation factor TFIID 250 kDa subunit) (TAF(II)250) (TAFII-250) (TAFII250)', 'uniprotEntry': u'TAF1_HUMAN', 'id': u'p21675', 'uniprotgenenames': u'TAF1 BA2R CCG1 CCGS TAF2A'}, 'T295': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T294': {'uniprotID': '', 'textname': 'S3e', 'id': 's3e'}, 'T297': {'uniprotID': u'Q9HD26', 'textname': 'Fig', 'uniprotprotnames': u'Golgi-associated PDZ and coiled-coil motif-containing protein (CFTR-associated ligand) (Fused in glioblastoma) (PDZ protein interacting specifically with TC10) (PIST)', 'uniprotEntry': u'GOPC_HUMAN', 'id': u'q9hd26', 'uniprotgenenames': u'GOPC CAL FIG'}, 'T296': {'uniprotID': '', 'textname': 'histone core region', 'id': 'histonecoreregion'}, 'T291': {'uniprotID': '', 'textname': 'molecule-D H2B', 'id': 'moleculedh2b'}, 'T290': {'uniprotID': u'A5YKK6', 'textname': 'not1', 'uniprotprotnames': u'CCR4-NOT transcription complex subunit 1 (CCR4-associated factor 1) (Negative regulator of transcription subunit 1 homolog) (NOT1H) (hNOT1)', 'uniprotEntry': u'CNOT1_HUMAN', 'id': u'a5ykk6', 'uniprotgenenames': u'CNOT1 CDC39 KIAA1007 NOT1 AD-005'}, 'T293': {'uniprotID': u'Q9HD26', 'textname': 'Fig', 'uniprotprotnames': u'Golgi-associated PDZ and coiled-coil motif-containing protein (CFTR-associated ligand) (Fused in glioblastoma) (PDZ protein interacting specifically with TC10) (PIST)', 'uniprotEntry': u'GOPC_HUMAN', 'id': u'q9hd26', 'uniprotgenenames': u'GOPC CAL FIG'}, 'T292': {'uniprotID': '', 'textname': 'V21', 'id': 'v21'}, 'T422': {'chebi': None, 'textname': 'N'}, 'T208': {'uniprotID': '', 'textname': 'Fig. 2e', 'id': 'fig.2e'}, 'T209': {'uniprotID': '', 'textname': 'E56', 'id': 'e56'}, 'T423': {'uniprotID': '', 'textname': 'NCPs', 'id': 'ncp'}, 'T202': {'uniprotID': '', 'textname': 'Fig. 2g', 'id': 'fig.2g'}, 'T203': {'uniprotID': '', 'textname': 'H4 molecule', 'id': 'h4molecule'}, 'T200': {'uniprotID': '', 'textname': 'H4 molecules', 'id': 'h4molecule'}, 'T201': {'uniprotID': '', 'textname': 'molecule-B H4', 'id': 'moleculebh4'}, 'T206': {'chebi': 'CHEBI:40574', 'textname': 'acetyl'}, 'T207': {'uniprotID': '', 'textname': 'NCPs', 'id': 'ncp'}, 'T204': {'uniprotID': '', 'textname': 'NCPs', 'id': 'ncp'}, 'T205': {'uniprotID': u'P08779', 'textname': 'K16', 'uniprotprotnames': u'Keratin, type I cytoskeletal 16 (Cytokeratin-16) (CK-16) (Keratin-16) (K16)', 'uniprotEntry': u'K1C16_HUMAN', 'id': u'p08779', 'uniprotgenenames': u'KRT16 KRT16A'}, 'T421': {'uniprotID': '', 'textname': 'TAF1-BrD', 'id': 'taf1brd'}, 'T426': {'uniprotID': '', 'textname': 'H4 N-terminal tail', 'id': 'h4nterminaltail'}, 'T427': {'uniprotID': '', 'textname': 'TAF1-BrD', 'id': 'taf1brd'}, 'T590': {'chebi': None, 'textname': 'Mn++'}, 'T424': {'uniprotID': '', 'textname': 'TAF1-BrD', 'id': 'taf1brd'}, 'T425': {'uniprotID': '', 'textname': 'TAF1-BrD', 'id': 'taf1brd'}, 'T116': {'uniprotID': u'P62805', 'textname': 'histone H4', 'uniprotprotnames': u'Histone H4', 'uniprotEntry': u'H4_HUMAN', 'id': u'p62805', 'uniprotgenenames': u'HIST1H4A H4/A H4FA; HIST1H4B H4/I H4FI; HIST1H4C H4/G H4FG; HIST1H4D H4/B H4FB; HIST1H4E H4/J H4FJ; HIST1H4F H4/C H4FC; HIST1H4H H4/H H4FH; HIST1H4I H4/M H4FM; HIST1H4J H4/E H4FE; HIST1H4K H4/D H4FD; HIST1H4L H4/K H4FK; HIST2H4A H4/N H4F2 H4FN HIST2H4; HIST2H4B H4/O H4FO; HIST4H4'}, 'T117': {'uniprotID': u'Q9HD26', 'textname': 'Fig', 'uniprotprotnames': u'Golgi-associated PDZ and coiled-coil motif-containing protein (CFTR-associated ligand) (Fused in glioblastoma) (PDZ protein interacting specifically with TC10) (PIST)', 'uniprotEntry': u'GOPC_HUMAN', 'id': u'q9hd26', 'uniprotgenenames': u'GOPC CAL FIG'}, 'T114': {'uniprotID': u'Q8WVN6', 'textname': 'K12', 'uniprotprotnames': u'Secreted and transmembrane protein 1 (Protein K-12)', 'uniprotEntry': u'SCTM1_HUMAN', 'id': u'q8wvn6', 'uniprotgenenames': u'SECTM1 K12'}, 'T115': {'uniprotID': u'P08779', 'textname': 'K16', 'uniprotprotnames': u'Keratin, type I cytoskeletal 16 (Cytokeratin-16) (CK-16) (Keratin-16) (K16)', 'uniprotEntry': u'K1C16_HUMAN', 'id': u'p08779', 'uniprotgenenames': u'KRT16 KRT16A'}, 'T112': {'uniprotID': '', 'textname': 'K5', 'id': 'k5'}, 'T113': {'uniprotID': u'P05787', 'textname': 'K8', 'uniprotprotnames': u'Keratin, type II cytoskeletal 8 (Cytokeratin-8) (CK-8) (Keratin-8) (K8) (Type-II keratin Kb8)', 'uniprotEntry': u'K2C8_HUMAN', 'id': u'p05787', 'uniprotgenenames': u'KRT8 CYK8'}, 'T110': {'uniprotID': '', 'textname': 'RF1', 'id': 'rf1'}, 'T111': {'uniprotID': u'P62805', 'textname': 'histone H4 proteins', 'uniprotprotnames': u'Histone H4', 'uniprotEntry': u'H4_HUMAN', 'id': u'p62805', 'uniprotgenenames': u'HIST1H4A H4/A H4FA; HIST1H4B H4/I H4FI; HIST1H4C H4/G H4FG; HIST1H4D H4/B H4FB; HIST1H4E H4/J H4FJ; HIST1H4F H4/C H4FC; HIST1H4H H4/H H4FH; HIST1H4I H4/M H4FM; HIST1H4J H4/E H4FE; HIST1H4K H4/D H4FD; HIST1H4L H4/K H4FK; HIST2H4A H4/N H4F2 H4FN HIST2H4; HIST2H4B H4/O H4FO; HIST4H4'}, 'T448': {'chebi': 'CHEBI:16443', 'textname': 'TAG'}, 'T449': {'chebi': None, 'textname': 'Invitrogen'}, 'T118': {'uniprotID': '', 'textname': 'S1b', 'id': 's1b'}, 'T119': {'uniprotID': u'Q9HD26', 'textname': 'Fig', 'uniprotprotnames': u'Golgi-associated PDZ and coiled-coil motif-containing protein (CFTR-associated ligand) (Fused in glioblastoma) (PDZ protein interacting specifically with TC10) (PIST)', 'uniprotEntry': u'GOPC_HUMAN', 'id': u'q9hd26', 'uniprotgenenames': u'GOPC CAL FIG'}, 'T289': {'uniprotID': '', 'textname': 'H2B', 'id': 'h2b'}, 'T417': {'uniprotID': '', 'textname': 'K5', 'id': 'k5'}, 'T503': {'uniprotID': '', 'textname': 'NCPs', 'id': 'ncp'}, 'T8': {'uniprotID': u'P05787', 'textname': 'K8', 'uniprotprotnames': u'Keratin, type II cytoskeletal 8 (Cytokeratin-8) (CK-8) (Keratin-8) (K8) (Type-II keratin Kb8)', 'uniprotEntry': u'K2C8_HUMAN', 'id': u'p05787', 'uniprotgenenames': u'KRT8 CYK8'}, 'T9': {'uniprotID': u'Q8WVN6', 'textname': 'K12', 'uniprotprotnames': u'Secreted and transmembrane protein 1 (Protein K-12)', 'uniprotEntry': u'SCTM1_HUMAN', 'id': u'q8wvn6', 'uniprotgenenames': u'SECTM1 K12'}, 'T6': {'uniprotID': '', 'textname': 'NCPs', 'id': 'ncp'}, 'T7': {'uniprotID': '', 'textname': 'H4-K5', 'id': 'h4k5'}, 'T4': {'uniprotID': '', 'textname': 'histone', 'id': 'histone'}, 'T5': {'uniprotID': '', 'textname': 'tetra-acetylated nucleosome core particles', 'id': 'tetraacetylatednucleosomecoreparticle'}, 'T2': {'chebi': 'CHEBI:25094', 'textname': 'lysine'}, 'T3': {'chebi': None, 'textname': 'N'}, 'T1': {'uniprotID': u'P62805', 'textname': 'histone H4', 'uniprotprotnames': u'Histone H4', 'uniprotEntry': u'H4_HUMAN', 'id': u'p62805', 'uniprotgenenames': u'HIST1H4A H4/A H4FA; HIST1H4B H4/I H4FI; HIST1H4C H4/G H4FG; HIST1H4D H4/B H4FB; HIST1H4E H4/J H4FJ; HIST1H4F H4/C H4FC; HIST1H4H H4/H H4FH; HIST1H4I H4/M H4FM; HIST1H4J H4/E H4FE; HIST1H4K H4/D H4FD; HIST1H4L H4/K H4FK; HIST2H4A H4/N H4F2 H4FN HIST2H4; HIST2H4B H4/O H4FO; HIST4H4'}, 'T416': {'uniprotID': '', 'textname': 'H4', 'id': 'h4'}, 'T467': {'chebi': 'CHEBI:42820', 'textname': 'guanidine'}, 'T509': {'chebi': 'CHEBI:6636', 'textname': 'MgCl2'}, 'T508': {'chebi': 'CHEBI:6636', 'textname': 'MgCl2'}, 'T545': {'uniprotID': '', 'textname': 'histone proteins', 'id': 'histone'}}}
	entity_recognition.process_a1_changed(path, pmid, dictmchem, dicoPMID_idDoc, '', compareDB)
	#entity_recognition.process_a2(path, pmid, 'search_Type', '')