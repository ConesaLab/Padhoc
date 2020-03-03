#!/usr/bin/env python

'''
Download the whole neo4j database
and match using the comparison algorithm to obtain the 
percentage of the entity that matches
'''


import time, sys, re
from sys import argv

from chebi_from_string import chebi_from_string
from uniprot_queries import uniprot_queries
from brenda_annotation import brenda_annotation
from modules import calculateSimilarity
from graph_database import graph_database


class compare_string_to_db():
	'''
	'''

	def __init__(self, neo4j):

		self.neo4j = neo4j
		self.proteinsDictio = {}
		self.compoundsDictio = {}


	def obtain_neo4j(self):
		'''
		Obtain neo4j db using a cypher request. Here we
		will only ask for the properties that we are interested,

		'''
		request = 'MATCH (n:Compound) RETURN n.chebiID as id, n.compoundName as name, n.textname as text'
		compoundResult = self.neo4j.session.run(request)
		request = 'MATCH (n:Protein) RETURN n.uniprotID as id, n.synonyms as names, n.uniprotProteinNames as uniNames'
		proteinResult = self.neo4j.session.run(request)
		return compoundResult, proteinResult

	def obtain_queried_neo4j(self, protein_id):
		'''
		'''
		request = 'MATCH (n:Protein) WHERE n.id="%s" RETURN n.uniprotID as id, n.synonyms as names, n.uniprotProteinNames as uniNames'%(protein_id)
		proteinResult = self.neo4j.session.run(request)
		request = 'MATCH (y:Protein)-[r]-(n:Compound) WHERE y.id="%s" RETURN n.chebiID as id, n.compoundName as name, n.textname as text'%(protein_id)
		compoundResult = self.neo4j.session.run(request)
		return compoundResult, proteinResult


	def neo4j_to_dictio(self, protein_id = ""):
		'''
		Request the whole 
		'''
		if protein_id == "":
			compoundResult, proteinResult = self.obtain_neo4j()
		else:
			compoundResult, proteinResult = self.obtain_queried_neo4j(protein_id)
		for r in compoundResult:
			chebiId = r['id']; chebiName = r['name']
			if chebiName != None:
				self.compoundsDictio[chebiName] = chebiId 
		for r in proteinResult:
			uniID = r['id']; synonyms = r['names'];
			uniprotNames = r['uniNames']
			if synonyms != None:
				for syn in synonyms:
					self.proteinsDictio[syn] = uniID
			elif uniprotNames != None:
				for name in uniprotNames:
					self.proteinsDictio[name] = uniID
		return None


	def delete_enzyme_action(self, entity):
		'''
		Delete the action of the enzyme that is being studied,
		this is done to compare the target of the enzyme
		'''
		entity = entity.lower()
		entity = re.sub('phosphatase.*', '', entity)
		entity = re.sub('carbamylase.*', '', entity)
		entity = re.sub('transaminase.*', '', entity)
		entity = re.sub('dehydrogenase.*', '', entity)
		entity = re.sub('synthetase.*', '', entity)
		entity = re.sub('synthase.*', '', entity)
		entity = re.sub('hydrolase.*', '', entity)
		entity = re.sub('lyase.*', '', entity)
		entity = re.sub('kinase.*', '', entity)
		entity = re.sub('deaminase.*', '', entity)
		entity = re.sub('oxidase.*', '', entity)
		entity = re.sub('lipase.*', '', entity)
		entity = re.sub('decarboxylase.*', '', entity)
		entity = re.sub('isomerase.*', '', entity)
		entity = re.sub('racemase.*', '', entity)
		entity = re.sub('oxidoreductase.*', '', entity)
		entity = re.sub('ligase.*', '', entity)
		entity = re.sub('transferase.*', '', entity)
		entity = re.sub('reductase.*', '', entity)
		entity = re.sub('Transcription factor.*', '', entity)
		entity = re.sub('dioxygenase', '', entity)
		entity = re.sub('epimerase.*', '', entity)
		return entity


	def search_in_dictio(self, entity, threshold = 0.8):
		'''
		Iterate over the dictionaries, and search by similarity the
		entity given as input
		'''
		compounds_list = []; proteins_list = []
		for k in self.compoundsDictio.keys():
			similarity = calculateSimilarity(self.clean(entity), self.clean(k))
#				print entity, k
			if similarity != 1.0:
				if (len(self.clean(entity).replace("ate", "")) <= 4) or (len(self.clean(k).replace("ate", "")) <= 4):
					similarity = similarity - 0.3
				if (len(self.clean(entity).replace("ylcoa", "")) <= 4) or (len(self.clean(k).replace("ylcoa", "")) <= 4):
					similarity = similarity - 0.3
			else:
				if (len(self.clean(entity).replace("ate", "")) <= 3) or (len(self.clean(k).replace("ate", "")) <= 3):
					similarity = similarity - 0.3
				if (len(self.clean(entity).replace("ylcoa", "")) <= 3) or (len(self.clean(k).replace("ylcoa", "")) <= 3):
					similarity = similarity - 0.3
			if similarity > threshold:
				compounds_list.append((k, similarity))
		# First delete the property of the enzyme, we want to compare the target
		for k in self.proteinsDictio.keys():
			if k != '':
				entity_mod = self.delete_enzyme_action(entity)
				k_mod = self.delete_enzyme_action(k)
		# Then compare everything, including the action
				if k_mod != '' and entity_mod != '':
					similarity = calculateSimilarity(self.clean(entity_mod), self.clean(k_mod))
					if similarity > threshold:
						similarity = calculateSimilarity(self.clean(entity), self.clean(k))
						if similarity != 1.0:
							if (len(self.clean(entity)) <= 4) or (len(self.clean(k)) <= 4):
								similarity = similarity - 0.3
						else:
							if (len(self.clean(entity)) <= 4) or (len(self.clean(k)) <= 4):
								similarity = similarity - 0.3
	#					print entity, k
						if similarity > threshold:
							proteins_list.append((k, similarity))
		return compounds_list, proteins_list


	def obtain_best_match(self, tuples_list):
		'''
		Order the elements of the list based on similarity,
		obtain the best match
		'''
		sorted_list = sorted(tuples_list,key=lambda x:(-x[1],x[0]))
		best_match = sorted_list[0]
		return best_match


	def compare_lists(self, compounds, proteins):
		'''
		Returns True if the compounds ist has a better best result,
		False if the best match is a protein
		'''
		# Order the lists based on the similarity
		best_comp = self.obtain_best_match(compounds)
		best_prot = self.obtain_best_match(proteins)
		if best_comp[1] > best_prot[1]:
			return True
		else:
			return False


	def obtain_list(self, elements_list, is_compound):
		'''
		Store in a list all the nodes found
		'''
		nodes_list = []
		for elem in elements_list:
			if is_compound == True:
				node_id = self.compoundsDictio[elem[0]]
				if node_id not in nodes_list:
					nodes_list.append(node_id)
			else:
				node_id = self.proteinsDictio[elem[0]]
				if node_id not in nodes_list:
					nodes_list.append(node_id)
		return nodes_list

	def similar_nodes(self, entity):
		'''
		'''
		start = time.time()
		ec_tester = ''.join([i for i in entity.replace('.','') if not i.isdigit()])
		if ec_tester != '':
			compounds_list, protein_list = self.search_in_dictio(entity)
		else:
			compounds_list = []
			protein_list = []
		if len(protein_list) > 0 and len(compounds_list) > 0:
			comparison = self.compare_lists(compounds_list, protein_list)
			if comparison == True:
				protein_list = []
			elif comparison == False:
				compounds_list = []
		if len(protein_list) == 0:
			#Recover the metabolite nodes
			nodes_ids = self.obtain_list(compounds_list, True)
			is_metabolite = True
		if len(compounds_list) == 0:
			#Recover the protein nodes
			nodes_ids = self.obtain_list(protein_list, False)
			is_metabolite = False
		if len(protein_list) == 0 and len(compounds_list) == 0:
			nodes_ids = []; is_metabolite = None
			#print 'No match found with database'
		if len(protein_list) > 20 or len(compounds_list) == 20:		#In case it unespecifically finds too many nodes from the database
			nodes_ids = []; is_metabolite = None
		return nodes_ids, is_metabolite


	def clean(self, word,specie=''):
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
		cleanw=re.sub('&','',cleanw)
		cleanw=re.sub(';','',cleanw)
		cleanw=cleanw.replace(",",'')
		cleanw=re.sub('<[^>]+>','',cleanw)
		cleanw=re.sub('^(r)-','',cleanw)
		cleanw=re.sub('^(s)-','',cleanw)
		cleanw=re.sub('^N-acetyl-','',cleanw)
		cleanw=re.sub('\)n$','',cleanw)
		cleanw=re.sub('\)m$','',cleanw)
		cleanw=re.sub('^\(\+?\-?\d{0,5}\)','',cleanw)
		cleanw=re.sub('\(','',cleanw)
		cleanw=re.sub('\)','',cleanw)
		cleanw=re.sub('\[','',cleanw)
		cleanw=re.sub('\]','',cleanw)
		cleanw=re.sub('\{','',cleanw)
		cleanw=re.sub('\}','',cleanw)
		cleanw=re.sub('^[0-9]+','',cleanw)
		cleanw=re.sub('^l-','',cleanw)
		cleanw=re.sub('^d-','',cleanw)
		cleanw=re.sub('^r-','',cleanw)
		cleanw=re.sub('^s-','',cleanw)
		cleanw=re.sub('^ec:','',cleanw)
		cleanw=re.sub('^ec ','',cleanw)
		cleanw=re.sub('^alpha-','',cleanw)
		cleanw=re.sub('^beta-','',cleanw)
		cleanw=re.sub('^[0-9]+','',cleanw)
		cleanw=re.sub(' genes{0,1}','',cleanw)
		cleanw=re.sub(' proteins{0,1}','',cleanw)
		cleanw=re.sub('beta:','',cleanw)
		cleanw=re.sub('gamma:','',cleanw)
		#cleanw=re.sub('\)','',cleanw)
		#remove words never used in the text
		cleanw=re.sub(' atom','',cleanw)
		cleanw=re.sub('s$','',cleanw) #remove plural
		cleanw=re.sub("ine$",'',cleanw)
		cleanw=re.sub("ase$",'',cleanw)
		cleanw=re.sub("ose$",'',cleanw)
		cleanw=cleanw.replace(' ','') #test 3 july
		cleanw=cleanw.replace('"','')
		cleanw=cleanw.replace('-','')
		cleanw=cleanw.replace('+','')
		cleanw=cleanw.replace("'",'')
		cleanw=cleanw.replace("icacid",'ate')
		cleanw=cleanw.replace("ate",'')
		cleanw=cleanw.replace("complex",'')
		cleanw = re.sub("\d+$", "", cleanw)
		for s in sp:
		    cleanw=re.sub(s,'',cleanw)#remove the specie name from the entity name (ex for Arabidopsis thaliana: remove arabidopsis and then remove thaliana)
		if cleanw=='':
		    cleanw='empty_string'
		return cleanw


if __name__ == '__main__':
	
	start = time.time()
	#file = open(argv[1]).read()

	chebi = chebi_from_string()
	chebi.chebi_connect()
	uniprot = uniprot_queries('Homo sapiens', '9606')

	brenda = brenda_annotation('salcagal@alumni.uv.es', 'salvacasani91') # create brenda_annotation object
	brenda.access_protocol() # access to brenda
	neo4j = graph_database(chebi, uniprot, brenda,'Homo sapiens','neo4j','salva') #
	neo4j.connect()

	compareDB = compare_string_to_db(neo4j)
	compareDB.neo4j_to_dictio()

	print compareDB.clean("PanD complex")
	print compareDB.similar_nodes('rpoB')
	sys.exit()
	print compareDB.clean("pantetheine 4'-phosphate")
	print compareDB.search_in_dictio("pantetheine 4'-phosphate")
	sys.exit()
	print compareDB.clean('ATP')
	print compareDB.search_in_dictio('ATP')
	sys.exit()
	print compareDB.search_in_dictio('D-alanyl-D-alanine')
	sys.exit()

	a = []
	n = 0
	for line in file.split('\n')[1:]:
		print line.split('\t')[0]
		ids = compareDB.similar_nodes(line.split('\t')[0])[0]
		is_metabolite = compareDB.similar_nodes(line.split('\t')[0])[1]
		if is_metabolite == True:
			for elem in ids:
				if elem not in a:
					a.append(elem)
		else:
			print ids, is_metabolite
			n += 1

	print a
	print n
	sys.exit()

	print compareDB.clean('WT protein')
	print compareDB.search_in_dictio('WT protein')
	sys.exit()
#	print compareDB.clean('(6S)-5,6,7,8-tetrahydrofolic acid')
#	print compareDB.clean('tetrahydrofolate')
#	print calculateSimilarity(compareDB.clean('(6S)-5,6,7,8-tetrahydrofolic acid'), compareDB.clean('tetrahydrofolate'))
#	print calculateSimilarity(compareDB.clean('phnP'), compareDB.clean('PhP'))
	print calculateSimilarity(compareDB.clean("alcohol (NADP+)"), compareDB.clean("alcohol"))
	sys.exit()
	nodes_list = compareDB.similar_nodes("tetrahydrofolate")
	#nodes_list = compareDB.similar_nodes('N-succinyl-L,L-diaminopimelic acid desuccinylase')

	print nodes_list
	sys.exit()


	end = time.time()
	print end - start