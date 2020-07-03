#!/usr/bin/env python

from __future__ import division
import sys
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn import metrics

from ..graph_database import graph_database
from ..chebi_from_string import chebi_from_string
from ..compare_string_to_db import compare_string_to_db
from ..modules import calculateSimilarity
from ..uniprot_queries import uniprot_queries


class compress_graph():
	'''
	'''

	def __init__(self, neo4j):
		
		self.neo4j = neo4j
		self.rels_dictio = {}
		self.compounds_dictio = {}
		self.proteins_dictio = {}

		self.element_translator = {}


	def neo4j_to_dictionary(self, pattern):
		'''
		Create the dictionaries to extract the neo4j database
		into python object
		'''
		self.neo4j.session.run('MATCH (n) WHERE "uniprotID" IN keys(n) SET n.database = TRUE')
		self.neo4j.session.run('MATCH (n) WHERE "chebiID" IN keys(n) SET n.database = TRUE')
		result = self.neo4j.extract_pattern(pattern)
		for r in result:
			nID = r['nID']; nSent = r['nSentences']
			yID = r['yID']; ySent = r['ySentences']
			nLabel = r['nlabel']; yLabel = r['ylabel']
			rSent = r['rSentences']
			yProtNames = r['yuniprotProteinNames']; ySyn = r['ysyn']
			nProtNames = r['nuniprotProteinNames']; nSyn = r['nsyn']
			ySyn = self.create_syns_list(ySyn, yProtNames)
			nSyn = self.create_syns_list(nSyn, nProtNames)
			nCName = [r['nCompound']]; yCName = [r['yCompound']]

			if r['nECs'] == None:
				nECs = []
			else:
				nECs = r['nECs']
			if r['yECs'] == None:
				yECs = []
			else:
				yECs = r['yECs']
			if nID not in self.rels_dictio.keys():
				if rSent == None:
					rSents = set([])
				else:
					rSents = set([])
					for sentence in rSent:
						rSents.add(sentence.encode('utf8'))
				self.rels_dictio[nID] = {yID:rSents}
			else:
				if rSent != None:
					if yID in self.rels_dictio[nID].keys():
						for sentence in rSent:
							self.rels_dictio[nID][yID].add(sentence.encode('utf8'))
				else:
					if rSent == None:
						rSents = set([])
					else:
						rSents = set([])
						for sentence in rSent:
							rSents.add(sentence.encode('utf8'))
					self.rels_dictio[nID][yID] = set(rSents)

			# Check the type of the node
			self.determine_node(nID, nLabel, nSent, nSyn, nCName, nECs)
			self.determine_node(yID, yLabel, ySent, ySyn, yCName, yECs)

			self.fill_translator(nID, nLabel, nSyn, nCName)
			self.fill_translator(yID, yLabel, ySyn, yCName)
		return None


	def create_syns_list(self, synonyms, protein_names):
		'''
		Create a common list between synonyms and protein names
		'''
		if synonyms == None:
			synonyms = []
		if protein_names == None:
			protein_names = []
		total_synonyms = synonyms + protein_names
		if total_synonyms == []:
			total_synonyms = None
		return total_synonyms


	def fill_translator(self, node_id, label, synonyms, compoundName):
		'''
		Create a "translator" in which the element id is
		inputed, and outputs the list of names that it has
		'''
		if 'Compound' in label:
			if node_id not in self.element_translator.keys():
				self.element_translator[node_id] = compoundName
		if 'Protein' in label:
			if node_id not in self.element_translator.keys():
				self.element_translator[node_id] = synonyms
		return None


	def determine_node(self, node_id, label, sentences, syns, cName, ECs):
		'''
		Fill the corresponding dictionary depending 
		on the type of node
		'''
		if 'Compound' in label:
			if node_id not in self.compounds_dictio.keys():
				self.compounds_dictio[node_id] = {'sentences': set(sentences), 'compoundNames': set(cName)}
			else:
				for sentence in sentences:
					self.compounds_dictio[node_id]['sentences'].add(sentence)
				for name in cName:
					self.compounds_dictio[node_id]['compoundNames'].add(name)
		if 'Protein' in label:
			if node_id not in self.proteins_dictio.keys():
				self.proteins_dictio[node_id] = {'sentences': set(sentences), 'synonyms': set(syns), 'ECs': set(ECs)}
			else:
				for sentence in sentences:
					self.proteins_dictio[node_id]['sentences'].add(sentence)
				for syn in syns:
					self.proteins_dictio[node_id]['synonyms'].add(syn)
				for EC in ECs:
					self.proteins_dictio[node_id]['ECs'].add(EC)
		return None


	def compare_synonyms(self, synonymA, synonymB, compareDB):
		'''
		Compare two synonyms lists, output the maximum 
		score between the elements of the two lists.
		There are two typs of comparison, if the enzyme without the
		action are similar, then a high similarity will be outputed.

		If the compound is not the same, we output the mean between the
		compound similarity and the enzyme similarity 

		If one of the comparisons total, automatically outputs 1.0
		'''
		list_of_scores = []
		for synA in synonymA:
			ec_tester = ''.join([i for i in synA.replace('.','') if not i.isdigit()])
			if (synA != 'More') and (ec_tester != ''):
				for synB in synonymB:
					ec_tester = ''.join([i for i in synB.replace('.','') if not i.isdigit()])
					if (synB != 'More') and (ec_tester != ''):
						synA_mod = compareDB.delete_enzyme_action(synA)
						synB_mod = compareDB.delete_enzyme_action(synB)
						score1 = calculateSimilarity(compareDB.clean(synA_mod), compareDB.clean(synB_mod))
						if score1 > 0.7:
							score = calculateSimilarity(compareDB.clean(synA), compareDB.clean(synB))
						else:
							score = calculateSimilarity(compareDB.clean(synA), compareDB.clean(synB))
							score = (score + score1)/2
						if score == 1.0:
							if (len(compareDB.clean(synA_mod)) <= 2) or (len(compareDB.clean(synB_mod)) <= 2):
								score = score - 0.3
							else:
								return score
						else:
							if (len(compareDB.clean(synA_mod)) <= 4) or (len(compareDB.clean(synB_mod)) <= 4):
								score = score - 0.3
						list_of_scores.append(score)

		#Obtain the maximum similarity between the two entities
		similarity = sorted(list_of_scores, reverse = True)[0]
		return similarity


	def similarity_matrix(self, compareDB, metab_blacklist):
		'''
		Create a similarity matrix based on the names of all the compounds
		Uses the dictionaries from the object,
		must have run neo4jtodictionary first
		'''
		modified_metab_blacklist = [compareDB.clean(i) for i in metab_blacklist]
		protein_similarity_matrix = []
		compound_similarity_matrix = []
		compound_names = []; protein_names = []
		for elem in self.compounds_dictio.keys():
			nameA = self.element_translator[elem][0]
			row = []
			if compareDB.clean(nameA) not in modified_metab_blacklist:
				for comp_elem in self.compounds_dictio.keys():
					nameB = self.element_translator[comp_elem][0]
					if compareDB.clean(nameB) not in modified_metab_blacklist:
						similarity = calculateSimilarity(compareDB.clean(nameA), compareDB.clean(nameB))
						#Penalization if the word is very small
						if similarity != 1.0:
							if (len(compareDB.clean(nameA).replace("ate", "")) <= 4) or (len(compareDB.clean(nameB).replace("ate", "")) <= 4):
								similarity = similarity - 0.3
						else:
							if (len(compareDB.clean(nameA)) <= 2) or (len(compareDB.clean(nameB)) <= 2):
								similarity = similarity - 0.3
						row.append(similarity)
			if len(row) > 0:
				compound_similarity_matrix.append(row)
				compound_names.append(elem)
		for elem in self.proteins_dictio.keys():
			synonymsA = self.element_translator[elem]
			row = []
			for elem2 in self.proteins_dictio.keys():
				synonymsB = self.element_translator[elem2]
				similarity = self.compare_synonyms(synonymsA, synonymsB, compareDB)
				row.append(similarity)
			protein_similarity_matrix.append(row)
			protein_names.append(elem)
		return protein_similarity_matrix, compound_similarity_matrix, compound_names, protein_names


	def clustering_similarities(self, matrix, eps=0.5):
		'''
		Using clustering to group the nodes
		'''
		nArray = np.array(matrix)
		db = DBSCAN(eps=eps, min_samples=1).fit(nArray)
		core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
		core_samples_mask[db.core_sample_indices_] = True
		labels = db.labels_

		# Number of clusters in labels, ignoring noise if present.
		n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)

		print('Estimated number of clusters: %d' % n_clusters_)
		try:
			print("Silhouette Coefficient: %0.3f"
		      % metrics.silhouette_score(nArray, labels))
		except ValueError:
			pass
		return labels


	def match_elem_cluster(self, cluster, element):
		'''
		Match an element with its cluster and store
		in a dictionary
		'''
		element2cluster = {}
		for n in range(0,len(element)):
			element2cluster[element[n]] = cluster[n]
		return element2cluster


	def is_protein(self, element):
		'''
		True if the element is in the protein dictionary
		'''
		if element in self.proteins_dictio.keys():
			return True
		else:
			return False


	def add_node(self, node_id, compound2cluster, protein2cluster):
		'''
		Check if the given node is a protein or a component, 
		determine if it has been already added and, if so,
		just include the properties in the graph
		'''
		is_protein = self.is_protein(node_id)
		if is_protein == True:
			cluster = protein2cluster[node_id]
			sent = list(self.proteins_dictio[node_id]['sentences'])
			syns = list(self.proteins_dictio[node_id]['synonyms'])
			ECs = list(self.proteins_dictio[node_id]['ECs'])
			node_type = 'Protein'
			if cluster not in self.prot_processed_clusters.keys():
				self.neo4j.create_second_level_enzyme(node_id, syns, sent, ECs, cluster)
				self.prot_processed_clusters[cluster] = [node_id]
			else:
				previous_nodes = self.prot_processed_clusters[cluster]
				self.neo4j.modify_second_level_enzyme(previous_nodes[0], node_id, sent, syns, ECs)
				if node_id not in self.prot_processed_clusters[cluster]:
					self.prot_processed_clusters[cluster].append(node_id)

		elif is_protein == False:
			sent = list(self.compounds_dictio[node_id]['sentences'])
			name = list(self.compounds_dictio[node_id]['compoundNames'])
			node_type = 'Compound'
			try:
				cluster = compound2cluster[node_id]
				if cluster not in self.comp_processed_clusters.keys():
					self.neo4j.create_second_level_compound(node_id, name, sent, cluster)
					self.comp_processed_clusters[cluster] = [node_id]
				else:
					previous_nodes = self.comp_processed_clusters[cluster]
					self.neo4j.modify_second_level_compound(previous_nodes[0], node_id, name, sent)
					if node_id not in self.comp_processed_clusters[cluster]:
						self.comp_processed_clusters[cluster].append(node_id)
			except KeyError:
				node_type = None
		return node_type

	def create_new_level_graph(self, compound2cluster, protein2cluster):
		'''
		Iterate through the relationships and create/fill nodes
		and relationships 
		'''
		self.comp_processed_clusters = {}; self.prot_processed_clusters = {}

		for k, v in self.rels_dictio.iteritems():
			k_type = self.add_node(k, compound2cluster, protein2cluster)
			for nodeB in v.keys():
				#add node
				node_type = self.add_node(nodeB, compound2cluster, protein2cluster)
				#add relationship -> need to consider to modify an existing relationship
				sentences = self.rels_dictio[k][nodeB]
				if not (k_type == 'Compound' and node_type == 'Compound'):
					if not (k == nodeB):
						self.neo4j.create_compressed_relationship(k, nodeB, k_type, node_type, sentences)
		return None



if __name__ == '__main__':
	
	chebi = chebi_from_string()
	chebi.chebi_connect()
	uniprot = uniprot_queries('Escherichia coli', '9606')

	gd = graph_database(chebi, uniprot, '', 'Escherichia coli', 'neo4j', 'salva') 
	gd.connect()

	compareDB = compare_string_to_db(gd)

	pattern = "MATCH (n)-[r]->(y) WHERE ('database' IN keys(n) AND 'textname' IN keys(n))\
	     AND ('database' IN keys(y) AND 'textname' IN keys(y))\
	      RETURN n.id as nID, n.sentences as nSentences, y.id as yID, y.sentences as ySentences,\
	      r.sentences as rSentences, n.ECs as nECs, y.ECs as yECs, n.compoundName as nCompound,\
	      y.compoundName as yCompound, n.synonyms as nsyn, y.synonyms as ysyn, labels(n) as nlabel, labels(y) as ylabel,\
	      y.uniprotProteinNames as yuniprotProteinNames, n.uniprotProteinNames as nuniprotProteinNames"

	# pattern = "MATCH (n)-[r]-(y)-[s]-(z) WHERE n.id IN ['CHEBI:21547', 'CHEBI:16953'] AND n<>z\
	# 	  RETURN n.id as nID, n.sentences as nSentences, y.id as yID, y.sentences as ySentences,\
	#       r.sentences as rSentences, n.ECs as nECs, y.ECs as yECs, n.compoundName as nCompound,\
	#       y.compoundName as yCompound, n.synonyms as nsyn, y.synonyms as ysyn, labels(n) as nlabel, labels(y) as ylabel,\
	#       y.uniprotProteinNames as yuniprotProteinNames, n.uniprotProteinNames as nuniprotProteinNames"

	metab_blacklist = ['water', 'H2O', 'proton', 'H+', 'CO2', 'carbon dioxide',\
		'CMP', 'AMP', 'PPi', 'pirophosphate', 'diphosphate', 'polyphosphate', 'triphosphate',\
		'Dicarboxylic acid dianion', 'NAD', 'NADP', 'NADH', 'NADPH', 'O2', 'Pi', 'ADP', 'ATP',\
		'3-5ADP', "adenosine 5'-monophosphate", 'NH4+', 'ammonium', 'GTP', 'GDP', 'CTP',\
		 'purine', 'pirimidine', 'guanosine phosphate', 'magnesium(2+)', 'ethanol', 'hydrogen']

	#metab_blacklist = []

	compress = compress_graph(gd)
	compress.neo4j_to_dictionary(pattern)

	protein_similarity_matrix, compound_similarity_matrix, compound_names, protein_names = compress.similarity_matrix(compareDB, metab_blacklist)

	labels = compress.clustering_similarities(compound_similarity_matrix, eps=1.0)

	compound2cluster = compress.match_elem_cluster(labels, compound_names)

	labels = compress.clustering_similarities(protein_similarity_matrix, eps=1.0)
	protein2cluster = compress.match_elem_cluster(labels, protein_names)

	compress.create_new_level_graph(compound2cluster, protein2cluster)