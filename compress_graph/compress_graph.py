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
		result = gd.extract_pattern(pattern)
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

	# pattern = "MATCH (n)-[r]-(y)-[s]-(z) WHERE n.id IN ['CHEBI:21547', 'CHEBI:16953', 'CHEBI:21647', 'CHEBI:16567', 'CHEBI:23628', 'CHEBI:17549', 'CHEBI:12109', 'CHEBI:724125', 'CHEBI:21601', 'CHEBI:15748', 'CHEBI:85313', 'CHEBI:18391', 'CHEBI:16737', 'CHEBI:28393', 'CHEBI:22102', 'CHEBI:25001', 'CHEBI:17268', 'CHEBI:10983', 'CHEBI:11805', 'CHEBI:37054', 'CHEBI:16724', 'CHEBI:18072', 'CHEBI:18064', 'CHEBI:18050', 'CHEBI:62423', 'CHEBI:27823', 'CHEBI:27978', 'CHEBI:15621', 'CHEBI:17612', 'CHEBI:41941', 'CHEBI:28497', 'CHEBI:16313', 'CHEBI:18112', 'CHEBI:91019', 'CHEBI:38622', 'CHEBI:50428', 'CHEBI:44303', 'CHEBI:30662', 'CHEBI:18186', 'CHEBI:50607', 'CHEBI:18332', 'CHEBI:15966', 'CHEBI:16810', 'CHEBI:29985', 'CHEBI:30921', 'CHEBI:15903', 'CHEBI:17009', 'CHEBI:17234', 'CHEBI:42548', 'CHEBI:29062', 'CHEBI:86262', 'CHEBI:17196', 'CHEBI:21553', 'CHEBI:17786', 'CHEBI:40410', 'CHEBI:21565', 'CHEBI:11951', 'CHEBI:63815', 'CHEBI:35932', 'CHEBI:16414', 'CHEBI:29995', 'CHEBI:29990', 'CHEBI:17053', 'CHEBI:32508', 'CHEBI:15428', 'CHEBI:16180', 'CHEBI:17310', 'CHEBI:17405', 'CHEBI:33409', 'CHEBI:25094', 'CHEBI:16828', 'CHEBI:17191', 'CHEBI:15603', 'CHEBI:43433', 'CHEBI:17368', 'CHEBI:17016', 'CHEBI:49031', 'CHEBI:30854'] AND z.id IN ['CHEBI:21547', 'CHEBI:16953', 'CHEBI:21647', 'CHEBI:16567', 'CHEBI:23628', 'CHEBI:17549', 'CHEBI:12109', 'CHEBI:724125', 'CHEBI:21601', 'CHEBI:15748', 'CHEBI:85313', 'CHEBI:18391', 'CHEBI:16737', 'CHEBI:28393', 'CHEBI:22102', 'CHEBI:25001', 'CHEBI:17268', 'CHEBI:10983', 'CHEBI:11805', 'CHEBI:37054', 'CHEBI:16724', 'CHEBI:18072', 'CHEBI:18064', 'CHEBI:18050', 'CHEBI:62423', 'CHEBI:27823', 'CHEBI:27978', 'CHEBI:15621', 'CHEBI:17612', 'CHEBI:41941', 'CHEBI:28497', 'CHEBI:16313', 'CHEBI:18112', 'CHEBI:91019', 'CHEBI:38622', 'CHEBI:50428', 'CHEBI:44303', 'CHEBI:30662', 'CHEBI:18186', 'CHEBI:50607', 'CHEBI:18332', 'CHEBI:15966', 'CHEBI:16810', 'CHEBI:29985', 'CHEBI:30921', 'CHEBI:15903', 'CHEBI:17009', 'CHEBI:17234', 'CHEBI:42548', 'CHEBI:29062', 'CHEBI:86262', 'CHEBI:17196', 'CHEBI:21553', 'CHEBI:17786', 'CHEBI:40410', 'CHEBI:21565', 'CHEBI:11951', 'CHEBI:63815', 'CHEBI:35932', 'CHEBI:16414', 'CHEBI:29995', 'CHEBI:29990', 'CHEBI:17053', 'CHEBI:32508', 'CHEBI:15428', 'CHEBI:16180', 'CHEBI:17310', 'CHEBI:17405', 'CHEBI:33409', 'CHEBI:25094', 'CHEBI:16828', 'CHEBI:17191', 'CHEBI:15603', 'CHEBI:43433', 'CHEBI:17368', 'CHEBI:17016', 'CHEBI:49031', 'CHEBI:30854'] AND n<>z\
	# 	  RETURN n.id as nID, n.sentences as nSentences, y.id as yID, y.sentences as ySentences,\
	#       r.sentences as rSentences, n.ECs as nECs, y.ECs as yECs, n.compoundName as nCompound,\
	#       y.compoundName as yCompound, n.synonyms as nsyn, y.synonyms as ysyn, labels(n) as nlabel, labels(y) as ylabel,\
	#       y.uniprotProteinNames as yuniprotProteinNames, n.uniprotProteinNames as nuniprotProteinNames"

	metab_blacklist = ['water', 'H2O', 'proton', 'H+', 'CO2', 'carbon dioxide',\
		'CMP', 'AMP', 'PPi', 'pirophosphate', 'diphosphate', 'polyphosphate', 'triphosphate',\
		'Dicarboxylic acid dianion', 'NAD', 'NADP', 'NADH', 'NADPH', 'O2', 'Pi', 'ADP', 'ATP',\
		'3-5ADP', "adenosine 5'-monophosphate", 'NH4+', 'ammonium', 'GTP', 'GDP', 'CTP',\
		 'purine', 'pirimidine', 'guanosine phosphate', 'magnesium(2+)', 'ethanol', 'hydrogen']

	metab_blacklist = []

	compress = compress_graph(gd)
	compress.neo4j_to_dictionary(pattern)

	protein_similarity_matrix, compound_similarity_matrix, compound_names, protein_names = compress.similarity_matrix(compareDB, metab_blacklist)

#	for elem in compress.proteins_dictio.values():
#		print elem['synonyms']

#	print len(protein_similarity_matrix), len(protein_similarity_matrix[0])

#	print compound_similarity_matrix
#	for x in [0.5, 0.7, 0.8, 0.9]:
	labels = compress.clustering_similarities(compound_similarity_matrix, eps=1.0)

#	print '\n\n\n\n'
#	print labels

	compound2cluster = compress.match_elem_cluster(labels, compound_names)

#	for x in [0.5, 0.7, 0.8, 0.9]:
	labels = compress.clustering_similarities(protein_similarity_matrix, eps=1.0)
	protein2cluster = compress.match_elem_cluster(labels, protein_names)

#	print '\n\n\n\n'
#	print labels
#	print protein2cluster

#	compound2cluster = {u'CHEBI:17790': 0, u'CHEBI:58476': 1, u'CHEBI:29991': 2, u'CHEBI:15603': 19, u'CHEBI:15525': 7, u'CHEBI:15687': 4, u'CHEBI:36655': 5, u'CHEBI:16414': 22, u'CHEBI:35366': 6, u'CHEBI:33019': 17, u'CHEBI:16704': 45, u'CHEBI:17381': 8, u'CHEBI:17668': 9, u'CHEBI:13850': 10, u'CHEBI:16708': 11, u'CHEBI:17418': 12, u'CHEBI:29032': 33, u'CHEBI:15635': 14, u'CHEBI:29985': 15, u'CHEBI:17326': 9, u'CHEBI:47191': 16, u'CHEBI:16838': 17, u'CHEBI:24040': 18, u'CHEBI:17094': 3, u'CHEBI:15346': 20, u'CHEBI:15343': 21, u'CHEBI:32988': 23, u'CHEBI:24898': 19, u'CHEBI:30089': 24, u'CHEBI:27266': 22, u'CHEBI:40574': 25, u'CHEBI:67016': 14, u'CHEBI:15351': 26, u'CHEBI:28644': 27, u'CHEBI:16842': 28, u'CHEBI:38083': 29, u'CHEBI:15428': 30, u'CHEBI:16530': 1, u'CHEBI:16335': 11, u'CHEBI:26806': 31, u'CHEBI:15980': 32, u'CHEBI:7916': 33, u'CHEBI:11812': 1, u'CHEBI:14737': 32, u'CHEBI:43103': 34, u'CHEBI:16856': 35, u'CHEBI:37075': 9, u'CHEBI:29805': 5, u'CHEBI:18697': 32, u'CHEBI:17984': 36, u'CHEBI:11561': 3, u'CHEBI:17625': 37, u'CHEBI:16753': 38, u'CHEBI:58450': 16, u'CHEBI:1989': 14, u'CHEBI:23019': 39, u'CHEBI:15769': 40, u'CHEBI:15763': 41, u'CHEBI:20612': 14, u'CHEBI:51277': 42, u'CHEBI:33704': 43, u'CHEBI:16862': 9, u'CHEBI:30849': 44, u'CHEBI:33709': 43, u'CHEBI:53152': 7, u'CHEBI:17188': 9, u'CHEBI:7754': 13, u'CHEBI:16454': 33, u'CHEBI:28971': 46, u'CHEBI:18359': 47, u'CHEBI:27594': 48, u'CHEBI:73651': 49, u'CHEBI:48942': 12, u'CHEBI:32816': 50, u'CHEBI:12109': 51, u'CHEBI:17196': 52, u'CHEBI:15570': 22, u'CHEBI:16189': 53, u'CHEBI:17191': 19, u'CHEBI:17053': 2, u'CHEBI:28340': 54, u'CHEBI:17295': 34, u'CHEBI:18347': 19, u'CHEBI:16467': 55, u'CHEBI:17561': 56, u'CHEBI:27585': 57, u'CHEBI:17115': 58, u'CHEBI:18021': 59, u'CHEBI:15740': 60, u'CHEBI:33838': 61, u'CHEBI:16000': 62, u'CHEBI:16139': 47, u'CHEBI:16004': 24, u'CHEBI:64552': 63, u'CHEBI:13534': 47, u'CHEBI:27575': 16, u'CHEBI:15468': 64, u'CHEBI:17972': 9, u'CHEBI:12071': 14, u'CHEBI:12936': 65, u'CHEBI:57845': 66, u'CHEBI:16646': 67, u'CHEBI:16763': 27, u'CHEBI:17360': 68, u'CHEBI:18036': 17, u'CHEBI:18409': 1, u'CHEBI:16018': 47, u'CHEBI:57774': 69, u'CHEBI:16398': 70, u'CHEBI:17141': 56, u'CHEBI:15905': 71, u'CHEBI:16958': 22, u'CHEBI:11851': 1, u'CHEBI:57308': 66, u'CHEBI:15361': 50, u'CHEBI:25115': 72, u'CHEBI:37240': 16, u'CHEBI:73313': 73}
#	protein2cluster = {u'P22259': 0, u'P31057': 1, u'P08373': 2, u'P0AF03': 3, u'P11868': 4, u'P42597': 5, u'P37330': 6, u'P07623': 3, u'P0A7D7': 7, u'P60560': 8, u'P05041': 9, u'P0AAI5': 10, u'P0A6X1': 11, u'P27550': 12, u'P0A6A8': 13, u'P62615': 14, u'P0AB80': 15, u'P0A6I6': 16, u'P0A9B6': 17, u'P08839': 18, u'P0A6I3': 19, u'P13445': 20, u'P0A6I9': 21, u'P42588': 22, u'P0A9G6': 57, u'P0A7M2': 24, u'P06983': 25, u'P32670': 18, u'P06709': 26, u'P0ABS8': 27, u'P21170': 28, u'P0A790': 29, u'P08401': 30, u'P21179': 31, u'P04693': 15, u'P0A953': 10, u'P0A884': 32, u'P0AB77': 38, u'P23845': 3, u'P31663': 19, u'P0A800': 27, u'Q59385': 35, u'P69902': 36, u'P0ABQ0': 37, u'P25437': 33, u'P07003': 39, u'P00895': 34, u'P75823': 40, u'P0A9M8': 74, u'P0AB65': 41, u'P16692': 70, u'P04395': 43, u'P04995': 44, u'P27248': 45, u'P0AG11': 27, u'P0A8V2': 27, u'P75691': 46, u'P76539': 47, u'P0A9Q7': 46, u'P52612': 35, u'P13009': 48, u'P17445': 49, u'P24224': 50, u'P80668': 51, u'P26281': 52, u'P14377': 30, u'P21515': 44, u'P00634': 41, u'P42630': 53, u'P0A825': 45, u'P0A7F6': 55, u'P00926': 53, u'P0A6L4': 23, u'P0A8N5': 58, u'P00864': 0, u'P0A6T3': 59, u'P0A9V8': 60, u'P0AFG8': 39, u'Q46893': 56, u'P60720': 36, u'P11071': 57, u'P07464': 38, u'P52097': 61, u'P09126': 62, u'P0ADF8': 63, u'P13029': 64, u'P37613': 65, u'P0AG30': 66, u'P0A752': 3, u'P52643': 67, u'P09323': 18, u'P76558': 68, u'P0AGJ7': 45, u'P0A910': 69, u'P16690': 70, u'P00936': 71, u'P0ACB2': 42, u'P04968': 53, u'P45568': 72, u'P69786': 18, u'P0ADA1': 73, u'P0A6R0': 10, u'P31441': 54, u'P05100': 43, u'P23538': 75, u'P30870': 22, u'P00861': 55, u'P0A7A9': 5, u'P0ABE9': 76, u'P05791': 77, u'P05793': 78, u'P0ACC7': 79, u'P27306': 80, u'P06968': 81, u'P42604': 82, u'P23908': 89, u'P0AFR2': 84, u'P0AD42': 5, u'P0A8K1': 85, u'P05055': 86, u'P09546': 87, u'P0A9J4': 88, u'P00509': 22, u'P37623': 50, u'P18843': 83, u'P0AFL6': 90}

	compress.create_new_level_graph(compound2cluster, protein2cluster)