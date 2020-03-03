#!/usr/bin/env python

'''
This script takes the txt file from intact that describes the interactions
and retrieves the identifiers and the confidence of the interaction,
which allows the user to filter the confidence
'''


import time, sys, re
from graph_database import graph_database
from chebi_from_string import chebi_from_string
from brenda_annotation import brenda_annotation
from uniprot_queries import uniprot_queries


class intact_database():
	
	'''
	Retrieve protein-protein interactions from the 
	IntAct database for a given organism
	'''

	def __init__(self, neo4j, organism, uniprot):

		self.neo4j = neo4j
		self.organism = organism
		self.uniprot = uniprot


	def process_intact(self, confidence_thr = 0.5):
		'''
		Process the intact database file and extract the
		relevant features/interactions
		'''
		file = open('data/intact.txt').read().split('\n')
		print file[0]
		for line in file[1:]:
			line = line.split('\t')
			try:
				species = line[28].split('|')[0].split(':')[1].split('(')[0]
				if self.organism == species:
					num = float(line[14].split('intact-miscore:')[1])
					if num >= confidence_thr:
						exp_tech = line[6].split('(')[1].split(')')[0]
						int_type = line[11].split('(')[1].split(')')[0]
						protein1 = line[0].split(':')[1]
						protein2 = line[1].split(':')[1]
						a = self.uniprot.uni_search(protein1)
						b = self.uniprot.uni_search(protein2)
						if (a[0] == protein1 and b[0] == protein2):
							self.add_protein(protein1, a[1], a[2], a[3])
							self.add_protein(protein2, b[1], b[2], b[3])
							self.add_relationship(protein1, protein2, exp_tech, int_type, num)
			except IndexError:
				pass
		return None


	def add_protein(self, uniprotID, uniprotName, uniGenes, uniProteins):
		'''
		Add the protein to the neo4j database
		'''
		is_node = self.neo4j.check_protein(uniprotID)
		if is_node == False:
			protList = []; geneList = []
			for prot in uniProteins.split('('):
				protList.append(re.sub('\)', '', prot).encode('utf8'))
			for gene in uniGenes.split():
				geneList.append(re.sub(':', '', gene).encode('utf8'))
			query = "CREATE (n:Protein {uniprotID:'%s', id:'%s', uniProtEntryName: '%s',\
			 uniprotGenesNames: %s, uniprotProteinNames: %s, specie: %s})"%\
			 (uniprotID, uniprotID, uniprotName, geneList, protList, self.organism)
			self.neo4j.session.run(query)
		return None


	def add_relationship(self, protein1, protein2, technique, int_type, confidence):
		'''
		Add the relationship to the neo4j database,
		including the properties
		'''
		species = set([self.organism.encode('utf8')]); experimental_score = set([confidence])
		int_type = set([int_type.encode('utf8')]); technique = set([technique.encode('utf8')])
		is_rel, r = self.neo4j.check_relationship(protein1, protein2, 'IntAct_interaction')
		if is_rel == True:
			for result in r:
				experimental_score.update(result['r']['experimental_score'])
				species.update(a.encode('utf-8') for a in result['r']['species'])
				int_type.update(a.encode('utf-8') for a in result['r']['interaction_type'])
				technique.update(a.encode('utf-8') for a in result['r']['detection_method'])
			query = 'MATCH (n)-[r]-(y) WHERE n.id = "%s" AND y.id ="%s" SET \
			 r.experimental_score = %s, r.species = %s,\
			  r.interaction_method = %s, r.detection_method = %s'%\
			 (protein1, protein2, list(experimental_score), list(species),\
			  list(int_type), list(technique))
			self.neo4j.session.run(query)
		else:
			query = "MATCH (n), (y) WHERE n.id = '%s' AND y.id = '%s'\
			 CREATE (n)-[r:IntAct_interaction {species:%s, experimental_score:%s,\
			  interaction_type:%s, detection_method:%s}]->(y)"%\
			 (protein1, protein2, list(species), list(experimental_score),\
			 list(int_type), list(technique))
			self.neo4j.session.run(query)
		return None


if __name__ == '__main__':

	start = time.time()

	email = 'salcagal@alumni.uv.es'; brendapass = 'salvacasani91'
	file = 'data/9606.protein.links.detailed.v10.5.txt'

	chebi = chebi_from_string()
	chebi.chebi_connect()
	brenda = brenda_annotation(email,brendapass)
	brenda.access_protocol()
	neo4j = graph_database(chebi, brenda,'Homo sapiens','neo4j','salva')
	neo4j.connect()
	uniprot = uniprot_queries('Homo sapiens', '9606')

	intact = intact_database(neo4j, '9606', uniprot)
	intact.process_intact()


	end = time.time()
	print 'The script took: ' + str(end - start)