#!/usr/bin/env python

'''
Use the TF information curated from Marta 
in order to create the relationship TF-gene 
in our database.
'''

import time, sys, re, urllib2
from graph_database import graph_database
from chebi_from_string import chebi_from_string
from brenda_annotation import brenda_annotation
from uniprot_queries import uniprot_queries
from intact_database import intact_database


class TF_data():
	'''
	Class to process the text file that state the 
	relationship between the TF and the regulated gene
	'''

	def __init__(self, neo4j, uniprot):

		self.neo4j = neo4j
		self.uniprot = uniprot
		self.organism = 'Homo sapiens'


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


	def add_relationship(self, protein1, protein2):
		'''
		Add the relationship to the neo4j database,
		including the properties
		'''
		species = set([self.organism.encode('utf8')])
		is_rel, r = self.neo4j.check_relationship(protein1, protein2, 'TF_regulation')
		if is_rel == False:
			query = "MATCH (n), (y) WHERE n.id = '%s' AND y.id = '%s'\
			 CREATE (n)-[r:TF_regulation {species:%s}]->(y)"%\
			 (protein1, protein2, list(species))
			self.neo4j.session.run(query)
		return None

	def process_files(self):
		'''
		Process the TF gene relationship data file
		'''
		file = open('data/TFdata/proximal_and_distal_tf.txt').read()
		for line in file.split('\n')[:-1]:
			protein1 = line.split('\t')[0]
			protein2 = line.split('\t')[1]
			a = self.uniprot.query_id(protein1)
			b = self.uniprot.query_id(protein2)
			if (a[0] != '' and b[0] != ''):
				self.add_protein(a[0], a[1], a[2], a[3])
				self.add_protein(b[0], b[1], b[2], b[3])
				self.add_relationship(a[0], b[0])
		return None


if __name__ == '__main__':

	start = time.time()

	email = 'salcagal@alumni.uv.es'; brendapass = 'salvacasani91'

	chebi = chebi_from_string()
	chebi.chebi_connect()
	brenda = brenda_annotation(email,brendapass)
	brenda.access_protocol()
	neo4j = graph_database(chebi, brenda,'Homo sapiens','neo4j','salva')
	neo4j.connect()
	uniprot = uniprot_queries('Homo sapiens', '9606')

	tf = TF_data(neo4j, uniprot)
	tf.process_files()

	end = time.time()
	print 'The script took: ' + str(end - start)