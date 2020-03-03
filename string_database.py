#!/usr/bin/env python

'''
Add the pazar database to the Neo4j graph database

Here we still need to add some filters to choose the 
correct files to process all this data based on the organism
'''

import glob, sys, time, re
from uniprot_queries import uniprot_queries
from graph_database import graph_database
from chebi_from_string import chebi_from_string
from brenda_annotation import brenda_annotation


class string_database():

	def __init__(self, neo4j, organism, uniprot):

		self.neo4j = neo4j
		self.organism = organism
		self.uniprot = uniprot


	def process_dbFile(self, dbFile, filteringScore):
		'''
		Use the file with the interactions from the string database
		and filter the interesting relationships
		'''
		with open(dbFile) as file:
			for line in file.read().split('\n')[1:-1]:
				line = line.split()
				experimental_score = int(line[6])
				if experimental_score > filteringScore:
					#protein1 = line[0].split('.')[1]
					#protein2 = line[1].split('.')[1]
					gene1Page = self.uniprot.mapping_id(line[0], 'STRING_ID')
					gene2Page = self.uniprot.mapping_id(line[1], 'STRING_ID')
					if (gene1Page[0] != '' and gene2Page[0] != ''):
						
						self.add_protein(gene1Page[0], gene1Page[1], gene1Page[2], gene1Page[3])
						self.add_protein(gene2Page[0], gene2Page[1], gene2Page[2], gene2Page[3])
						self.add_relationship(gene1Page[0], gene2Page[0], experimental_score)
		return None




	def add_protein(self, uniprotID, uniprotName, uniGenes, uniProteins):
		'''
		Add the protein to the neo4j database as a protein instance
		'''
		is_node = self.neo4j.check_protein(uniprotID)
		if is_node == False:
			protList = []; geneList = []
			for prot in uniProteins.split('('):
				protList.append(re.sub('\)', '', prot).encode('utf8'))
			for gene in uniGenes.split():
				geneList.append(re.sub(':', '', gene).encode('utf8'))
			query = "CREATE (n:Protein {uniprotID:'%s', id:'%s', uniProtEntryName: '%s',\
			 uniprotGenesNames: %s, uniprotProteinNames: %s, specie: '%s'})"%\
			 (uniprotID, uniprotID, uniprotName, geneList, protList, self.organism)
			self.neo4j.session.run(query)
		return None


	def add_relationship(self, protein1, protein2, experimental_score):
		'''
		Add the interaction to the neo4j database as a 
		'''
		species = set([self.organism]); experimental_score = set([experimental_score])
		is_rel, r = self.neo4j.check_relationship(protein1, protein2, 'StringDB_interaction')
		if is_rel == True:
			for result in r:
				experimental_score.update(a for a in result['r']['experimental_score'])
				species.update(a.encode('utf-8') for a in result['r']['species'])
			query = 'MATCH (n)-[r]-(y) WHERE n.id = "%s" AND y.id ="%s" SET \
			 r.experimental_score = %s, r.species = %s'%\
			 (protein1, protein2, list(experimental_score), list(species))
			self.neo4j.session.run(query)
		else:
			query = "MATCH (n), (y) WHERE n.id = '%s' AND y.id = '%s'\
			 CREATE (n)-[r:StringDB_interaction {species:%s, experimental_score:%s}]->(y)"%\
			 (protein1, protein2, list(species), list(experimental_score))
			self.neo4j.session.run(query)
		return None


if __name__ == '__main__':

	start = time.time()

	email = 'salcagal@alumni.uv.es'; brendapass = 'salvacasani91'
	file = 'data/string_data/9606.protein.links.detailed.v10.5.txt'

	chebi = chebi_from_string()
	chebi.chebi_connect()
	brenda = brenda_annotation(email,brendapass)
	brenda.access_protocol()
	neo4j = graph_database(chebi, brenda,'Homo sapiens','neo4j','salva')
	neo4j.connect()
	uniprot = uniprot_queries('Homo sapiens', '9606')

	string = string_database(neo4j, 'Homo sapiens', uniprot)
	string.process_dbFile(file, 500)

	end = time.time()
	print 'The script took: ' + str(end - start)