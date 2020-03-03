#!/usr/bin/env python

'''
This script takes the txt file from intact that describes the interactions
and retrieves the identifiers and the confidence of the interaction,
which allows the user to filter the confidence
'''

import time, sys, re, urllib2
from graph_database import graph_database
from chebi_from_string import chebi_from_string
from brenda_annotation import brenda_annotation
from uniprot_queries import uniprot_queries
from intact_database import intact_database


class omnipath():
	'''
	Obtain the interactions present at the Omnipath signaling
	pathways database.
	Remember this database is restricted to human.
	'''

	def __init__(self, neo4j, uniprot):

		self.neo4j = neo4j
		self.uniprot = uniprot


	def obtain_files(self):
		'''
		Obtain the files that will be used for
		the interactions processing
		'''
		interactionFile ='http://omnipathdb.org/interactions'
		response = urllib2.urlopen(interactionFile)
		interactionFile = response.read()
		PTMfile = 'http://omnipathdb.org/ptms'
		response = urllib2.urlopen(PTMfile)
		PTMfile = response.read()
		return interactionFile, PTMfile


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
		species = set(['Homo sapiens'])
		is_rel, r = self.neo4j.check_relationship(protein1, protein2, 'Omnipath_interaction')
		if is_rel == True:
			for result in r:
				species.update(a.encode('utf-8') for a in result['r']['species'])
			query = 'MATCH (n)-[r]-(y) WHERE n.id = "%s" AND y.id ="%s" SET \
			 r.species = %s'%(protein1, protein2, list(species))
			self.neo4j.session.run(query)
		else:
			query = "MATCH (n), (y) WHERE n.id = '%s' AND y.id = '%s'\
			 CREATE (n)-[r:Omnipath_interaction {species:%s}]->(y)"%\
			 (protein1, protein2, list(species))
			self.neo4j.session.run(query)
		return None


	def add_ptm_rel(self, protein1, protein2, ptm):
		'''
		Add a PTM as a property of the relationship between
		the given proteins
		'''
		species = set(['Homo sapiens']); ptm = set([ptm])
		is_rel, r = self.neo4j.check_relationship(protein1, protein2, 'Omnipath_interaction')
		if is_rel == True:
			for result in r:
				species.update(a.encode('utf-8') for a in result['r']['species'])
				if result['r']['ptm'] != None:
					ptm.update(a.encode('utf-8') for a in result['r']['ptm'])
			query = 'MATCH (n)-[r]-(y) WHERE n.id = "%s" AND y.id ="%s" SET \
			 r.species = %s, r.ptm = %s'%(protein1, protein2, list(species), list(ptm))
			self.neo4j.session.run(query)
		else:
			query = "MATCH (n), (y) WHERE n.id = '%s' AND y.id = '%s'\
			 CREATE (n)-[r:Omnipath_interaction {species:%s, ptm:%s}]->(y)"%\
			 (protein1, protein2, list(species), list(ptm))
			self.neo4j.session.run(query)
		return None


	def process_omnipath(self):
		'''
		Go over the files in order to get the proteins that 
		interact, as well as the PTMs from the second file
		'''
		interactions, ptm = self.obtain_files()
		#Manage interactions file
		for line in interactions.split('\n')[1:-1]:
			protein1 = line.split('\t')[0]
			protein2 = line.split('\t')[1]
			a = self.uniprot.uni_search(protein1)
			b = self.uniprot.uni_search(protein2)
			if (a[0] == protein1 and b[0] == protein2):
				self.add_protein(protein1, a[1], a[2], a[3])
				self.add_protein(protein2, b[1], b[2], b[3])
				self.add_relationship(protein1, protein2)
		#Manage PTms file
		for line in ptm.split('\n')[1:-1]:
			protein1 = line.split('\t')[0]
			protein2 = line.split('\t')[1]
			modification = line.split('\t')[-1]
			a = self.uniprot.uni_search(protein1)
			b = self.uniprot.uni_search(protein2)
			if (a[0] == protein1 and b[0] == protein2):
				self.add_protein(protein1, a[1], a[2], a[3])
				self.add_protein(protein2, b[1], b[2], b[3])
				self.add_ptm_rel(protein1, protein2, modification)
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

	omnipath = omnipath(neo4j, uniprot)
	omnipath.process_omnipath()

	end = time.time()
	print 'The script took: ' + str(end - start)