#!/usr/bin/env python

'''
Add the pazar database to the Neo4j graph database
'''

import glob, sys, time, re
from uniprot_queries import uniprot_queries
from graph_database import graph_database
from chebi_from_string import chebi_from_string
from brenda_annotation import brenda_annotation

class pazar_database():

	def __init__(self, organism, neo4j, uniprot_query):

		self.organism = organism
		self.rel_dict = {}
		self.neo4j = neo4j
		self.uniprot_query = uniprot_query
		self.information_dictio = {}


	def get_data(self):
		'''
		Obtain all the files to be processed.
		Carefull: The directories are harcodded.
		Better not to move files/folders
		'''
		directory = './data/pazar_data/'
		files = glob.glob(directory + '*')
		return files


	def fill_info_dict(self, protein, name, genes, proteins):
		'''
		Insert the information from the database
		to get all the info from a protein
		'''
		if protein not in self.information_dictio.keys():
			protList = []; geneList = []
			for prot in proteins.split('('):
				protList.append(re.sub('\)', '', prot).encode('utf8'))
			for gene in genes.split():
				geneList.append(re.sub(':', '', gene).encode('utf8'))
			self.information_dictio[protein] = {'Entry':name, 'GeneNames':geneList, 'ProteinNames':protList}
		return None


	def process_files(self):
		'''
		Obtain the TF and target gene according to 
		the pazar DB files
		'''
		files = self.get_data()
		#files = ['/home/salva/extract_path_neo4j/data/pazar_data/pazar_ABS_20120522.csv']
		for file in files:
			lines = [line for line in open(file) if self.organism in line]
			for line in lines:
				line = line.split('\t')
				if line[8] == self.organism:
					TF = line[2]; gene = line[4]
					pmid = line[-2]; method = line[-1]
					TFpage = self.uniprot_query.mapping_id(TF, 'ID')
					TF = TFpage[0]; entry = TFpage[1]; genes = TFpage[2]; proteins = TFpage[3]
					self.fill_info_dict(TF, entry, genes, proteins)
					genePage = self.uniprot_query.mapping_id(gene, 'ENSEMBL_ID')
					gene = genePage[0]; entry = genePage[1]; genes = genePage[2]; proteins = genePage[3]
					self.fill_info_dict(gene, entry, genes, proteins)
					if (TF != '' and gene != ''):
						if TF not in self.rel_dict.keys():
							self.rel_dict[TF] = {}
						if gene in self.rel_dict[TF].keys():
							self.rel_dict[TF][gene]['method'].add(method)
							self.rel_dict[TF][gene]['pmid'].add(pmid)
						else:
							self.rel_dict[TF][gene] = {'method':set([method]), 'pmid':set([pmid])}
		return None


	def add_protein(self, protein, TF):
		'''
		Check if a protein exists and add it with the 
		required properties
		'''
		is_node = self.neo4j.check_protein(protein)
		if is_node == False:
			entry = self.information_dictio[protein]['Entry']
			genes = self.information_dictio[protein]['GeneNames']
			proteins = self.information_dictio[protein]['ProteinNames']
			query = "CREATE (n:Protein {uniprotID:'%s', id:'%s', uniProtEntryName: '%s',\
			 uniprotGenesNames: %s, uniprotProteinNames: %s, specie: %s})"%\
			 (protein, protein, entry, genes, proteins, self.organism)
			self.neo4j.session.run(query)
		if TF == True:
			query = "MATCH (n) WHERE n.id = '%s' SET n.TF = True"%\
			(protein)
			self.neo4j.session.run(query)
		return None


	def add_relationship(self, nodeA, nodeB, pmid, method):
		'''
		Check if a relationship exists and add it if necessary
		fill the proper properties
		'''
		species = set([self.organism])
		is_rel, r = self.neo4j.check_relationship(nodeA, nodeB, 'Pazar_relationship')
		if is_rel == True:
			for result in r:
				method.update(a.encode('utf-8') for a in result['r']['method'])
				pmid.update(a.encode('utf-8') for a in result['r']['pmid'])
				species.update(a.encode('utf-8') for a in result['r']['species'])
			query = 'MATCH (n)-[r]-(y) WHERE n.id = "%s" AND y.id ="%s" SET \
			 r.method = %s, r.pmid = %s, r.species = %s'%\
			 (nodeA, nodeB, list(method), list(pmid), list(species))
			self.neo4j.session.run(query)
		else:
			query = "MATCH (n), (y) WHERE n.id = '%s' AND y.id = '%s'\
			 CREATE (n)-[r:Pazar_relationship {species:%s, pmid:%s, method:%s}]->(y)"%\
			 (nodeA, nodeB, list(species), list(pmid), list(method))
			self.neo4j.session.run(query)
		return None


	def fill_database(self):
		'''
		Fill the Neo4j database using the information from the 
		dictionary constructed from the pazar data
		'''
		for k, v in self.rel_dict.iteritems():
			self.add_protein(k, True)
			for gene in v.keys():
				method = self.rel_dict[k][gene]['method']
				pmid = self.rel_dict[k][gene]['pmid']
				self.add_protein(gene, False)
				self.add_relationship(k, gene, pmid, method)
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
	uniprot_query = uniprot_queries('Homo sapiens', '9606')

	pazar = pazar_database('Homo sapiens', neo4j, uniprot_query)
	pazar.process_files()
	
	pazar.fill_database()

	end = time.time()
	print 'The script took: ' + str(end - start)