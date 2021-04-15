#!/usr/bin/env python


from graph_database import graph_database
from uniprot_queries import uniprot_queries
from chebi_from_string import chebi_from_string
from urllib2 import Request, urlopen, URLError
import xml.etree.ElementTree as ET
import sys, os, re, glob, multiprocessing
from activepool import ActivePool
import httplib

class inparanoid():
	'''
	Search for orthologs of the species proteins
	'''

	def __init__(self, gd, fasta_repo, blastp_repo, inparanoid_repo):

		self.gd = gd
		self.fasta_repo = fasta_repo
		self.blastp_repo = blastp_repo
		self.inparanoid_repo = inparanoid_repo


	def research_orthologs(self, uniID, species):
		'''
		research orthologs of a specific UniProtID in the inparanoid database available on the web
        save the info in the neo4j database
		'''
		url="http://inparanoid.sbc.su.se/cgi-bin/gene_search.cgi?idtype=all;all_or_selection=all;scorelimit=0.05;rettype=xml;id=" + uniID
		request = Request(url)
		for n in range(10):
			try:
				response = urlopen(request)
				txml = response.read()
				root = ET.fromstring(txml)
				for cluster in root.iter('cluster'):
					c = cluster.attrib['nr']			#Cluster ID
					for p in cluster.iter('protein'):
						pid = p.attrib['prot_id']
						spec = p.attrib['speclong']
						if pid != uniID:
							if spec in species:
								self.gd.add_orthology_relationship(pid, uniID, c)
				break
			except URLError, e:
				print 'Warning : No ortholog for %s ' % uniID
				break
			except httplib.IncompleteRead:
				pass
		return None


	def extract_fasta(self, specie):
		'''
		Extract the fasta files for all the proteins present in the database
		'''
		uniprot_IDs = self.gd.obtain_uniIDs_from_specie(specie)
		url = 'https://www.uniprot.org/uniprot/'
		for uniID in uniprot_IDs:
			if uniID != None:
				uni_url = url + uniID + '.fasta'
				request = Request(uni_url)
				for n in range(10):
					try:
						response = urlopen(request, timeout=150).read()
						specie_name = specie.replace(' ', '_')
						with open(self.fasta_repo + '/' +  specie_name + '.fasta', 'a') as uniID_file:
							uniID_file.write(response)
						break
					except URLError, httplib.BadStatusLine:
						pass
		return None


	def run_makeBlastDb(self, file, s, pool):
		'''
		'''
		name = multiprocessing.current_process().name
		with s:
			pool.makeActive(name)
			sys.stderr.write("#################################")
			print "Now running: %s" %str(pool)
			runmakedb = "makeblastdb -in " + file + " -dbtype prot -title " + file + "_DB"
			print runmakedb
			os.system(runmakedb)
			pool.makeInactive(name)
		return None


	def run_blastp(self, file, db, sm, pool):
		'''
		Run blast mith multiprocessing for a search
		on a species database
		'''
		name = multiprocessing.current_process().name
		print file, db
		specie1 = file.split('/')[-1].split('.')[0]
		specie2 = db.split('/')[-1].split('.')[0]
		print specie1, specie2
		with sm:
			pool.makeActive(name)
			sys.stderr.write("#################################")
			print "Now running: %s" %str(pool)
			rundoblast = "blastp -query " + file + " -db " + db + " -outfmt \"7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen\" -seg yes -out "+self.blastp_repo+'/'+specie1+"-"+specie2
			print rundoblast
			os.system(rundoblast)
			pool.makeInactive(name)
		return None


	def doblast(self, species, cpu=1):
		'''
		Run blast on all the fasta files extracted from the 
		proteins present in the database
		Need to install blast: sudo apt-get install ncbi-blast+
		'''
		print "Run the blasts (%s CPUs)" %str(cpu)
		#curentpath=os.getcwd()
		#os.chdir(self.fasta_repo)
		listfa = glob.glob(self.fasta_repo + '/*.fasta')
		print listfa

		jobs = []
		pool = ActivePool()
		sm = multiprocessing.Semaphore(cpu)

		for file in listfa:
			if not os.path.isfile(file+".phr"):
				p = multiprocessing.Process(target=self.run_makeBlastDb, args=(file, sm, pool))
				jobs.append(p)

		for j in jobs:
			j.start()

		for j in jobs:
			j.join()
			sys.stderr.write('Now running: %s\n' % str(pool))


		jobs = []
		pool = ActivePool()
		sm   = multiprocessing.Semaphore(cpu)

		for specie1 in species:
			i = self.fasta_repo + '/' + specie1.replace(' ', '_') + '.fasta'
			for specie2 in species:
				j = self.fasta_repo + '/' + specie2.replace(' ', '_')	+ '.fasta'
#				if specie1 != specie2:
				check_rel = self.gd.check_species_relationship(specie1, specie2, 'Inparanoid_orthology')
				if check_rel == False:
					p = multiprocessing.Process(target = self.run_blastp, args = (j, i, sm, pool))
					jobs.append(p)
				else:
					print "Ortholog relations between %s and %s are already in the database" %(specie1, specie2)

		for j in jobs:
			j.start()

		for j in jobs:
			j.join()
			sys.stderr.write('Now running: %s\n' % str(pool))
		return None


	def doinparanoid(self, cpu=1):
		'''
		Run inparanoid on the blast results
		'''
		print "Run inparanoid"
		runinp="perl Inparanoid/SCRIPTS_INPARANOID.pl "+self.blastp_repo + " " + self.fasta_repo + " " + self.inparanoid_repo + " " + str(cpu)
		print runinp
		os.system(runinp)
		return None


	def parse_inparanoid(self, fname):
		'''
		Parse the inparanoid results, stored in sql tables
		'''
		orthologs = {}
		previous_group = 'NA'
		with open(fname) as file:
			file = file.read()
			for content in file.split('\n')[:-1]:
				tmp = content.split('\t')
				if not previous_group == tmp[0]:
					orthologs[tmp[0]] = {}
					previous_group = tmp[0]
				orthologs[tmp[0]][tmp[4]] = tmp[2]
		return orthologs


	def add_blast_orthology_relationship(self, orthologs):
		'''
		Add the orthology relationships present in the orthology dictionary
		into the neo4j graph database
		'''
		for og in orthologs.keys():
			for seq in orthologs[og]:
				seq1 = seq.rstrip()
				tabseq1 = seq1.split('|')
				specie1 = orthologs[og][seq].replace('.fasta', '').replace('_', ' ')
				for seq_2 in orthologs[og]:
					seq2 = seq_2.rstrip()
					tabseq2 = seq2.split('|')
					specie2 = orthologs[og][seq_2].replace('.fasta', '').replace('_', ' ')
					if seq1 != seq2:
						is_protein1 = self.gd.check_protein(tabseq1[1])
						is_protein2 = self.gd.check_protein(tabseq2[1])
						if (is_protein1 == True) and (is_protein2 == True):
							is_relationship = self.gd.check_relationship_link(tabseq1[1], tabseq2[1], 'Orthology_relationship')
							if is_relationship == False:
								self.gd.create_blastp_orthology_relationship(tabseq1[1], tabseq2[1])
		return None



if __name__ == '__main__':

	
	inparanoid_repo = '/home/salva/citrus2/BrasionSynth/inportho/'

	species = ['Citrus Clementina', 'Citrus sinensis', 'Arabidopsis thaliana', 'Physcomitrella patens']
	uniprot = uniprot_queries('Arabidopsis thaliana', '3702')
	chebi = chebi_from_string()
	chebi.chebi_connect()

	gd = graph_database(chebi, uniprot, '', '', 'neo4j', 'salva')
	gd.connect()
	inparanoid = inparanoid(gd, '/home/salva/citrus2/data/', '/home/salva/citrus2/BrasionSynth/blastp/', '/home/salva/citrus2/BrasionSynth/inportho/')

	#for specie in species:
	#	inparanoid.extract_fasta(specie)
	inparanoid.doblast(species)
	inparanoid.doinparanoid()

	resfiles = glob.glob(inparanoid_repo + '/sqltable/sqltable*')
	for fname in resfiles:
		orthologs = inparanoid.parse_inparanoid(fname)
		inparanoid.add_blast_orthology_relationship(orthologs)
	
	inparanoid.research_orthologs('Q8LEB2', 'Arabidopsis thaliana')