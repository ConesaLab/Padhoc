#!/usr/bin/env python

import os, glob, sys
from sys import argv


class run_tees():

	def __init__(self, out_repository):
            '''
            Initiate class run_tees
            input: tees repository
            '''
            self.repository = out_repository

	
	def split_input(self, randomString, tees_rep, corpora, in_rep, n=50):
		'''
		Split input files by n, so the tees procedure isn't overloaded
		randomString: the value of the randomString caracterizing this run of the pipeline
		tees_rep: the repository with the classify file 
		corpora: the tested corpora
		in_rep: repository with all the article files (txt)
		n: number of files in one split.
		'''
		corpora_name = corpora
		if len(corpora.split('/')) > 1:
			corpora_name = corpora.split('/')[-2]
		corpora_repository = self.repository + '/' + corpora_name
		if not os.path.exists(self.repository): # test if the general output repository for TEES exist
			os.mkdir(self.repository)
		if not os.path.exists(corpora_repository):
			print corpora_repository
			os.mkdir(corpora_repository)
		direct = os.getcwd()
		files=glob.glob(in_rep + '/*.txt')
		p = 0 ; t = 0
		tmp = '/tmp/txt' + randomString + '_P' + str(p) + '_' + corpora_name
		out_rep = '/' + 'a1a2_' + randomString + '_P' + str(p)

		print tmp
		while os.path.isdir(tmp):
			p=p+1
			tmp = '/tmp/txt' + randomString + '_P' + str(p) + '_' + corpora_name
			out_rep = '/' + 'a1a2_' + randomString + '_P' + str(p)
		os.mkdir(tmp)
		
		for f in files:
			t += 1
			toRun = ("cp %(file)s %(tempRep)s" % {'file':f, 'tempRep':tmp})
			print toRun
			os.system(toRun)
			if t%n==0:
				self.run(tees_rep, corpora, tmp, out_rep, corpora_name)
				print t, len(files)
				print 'asfdasfdasfdasfdasfdasfdasfdasfdasfda'
				p += 1
				tmp = '/tmp/txt' + randomString + '_P' + str(p) + '_' + corpora_name
				print tmp
				out_rep = '/' + 'a1a2_' + randomString + '_P' + str(p)
				os.mkdir(tmp)

		#To also run the last set of articles that do not sum up to the threshold n
		print os.listdir(tmp)
		nbflastpart=len([name for name in os.listdir(tmp) if os.path.isfile(tmp+'/'+name)])
		print tmp
		print nbflastpart
		if nbflastpart>0:
			self.run(tees_rep,corpora,tmp,out_rep, corpora_name)
		return None


	def run(self, tees_rep, corpora, in_rep, out_rep, corpora_name):
		'''
		Run the TEES program using the corpora of interest
		'''
		direct = os.getcwd()
		os.chdir(tees_rep)
		output = self.repository + corpora_name + '/' + out_rep
		print '\n\n\n\n\n\n\n\n'
		print output
		os.mkdir(output)
		prefix = output[-2:]
		print '\n\n\n\n\n',output, prefix
		toRun = "python2 classify.py -m %(model)s -i %(input)s -o %(out_rep)s"%{'model':corpora, 'input':in_rep, 'out_rep':output + '/'}
		#print toRun
		os.system(toRun)
		os.chdir(output)
		toRun = 'tar -xvzf -events.tar.gz'
		os.system(toRun)
		os.chdir(direct)
		# CECILE: I don't think that we run tmChem for the other TEES corpora (they are for proteins no?)
		toRun = 'cp %(tees_rep)s/sentences.pubtator.tmChem %(out_rep)s' % {'tees_rep':tees_rep, 'out_rep': output}
		os.system(toRun)
		return None


if __name__ == '__main__':

	in_rep = argv[1]; out_rep = argv[2]
	tees_rep = '/home/salva/TEES-master/'
	randomString = 'hkhaskjsajkjsa'
	tees = run_tees(out_rep)
	tees.split_input(randomString, tees_rep, 'GE11', in_rep)