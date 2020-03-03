#!/usr/bin/env python

import os, glob
import sys, re

class run_vagrant_metrecon():

	def __init__(self, in_rep, out_rep, vagrant_rep):

		self.in_rep = in_rep			# Location of the text files
		self.out_rep = out_rep			# Location of the future output
		self.vagrant_rep = vagrant_rep

	def split_files(self, randomString, n=50):
		'''
		Split the articles in groups of n files and
		run metrecon on them
		'''
		current_path = os.getcwd()
		os.chdir(self.in_rep)
		files = glob.glob('*.txt')		# Text files
		t = 0; part = 1
		
		text_subset = randomString + '_P' + str(part)
		temp_subset = '/tmp/' + text_subset
		
		if not os.path.exists(self.out_rep):
			os.mkdir(self.out_rep)
		
		while os.path.isdir(temp_subset):
			part=part+1
			text_subset = randomString + '_P' + str(part)
			temp_subset = '/tmp/' + text_subset

		# TODO: Ask Cecile if the while loop should be included. 
		# Isn't it covered due to the randomString?: CECILE: the aim of the while loop is to find the first new part available, it allow to run this function twice with two different cut values (in order to deal with metrecon errors)
		# may be not necessary if after the run the tmp file is destroy

		os.mkdir(temp_subset)			# Create temporary folder to store subset of files

		# TODO: Is the dictionary with the a1a2 files needed? What for?
		# CECILE: When you run this function for the second time, this dictionary store the files already done (a1/a2 available)
		# It allow to redo only the files with errors.

		ptmp2=self.out_rep+"/a1a2_*/*"
		filesdone=glob.glob(ptmp2)
		dicofilesdone={} #dictionary with all the a1 and a2 files already done
		for i in filesdone:
		    tmp=i.split('/')
		    tmp[-1]=re.sub('a2$','txt',tmp[-1])
		    dicofilesdone[tmp[-1]]=""
		    print tmp[-1]
		
		for f in files:
			if not f in dicofilesdone.keys(): #file not already analysed by metrecon
				t += 1
				toRun = "cp %(file)s %(temp_rep)s" % {'file': f, 'temp_rep': temp_subset}
				os.system(toRun)			# Copy text files to the temp subfolder
				if t%n == 0:
					coprepvag = 'cp -r %(text_path)s %(vagrant_path)s' \
					% {'text_path':temp_subset, 'vagrant_path': self.vagrant_rep}
					os.system(coprepvag)	# Copy text subset to vagrant machine
					vagrant_out_subset = 'a1a2_' + randomString + '_P' + str(part)
					os.mkdir(self.vagrant_rep + vagrant_out_subset)
					self.run(text_subset, vagrant_out_subset)
					part = part + 1
					text_subset = randomString + '_P' + str(part)
					temp_subset = '/tmp/' + text_subset
					os.mkdir(temp_subset)

			#In case there are more files than a multiple of n
			# CECILE: I think that this will not work with the possibility to run twice the function (to rerun metrecon errors)
				#elif (part == ((len(files)//n) + 1)) and (t%n > 0):
				#	if t == len(files):
				#		coprepvag = 'cp -r %(text_path)s %(vagrant_path)s' \
				#		% {'text_path':temp_subset, 'vagrant_path': self.vagrant_rep}
				#		os.system(coprepvag)	# Copy text subset to vagrant machine
				#		vagrant_out_subset = 'a1a2_' + randomString + '_P' + str(part)
				#		self.run(text_subset, vagrant_out_subset)
		
		#test if the last part repository is empty (less or exactly cut files)
		nbflastpart=len([name for name in os.listdir(temp_subset) if os.path.isfile(temp_subset+'/'+name)])
		if nbflastpart>0:
			#for the last part
			coprepvag = 'cp -r %(text_path)s %(vagrant_path)s' \
			% {'text_path':temp_subset, 'vagrant_path': self.vagrant_rep}
			os.system(coprepvag)	# Copy text subset to vagrant machine
			vagrant_out_subset = 'a1a2_' + randomString + '_P' + str(part)
			self.run(text_subset, vagrant_out_subset)
		os.chdir(current_path)
		
		#removeFiles = 'rm -r ' + self.vagrant_rep + '/text' + randomString + '_P*' #correction cecile
		#os.system(removeFiles)
		removeFiles = 'rm -r ' + self.vagrant_rep + '/' + text_subset
		os.system(removeFiles)
		return None

	def run(self, in_rep, out_rep):
		'''
		Run metrecon
		'''
		current_path = os.getcwd()
		os.chdir(self.vagrant_rep)
		os.system("vagrant up")
		torun="vagrant ssh -c \"sh bin/classification -i /vagrant/"+in_rep+" -o /vagrant/"+out_rep+"\""
		os.system(torun); os.chdir(current_path)
		toRun = 'mv ' + self.vagrant_rep + '/' + out_rep + ' ' + self.out_rep
		os.system(toRun)
		return None

if __name__ == '__main__':

	in_rep = '/home/salva/ecoli_comparison/textCYDNKJ5D/'
	out_rep = '/home/salva/extractVagrant/'
	metrecon = '/home/salva/metrecon/'
	randomString = 'fsdafsdasfda'

	metrecon = run_vagrant_metrecon(in_rep, out_rep, metrecon)
	metrecon.split_files(randomString)
