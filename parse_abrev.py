#!/usr/bin/env python

import glob, os

class parse_abrev():
    def __init__(self,pathabrev):
        '''
        initializethe parse abrev class
        '''
        self.pathabrev=pathabrev
        self.dicabrev={}
        self.listfiles=[]

    def parse(self):
        '''
        Parse the abrev files and put the result in a dictionary
        Return the dictionary
        '''
        for file in glob.glob(self.pathabrev + '/*.tmp'):
            #print(file)
            self.listfiles.append(file)
        for fil in self.listfiles:
            base=os.path.basename(fil)
            filename=os.path.splitext(base)[0]
            #print(filename)
            #print('filename ' + filename)
            self.dicabrev[filename]={}
            with open(fil) as f:
                for line in f:
                    line = line.rstrip()
                    line = line.split('\t')
                    if(len(line)==2):
                        self.dicabrev[filename][line[0]]=line[1]
                        #print(filename+' '+line[0]+' '+line[1])
        return self.dicabrev
