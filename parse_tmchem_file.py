#!/usr/bin/env python

import os, re
from libchebipy import ChebiEntity

class parse_tmchemfile(): 
    def __init__(self,ptmchem,dictmchem,r):
        '''
        class to parse the tmchemfile that contain normalize names for the compounds
        dictmchem[a1a2_textHT0RNL37_P1.TEES.d5][compound]=chebiID
        '''
        self.file=ptmchem
        self.dictmchem=dictmchem
        self.tmpr=r
    
    def parse(self):
        '''
        parse tmchem file (sentence.pubtator.tmchem), save in a dictionary the sentence id, the name of the compound and the chebi name
        Save only the chebi ID, not the mesh
        '''
        with open(self.file) as f:
            for line in f:
                line=line.rstrip()
                line=line.split('\t')
                if(len(line)==6):
                    #print(line)
                    #print(line[0])
                    iddoc=line[0] #id of the sentence in TEES: ex TEES.d0.s1
                    #print(iddoc)
                    iddoc=os.path.splitext(iddoc)[0] #keep only the ID of the doc: ex TEES.d0
                    iddoc=self.tmpr+'.'+iddoc
                    #print(iddoc)
                    
                    # Some chebi are parent and some are children (duplicate and in that case, we recover and take the parent)
                    chebimesh=line[5]
                    chebi=re.search('(?<=CHEBI:)\w+',chebimesh)
                    if(chebi):
                        chebientity=ChebiEntity("CHEBI:"+chebi.group(0))
                        parent=chebientity.get_parent_id() # ex: CHEBI:22982 as for parent CHEBI:27732
                        if parent:
                            #print('take the parent: '+chebimesh+' '+parent)
                            chebimesh=parent
                    
                        if iddoc in self.dictmchem.keys():
                            #print line[3]
                            self.dictmchem[iddoc][line[3]]=chebimesh #iddoc: document ID (TEES), name as find in the text, chebimesh=chebi ID or mesh ID
                        else:
                            self.dictmchem[iddoc]={}
                            self.dictmchem[iddoc][line[3]]=chebimesh
        return self.dictmchem
