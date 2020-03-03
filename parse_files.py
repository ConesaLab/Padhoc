#!/usr/bin/env python

import xml.etree.ElementTree as ET
import sys

class parse_file():
    
    def __init__(self,nerfile,artTsentence,dicoPMID_idDoc,r):
        '''
        Initialize the parsing of the file placeholder-preprocessed.xml.gz-ner.xml.gz (WARNING, it's not a gz file...)
        Retreave the sentences of the T (compounds or proteins) 
        '''
        self.file=nerfile
        #self.artTsentence={}
        #self.dicoPMID_idDoc={}
        self.artTsentence=artTsentence
        self.dicoPMID_idDoc=dicoPMID_idDoc
        self.r=r
        
    def docTsent(self):
        tree = ET.parse(self.file)
        root = tree.getroot()
        n=1 #number of the element in the a1 file
        for doc in root: #documents (complete article/abstract)
            #print(doc.attrib['origId'])
            iddoctmp=self.r+'.'+doc.attrib['id'] 
            print iddoctmp
            print doc.attrib['origId']
            #self.dicoPMID_idDoc[doc.attrib['origId']]=[doc.attrib['id']]
            self.dicoPMID_idDoc[doc.attrib['origId']]=[iddoctmp]
            self.artTsentence[doc.attrib['origId']]={}
            for sentence in doc: #sentences inside each document
                #print sentence.attrib['id']
                for node in sentence: #compound or protein
                    #print node.attrib['id']
                    T='T'+str(n) # ID of the element in the a1 file
                    self.artTsentence[doc.attrib['origId']][T]=sentence.attrib['text'].replace('"',"'") # dictionary with keys [PMID][T0]=sentence
                    n=n+1
            n=1
            
        return self.artTsentence, self.dicoPMID_idDoc
