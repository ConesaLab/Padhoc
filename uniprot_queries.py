#!/usr/bin/env python

import re, time
from bioservices import UniProt
import zeep, sys
from zeep import Client
from modules import calculateSimilarity
import httplib
from sys import argv


def clean(word,specie=''):
    '''
    Clean words to allow them comparison.
    See if I need to add all the cases that I covered with perl before or if the score in the chebi comparison is enouth
    '''
    #lowercases
    cleanw=word.lower()
    specie=specie.lower()
    sp=specie.split(' ')
    specie=specie.strip()
    #no whitecharacters
    cleanw=cleanw.strip()
    #remove parenthesis
    cleanw=re.sub('^(r)-','',cleanw)
    cleanw=re.sub('^(s)-','',cleanw)
    cleanw=re.sub('\)n$','',cleanw)
    cleanw=re.sub('\)m$','',cleanw)
    cleanw=re.sub('^\(\+?\-?\d{0,5}\)','',cleanw)
    cleanw=re.sub('\(','',cleanw)
    cleanw=re.sub('\)','',cleanw)
    cleanw=re.sub('\[','',cleanw)
    cleanw=re.sub('\]','',cleanw)
    cleanw=re.sub('\{','',cleanw)
    cleanw=re.sub('\}','',cleanw)
    cleanw=re.sub('^l-','',cleanw)
    cleanw=re.sub('^d-','',cleanw)
    cleanw=re.sub('^r-','',cleanw)
    cleanw=re.sub('^s-','',cleanw)
    cleanw=re.sub('^ec:','',cleanw)
    cleanw=re.sub('^ec ','',cleanw)
    cleanw=re.sub('^alpha-','',cleanw)
    cleanw=re.sub('^beta-','',cleanw)
    cleanw=re.sub('^[0-9]+','',cleanw)
    cleanw=re.sub(' genes{0,1}','',cleanw)
    cleanw=re.sub(' proteins{0,1}','',cleanw)
    #cleanw=re.sub('\)','',cleanw)
    #remove words never used in the text
    cleanw=re.sub(' atom','',cleanw)
    cleanw=re.sub('s$','',cleanw) #remove plural
    cleanw=cleanw.replace(' ','') #test 3 july
    cleanw=cleanw.replace('"','')
    cleanw=cleanw.replace('-','')
    cleanw=cleanw.replace("'",'')
    cleanw=cleanw.replace("icacid",'ate')
    cleanw = re.sub("\d+$", "", cleanw)
    for s in sp:
        cleanw=re.sub(s,'',cleanw)#remove the specie name from the entity name (ex for Arabidopsis thaliana: remove arabidopsis and then remove thaliana)
    if cleanw=='':
        cleanw='empty_string'
    return cleanw

class uniprot_queries():
    '''
    Interrogation of the uniprot database
    '''

    def __init__(self,specie,taxid):

        self.query = ""
        self.specie = specie
        self.taxid = taxid
        self.dicotested = {} #avoid to test several time the same ID
        self.unip = UniProt()
              

    def create_query(self, query):
        '''
        Create the query that will be used to search the database
        '''
        self.query = query
        self.query = re.sub('\+','\\+',self.query)
        self.query = re.sub('"','',self.query)
        self.query = self.query.replace(':','')
        self.query = self.query.replace(';','')

        cq=''
        #change cq in function of the specie specifications
        if self.specie != '':
            if self.taxid != '':
                cq = self.query + ' AND organism: ' + self.specie + ' [' + self.taxid + ']'
            else:
                cq = self.query + ' AND organism: ' + self.specie
        else:
            cq = self.query
        return cq


    def uniprot_request(self, cq):
        '''
        Request to the uniprot database using the uniprot
        module from bioservices
        '''
        for n in range(10):
            try:
                cq=re.sub('^\[','',cq) # remove the [ when it start by it 
                d = self.unip.quick_search(cq,limit=1) #I take the best result of the match
                break
            except(zeep.exceptions.Fault, zeep.exceptions.TransportError, TypeError, AttributeError):
                d={}
                pass
        return d

    def query_id(self, query):
        '''
        Use bioservice
        research uniprot ID corresponding to a query
        RETURN: uniprotID, uniprot entry name, uniprot gene names
        '''
        cq = self.create_query(query)
        if not cq in self.dicotested.keys():#cq not already tested
            d = self.uniprot_request(cq)
            if(len(d.keys())>0): #at least one result
                unipID = d.keys()[0]
                for gene in d[unipID]['Gene names'].split(' '):
                    if clean(query) == clean(gene):
                        tabtosave=[unipID, d[unipID]['Entry name'], d[unipID]['Gene names'], d[unipID]['Protein names']]
                        self.dicotested[cq] = tabtosave
                        return unipID, d[unipID]['Entry name'], d[unipID]['Gene names'], d[unipID]['Protein names']
                    if len(query.split('_')) > 1:
                        query = query.split('_')[0]
                    if clean(query) == clean(gene):
                        tabtosave=[unipID, d[unipID]['Entry name'], d[unipID]['Gene names'], d[unipID]['Protein names']]
                        self.dicotested[cq] = tabtosave
                        return unipID, d[unipID]['Entry name'], d[unipID]['Gene names'], d[unipID]['Protein names']
                for protein in d[unipID]['Protein names'].split('('):
                    protein = protein.strip()
                    protein = re.sub("\)$",'',protein)
                    protein = protein.replace('+','')
                    if clean(query) == clean(protein):
                        tabtosave = [unipID, d[unipID]['Entry name'], d[unipID]['Gene names'], d[unipID]['Protein names']]
                        self.dicotested[cq] = tabtosave
                        return unipID, d[unipID]['Entry name'], d[unipID]['Gene names'], d[unipID]['Protein names']
                    #By Cecile:
                    if(re.search(clean(protein),clean(query))):
                        tabtosave=[unipID, d[unipID]['Entry name'], d[unipID]['Gene names'], d[unipID]['Protein names']]
                        self.dicotested[cq] = tabtosave
                        return unipID, d[unipID]['Entry name'], d[unipID]['Gene names'], d[unipID]['Protein names']
            tabtosave = ['','','','']
            self.dicotested[cq] = tabtosave
            return '', '', '', ''
        else:
            return self.dicotested[cq][0], self.dicotested[cq][1], self.dicotested[cq][2], self.dicotested[cq][3]


    def mapping_id(self, prot_id, database):
        '''
        Map a protein id to obtain the information from the
        uniprot database

        This function became deprecated after bioservices uniprot stopped working, 
        not allowing the retrieval of the mapping id.
        In turn, I created the mapping_id2 functiion to make API requests to uniprot db
        '''
        if not prot_id in self.dicotested.keys():
            #mapping = self.unip.mapping(fr=database, to='ID', query = prot_id) # No longer used since bioservices stoped working
            #print mapping
            #sys.exit()
            mapping = self.mapping_id2(prot_id, database)
            if len(mapping) > 1:
                for k in mapping:
                    value = k.split('\t')[1]
                    break
            else:
                tabtosave = ['','','','']
                self.dicotested[prot_id] = tabtosave
                return '', '', '', ''
            cq = self.create_query(value)
            r = self.unip.quick_search(value, limit=1)
            if len(r.keys()) > 0:
                unipID = r.keys()[0]
                tabtosave = [unipID, r[unipID]['Entry name'], r[unipID]['Gene names'], r[unipID]['Protein names']]
                self.dicotested[prot_id] = tabtosave
                return unipID, r[unipID]['Entry name'], r[unipID]['Gene names'], r[unipID]['Protein names']
            tabtosave = ['','','','']
            self.dicotested[cq] = tabtosave
            return '', '', '', ''
        else:
            return self.dicotested[prot_id][0], self.dicotested[prot_id][1], self.dicotested[prot_id][2], self.dicotested[prot_id][3]


    def mapping_id2(self, prot_id, database):
        '''
        '''
        import urllib,urllib2

        url = 'https://www.uniprot.org/uploadlists/'

        params = {
        'from': database,
        'to':'ID',
        'format':'tab',
        'query': prot_id
        }

        data = urllib.urlencode(params)
        request = urllib2.Request(url, data)
        contact = "" # Please set a contact email address here to help us debug in case of problems (see https://www.uniprot.org/help/privacy).
        request.add_header('User-Agent', 'Python %s' % contact)
        for n in range(10): # deal with problems to connect to the server
            try: 
                response = urllib2.urlopen(request)
                page = response.read(200000)
                page = page.split('\n')[1:]
                break
            except (urllib2.HTTPError, httplib.BadStatusLine):
                page = ''
                pass
        return page


    def uni_search(self, protein):
        '''
        Do the quick search but not trying to find the 
        '''
        if not protein in self.dicotested.keys():
            cq = self.create_query(protein)
            r = self.unip.quick_search(cq, limit=1)
            unipID = r.keys()[0]
            tabtosave = [unipID, r[unipID]['Entry name'], r[unipID]['Gene names'], r[unipID]['Protein names']]
            self.dicotested[protein] = tabtosave
            return unipID, r[unipID]['Entry name'], r[unipID]['Gene names'], r[unipID]['Protein names']
        else:
            return self.dicotested[protein][0], self.dicotested[protein][1], self.dicotested[protein][2], self.dicotested[protein][3]
        return None


if __name__ == '__main__':

    lista = []
    #uniprot = uniprot_queries('Homo sapiens', '9606')
    #print uniprot.query_id('rpoB')
    uniprot = uniprot_queries('', '')

    a = open(argv[1]).read()
    for gene in a.split('\n'):
        uniID = uniprot.query_id(gene)[0]
        if uniID != '':
            lista.append(uniID.encode('utf-8'))
    print lista
    sys.exit()

    uniprot = uniprot_queries('Homo sapiens', '9606')
    print uniprot.query_id('PanD complex')

    sys.exit()


    'EC-4.1.1.11'
    uniprot = uniprot_queries('Physcomitrella patens', '3218')
    print uniprot.query_id('(+)-endo-beta-bergamotene synthase ((2Z,6Z)-farnesyl diphosphate')

    sys.exit()
    uniprot = uniprot_queries('Homo sapiens', '9606')
    a,b,c,d = uniprot.query_id('3-PG')
    print a, b, c, d
    sys.exit()
    a,b,c,d = uniprot.query_id('alcohol dehydrogenase')
    a,b,c,d = uniprot.query_id('MYCD_HUMAN') 
    a, b, c, d = uniprot.mapping_id('MYCD_HUMAN', 'ID')
    sys.exit()

    a, b, c, d = uniprot.mapping_id('SRF_HUMAN', 'ID')