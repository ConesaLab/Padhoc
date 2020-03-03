
import os
from Bio import Entrez, Medline
from urllib2 import Request, urlopen, URLError
import json
import sys, httplib

class extract_articles():
    '''
    Class to extract the article (abstracts or complete text) from a query 
    Arguments:
    request, filelist, output repository, email, max number of articles, name output repository
    '''
    def __init__(self,req,filelist,rep,email,maxa,nametextrep):
        self.request = req
        self.email = email
        self.dirtxt = rep + '/' + nametextrep
        self.dirabrev = rep + '/' + 'abrev'
        os.mkdir(self.dirtxt)
        os.mkdir(self.dirabrev)
        self.filelist = filelist
        self.idlist = []
        self.maxa=maxa

    def pubmed_request(self):
        '''
        Execute a pubmed request and write the output pubmed id in the file outfile, 
        return the result of the request in xml format
        '''
        Entrez.email = self.email
        print self.request
        handle=Entrez.esearch(db="pubmed",
                            term=self.request,
                            retmax=self.maxa,
                            retmode='xml',
                            sort='relevance')
        record=Entrez.read(handle)
        self.idlist = record['IdList']
        if len(self.idlist) > 0:
            self.fetch_details_xml() # extract abstracts
            self.extract_text()      # extract full text
            self.extract_abrev()
        else:
            print 'No text for the path searched'
        return self.idlist
    
    def read_pmidlist(self,filepmid):
        '''
        Read a list of pmid of interrest and extract uniquely this publications
        '''
        with open(filepmid,"r") as ins:
            for line in ins:
                self.idlist.append(line.replace('\n', ''))
        self.fetch_details_xml() # extract abstracts
        self.extract_text()      # extract full text
        self.extract_abrev()     # extract the abreviations from the texts
        return self.idlist
    
    def fetch_details_xml(self):
        '''
        Retrieve title and abstract corresponding to a pubmedID
        INPUTS:
        id_list: list of the pubmed id
        repsave: directory where title and abstracts should be saved
        RETURN: nothing
        '''
        for n in range(10):
            try:
                ids = ','.join(self.idlist)
                Entrez.email = self.email
                handle = Entrez.efetch(db='pubmed',
                                    retmode='xml',
                                    id=ids)
                results = Entrez.read(handle)
            except httplib.IncompleteRead:
                pass
        for i in results['PubmedArticle']:
            f=open("%s/%s.txt" % (self.dirtxt, i['MedlineCitation']['PMID']),'w')
            towrite=i['MedlineCitation']['Article']['ArticleTitle']
            if towrite: # test if the string exist
                #print "title "+towrite
                if not towrite.isspace(): # test if the string is empty (only white characters)
                    #print "after isspace"
                    f.write("%s\n" % towrite.encode('ascii','ignore') )
                    #print towrite
            
            try:
                for j in i['MedlineCitation']['Article']['Abstract']['AbstractText']:
                    towrite=j
                    #print "abstract "+towrite
                    if not towrite.isspace():
                        #print "after isspace"
                        f.write("%s\n" % j.encode('ascii','ignore'))
            except KeyError:
                print("Warning: no abstract for %s" % i['MedlineCitation']['PMID'])
            f.close()
        return None

    def extract_text(self):
        '''
        Extract complete article from PMC (not available for all the PMID)
        rehand is a record handle from the pubmed request
        self.idlist is the list of the pmid result from the pubmed query
        return nothing
        '''
        for i in self.idlist:
            request = Request("https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/%s/ascii" % (i))
            try:
                response = urlopen(request)
                jsontext = response.read()
                jsonout = json.loads(jsontext)
                f=open("%s/%s.txt" % (self.dirtxt, i),'w')
                for j in jsonout['documents']: #for each element of the list inside documents (normally only one element)
                    for p in j['passages']: # for each element of the list passages 
                        if p['infons']['type'] != 'table':
                            if p['infons']['type'] != 'ref':
                                if p['infons']['section_type'] != 'METHODS':
                                    if not p['text'].isspace(): #test if the string contain only whitespaces
                                        f.write('%s\n' % p['text'])
                f.close()
                #print ("complete text for %s" % i)
            except URLError, e:
                #print 'warning : No complete text for %s ' % i
                pass
        return None

    def extract_abrev(self):
        '''
        Extract abreviation foreach text (abstract of complete text) and save them in a file
        In the repository self.dirabrev the abrev foreach files correspond to the rows of the file pmid.tmp 
        '''
        path=os.getcwd()
        repoPath = os.path.dirname(os.path.realpath("extract_articles.py"))
        os.chdir(repoPath+'/tools')
        for i in self.idlist:
            #toexecute='java ExtractAbbrev ../'+self.dirtxt+'/'+i+'.txt > ../'+self.dirabrev+'/'+i+'.tmp'
            #print(toexecute)
            if os.path.isdir('../'+self.dirtxt):# relative path to the repositories
                toexecute='java ExtractAbbrev ../'+self.dirtxt+'/'+i+'.txt > ../'+self.dirabrev+'/'+i+'.tmp'
            else: #complete path
                toexecute='java ExtractAbbrev '+self.dirtxt+'/'+i+'.txt > '+self.dirabrev+'/'+i+'.tmp'
            os.system(toexecute)
            #subprocess.check_call()
        os.chdir(path)
        return None


if __name__ == '__main__':

    # articles = extract_articles('Histone acetylation AND Homo Sapiens[MeSH Terms]', "", '/home/salva/pubmed_trial/', 'pereiracecile@live.fr', 1000, '')
    #listpmid = articles.pubmed_request()
    articles = extract_articles("", "", '/home/salva/panto_mm/', 'pereiracecile@live.fr', "", "")
    articles.read_pmidlist('/home/salva/pmid_trial.txt')