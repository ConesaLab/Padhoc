
import zeep, time
from zeep import Client
import re, requests
from libchebipy import ChebiEntity
from random import randint
from modules import calculateSimilarity

import sys

class chebi_from_string(): #TODO: test to see what we get in function of the input string. If input string== chebi name then keep chebi. Clean function to compare both?
    
    '''
    This class should be initialized just once, the compound name to test
    should be added in the chebi() function

    TO DO:
    1. Add a dictionary of short, common, reliable compounds (e.g. ATP, H2O)
    2. Check the clean function on the metabolites
    3. Delete either the clean or the normal at the analyzelistchebi() function
    '''

    def __init__(self,name=None):
        '''
        search in the chebi database a chebi number corresponding to the chebi name
        '''
        self.chebi_dictio = {}
        self.long=name
        self.chebiid=""
        self.nameid={} #name to id (include synonyms)
        self.idname={}
        self.idchebi={} 

    
    def clean(self, word,specie=''):
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
        cleanw=re.sub('<[^>]+>','',cleanw)
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
        for s in sp:
            cleanw=re.sub(s,'',cleanw)#remove the specie name from the entity name (ex for Arabidopsis thaliana: remove arabidopsis and then remove thaliana)
        if cleanw=='':
            cleanw='empty_string'
        return cleanw


    def chebi_connect(self):
        '''
        Connect to the chebi client
        '''
        wsdl='http://www.ebi.ac.uk/webservices/chebi/2.0/webservice?wsdl'
        transport = zeep.Transport(operation_timeout=300)
        self.client = zeep.Client(wsdl=wsdl, transport=transport)
        return None



    def chebi_from_dictio(self, dicofile):
        '''
        Cecile Pereira: 14 August #TODO TEST
        Read the dictionary file,
        Save the info of name, id, chebi in dictionaries.
        '''
        #print "Read chebi dictionary file"
        f=open(dicofile,'r')
        tmpid=""
        for line in f:
            if line.startswith('ID '):
                tmp=split(' ',line)
                tmpid=tmp[1]
                self.idname[tmpid]={}
            elif line.startswith('NA '):
                tmp=split(' ',line)
                tmpname=clean(tmp[1])
                self.nameid[tmpname]=tmpid
                self.idname[tmpid][tmpname]=''
            elif line.startswith('TM '):
                tmp=split(' ',line) #synonyms
                tmp=split('\@',tmp[1])
                tmp[0]=clean(tmp[0])
                self.nameid[tmp[0]]=tmpid
                self.idname[tmpid][tmp[0]]='' #all the names of the ID
            elif line.startswith('DB CHEB_'):
                tmp=split(' ',line)
                for t in tmp:
                    t=re.sub('^CHEB_','',t)
                    t=re.sub('^CHEBI:','',t)
                    t=re.sub('Not available','',t)
                    if not t == '':
                        chebi=t
                        if(chebi):
                            chebientity=ChebiEntity("CHEBI:"+chebi)
                            parent=chebientity.get_parent_id() # ex: CHEBI:22982 as for parent CHEBI:27732
                            if parent: #take the CHEBI parent...
                                self.idchebi[tmpid][parent]=''



    def chebi(self, compName, score=10):
        '''
        http://docs.python-zeep.org/en/master/
        http://www.ebi.ac.uk/chebi/webServices.do
        '''
        chebiId=None; chebiName=None
        if self.clean(compName) in self.chebi_dictio.keys(): #if the name was already found in chebi
            chebiId = self.chebi_dictio[self.clean(compName)][0] #modification cecile 6 July: add clean function 
            chebiName = self.chebi_dictio[self.clean(compName)][1]
        else:
            client=None
            for n in range(10): # deal with problems to connect to the server
                try:
                    #wsdl='http://www.ebi.ac.uk/webservices/chebi/2.0/webservice?wsdl'
                    #client = zeep.Client(wsdl=wsdl)
                    listchebiname=self.client.service.getLiteEntity(compName,"ALL NAMES",5,'ALL')
                    #print listchebiname
                    break
                except (zeep.exceptions.Fault, zeep.exceptions.TransportError, requests.exceptions.ConnectionError, requests.exceptions.ReadTimeout, TypeError):
                    listchebiname = None
                    pass
            if listchebiname!=None: #at least one result
                chebiId,chebiName=self.analyzelistchebi(listchebiname,score, compName)
                if chebiId!=None:
                    for n in range(10): # deal with problems to connect to the server
                        try:
                            chebisyn=self.client.service.getCompleteEntity(chebiId)
                            break
                        except (zeep.exceptions.Fault, zeep.exceptions.TransportError, requests.exceptions.ConnectionError, requests.exceptions.ReadTimeout):
                            chebisyn = None
                            pass
                    for s in chebisyn.Synonyms: # Cecile 25 aout: save the synonyms in the dictionary
                        self.chebi_dictio[self.clean(s.data)]=[chebiId,chebiName]
        return chebiId, chebiName #return also chebi ascii name
    

    def analyzelistchebi(self, listchebiname, score, compName):
        '''
        Analyze the result of the request to chebi
        First check if the compound is present among the results, 
        otherwise look for it in the synonyms, and if non of them worked, 
        obtian the first instance from the request.
        '''
        # Salva: Here the function is repeated twce, it may happen that the dictionary 
        # is being overloaded with unusefull information, either clean the compound before
        # or not, but we shoudln't do everything twice.
        for entity in listchebiname: # test if name identical (before cleanm it contains more information)
            chebiId = entity.chebiId
            chebiName = entity.chebiAsciiName
            if compName == chebiName:
                self.chebi_dictio[self.clean(compName)] = [chebiId, chebiName]
                #print "long eq chebiName 1 "+compName+" "+chebiName
                return chebiId, chebiName
            # Salva: the requests should include try except to avoid errors.
            
            for n in range(10):
                try:
                    chebisyn=self.client.service.getCompleteEntity(chebiId) #synonyms
                    break
                except zeep.exceptions.TransportError, zeep.exceptions.Fault:
                    chebisyn = []
                    pass

            for s in chebisyn.Synonyms:
                if s.data == compName:
                    #print "long eq chebiName 2"+compName+" "+s.data
                    self.chebi_dictio[self.clean(compName)] = [chebiId, chebiName]
                    return chebiId, chebiName
        for entity in listchebiname:#test if name after clean are identical
            chebiId = entity.chebiId
            chebiName = entity.chebiAsciiName
            #print "second for "+chebiName+' '+self.clean(compName)
            if self.clean(compName) == self.clean(chebiName):
                #print "clean long clean chebiname"+compName+" "+chebiName
                self.chebi_dictio[self.clean(compName)] = [chebiId, chebiName]
                return chebiId, chebiName
            for n in range(10):
                try:
                    chebisyn=self.client.service.getCompleteEntity(chebiId) #synonyms
                except (zeep.exceptions.TransportError, zeep.exceptions.Fault, requests.exceptions.ConnectionError, requests.exceptions.ReadTimeout):
                    chebisyn = []
                    pass
#            for s in chebisyn.Synonyms:
#                print s
#                if self.clean(s.data)== self.clean(compName):
                    #print "clean synonyms"+compName+" "+s.data
#                    self.chebi_dictio[self.clean(compName)] = [chebiId, chebiName]
#                    return chebiId, chebiName
        if listchebiname[0].searchScore >= score:
            chebiId = listchebiname[0].chebiId
            chebiName = listchebiname[0].chebiAsciiName
            self.chebi_dictio[self.clean(compName)] = [chebiId, chebiName]
            return chebiId, chebiName
        self.chebi_dictio[self.clean(compName)] = [None, None]
        return None,None


    def obtain_chebi_parents(self, chid, ontology_type = ''):
        '''
        Obtain the parents of a given chebi ID
        '''
        ontologies = []
        ontology = self.client.service.getOntologyParents(chid)
        if ontology_type != '':
            for ont in ontology:
                if ont['type'] == ontology_type:
                    ontologies.append(ont)
        else:
            ontologies = ontology
        return ontologies


    def obtain_chebi_children(self, chid, ontology_type = ''):
        '''
        Obtain the parents of a given chebi ID
        '''
        ontologies = []
        ontology = self.client.service.getOntologyChildren(chid)
        if ontology_type != '':
            for ont in ontology:
                if ont['type'] == ontology_type:
                    ontologies.append(ont)
        else:
            ontologies = ontology
        return ontologies


    def obtain_ontology_chid(self, chid):
        '''
        Obtain the desired ontologic order for the
        ChEBI ID
        '''
 #       print chid
        completeChid = self.client.service.getCompleteEntity(chid)
        if completeChid['charge'] != None:
            neutralChid, metabolites, charged = self.get_neutral_id(completeChid, [], [])
            if metabolites != False:
                neutralChid = self.obtain_lowest_charge(charged)
            generalChid = self.get_general_id(neutralChid)
            return generalChid
        else:
            specificChid = self.get_specific_id(completeChid)
            if specificChid == completeChid:
                return completeChid
            else:
                return self.obtain_ontology_chid(specificChid['chebiId'])


    def get_neutral_id(self, completeChid, metabolites = [], charged = []):
        '''
        Obtain the metabolite that is zero charged (or closest to 0)
        '''
        charge = completeChid['charge']; chebiID = completeChid['chebiId']
        charged.append((chebiID, charge)); metabolites.append(chebiID)
        # Create list with all the 'synonyms'
        if int(charge) == 0:
            metabolites = False
            return completeChid, metabolites, charged
        elif int(charge) < 0:
            chid = self.get_parent_from_complete(completeChid, 'is conjugate base of')
            chid = [i for i in chid if i not in metabolites]
            if len(chid) == 0:
                metabolites = True
                return completeChid, metabolites, charged
            else:
                for element in chid:
                    if element not in metabolites:
                        completeChid = self.client.service.getCompleteEntity(element)
                        if completeChid['charge'] != None:
                            return self.get_neutral_id(completeChid, metabolites, charged)
        elif int(charge) > 0:
            chid = self.get_parent_from_complete(completeChid, 'is conjugate acid of')
            chid = [i for i in chid if i not in metabolites]
            if len(chid) == 0:
                metabolites = True
                return completeChid, metabolites, charged
            else:
                for element in chid:
                    if element not in metabolites:
                        completeChid = self.client.service.getCompleteEntity(element)
                        if completeChid['charge'] != None:
                            return self.get_neutral_id(completeChid, metabolites, charged)


    def obtain_lowest_charge(self, visited_charges):
        '''
        Order the chIDs by their charge, giving priority to the
        positive numbers and closest to zero
        e.g. (+1, -1, +2, -2, ...)
        '''
        sorted_list = sorted(visited_charges, key=lambda x: abs(int(x[1])))
        completeChid = self.client.service.getCompleteEntity(sorted_list[0][0])
        return completeChid


    def get_parent_from_complete(self, completeChid, ontology_type, deleted_entities = []):
        '''
        Get the chebiId of the parent metabolite that has
        the described ontology type.
        '''
        ontology = []
        parent_ontology = completeChid['OntologyParents']
        for ont in parent_ontology:
            if ont['type'] == ontology_type:
                if ont['chebiId'] not in deleted_entities:
                    ontology.append(ont['chebiId'])
        return ontology


    def get_general_id(self, chid):
        '''
        Get the most general ID possible in order to avoid steroisomers
        '''
        forbidden = ['CHEBI:33704']
        ontology = self.get_parent_from_complete(chid, 'is a')
        if len(ontology) >= 1:
            new_chid = self.client.service.getCompleteEntity(ontology[0])
            if (new_chid['charge'] == None or new_chid['chebiId'] in forbidden):
                return chid
            else:
                return self.get_general_id(new_chid)
        else:
            return chid


    def get_specific_id(self, chid):
        '''
        Get the most general ID possible in order to avoid steroisomers
        '''
        ontology = self.get_child_from_complete(chid, 'is a')
        if len(ontology) >= 1:
            new_chid = self.client.service.getCompleteEntity(ontology[0])
            if new_chid['charge'] == None:
                return self.get_specific_id(new_chid)
            else:
                return new_chid
        else:
            return chid


    def get_child_from_complete(self, completeChid, ontology_type, deleted_entities = []):
        '''
        Get the chebiId of the parent metabolite that has
        the described ontology type.
        '''
        ontology = []
        child_ontology = completeChid['OntologyChildren']
        for ont in child_ontology:
            if ont['type'] == ontology_type:
                if ont['chebiId'] not in deleted_entities:
                    ontology.append(ont['chebiId'])
        return ontology



if __name__ == '__main__':
    
    chebi = chebi_from_string()
    chebi.chebi_connect()
    #a, b = chebi.chebi('NAD')
    #print a, b

#    a, b = chebi.chebi('aminopropylcadaverine')
#    print a, b

#   a = chebi.obtain_ontology_chid('CHEBI:64858')   
#    print a
    print chebi.chebi("4'-phosphopantetheine")
    
    print chebi.chebi("pantetheine 4'-phosphate")
    sys.exit()

    chid = chebi.obtain_ontology_chid('CHEBI:15763')
    print chid
    sys.exit()

    #chid = chebi.obtain_ontology_chid('CHEBI:15980')
    #chid = chebi.obtain_ontology_chid('CHEBI:29356')
    #chid = chebi.obtain_ontology_chid('CHEBI:26806')
    chid = chebi.obtain_ontology_chid('CHEBI:18009')
    #chid = chebi.obtain_ontology_chid('CHEBI:15468')
    chid = chebi.obtain_ontology_chid('CHEBI:15570')



    #chid = chebi.obtain_ontology_chid('CHEBI:18009')
    print chid
    sys.exit()



    print 'Pantoic acid'
    ont = chebi.obtain_chebi_parents('CHEBI:14737', '')
    print ont
    print
    print 'R-pantoate'
    ont = chebi.obtain_chebi_parents('CHEBI:15980', 'is conjugate base of')
    print ont
    print
    print 'R-pantoic acid'
    ont = chebi.obtain_chebi_parents('CHEBI:18697', '')
    print ont

    ont = chebi.obtain_chebi_parents('CHEBI:16027')
    print '\n\n'
    print ont

    ont = chebi.obtain_chebi_parents('CHEBI:40721')
    print '\n\n'
    print ont

    ont = chebi.obtain_chebi_parents('CHEBI:37096')
    print '\n\n'
    print ont

    print '------------------------------------------'
    print 
    print chebi.client.service.getCompleteEntity('CHEBI:16027')
    print '------------------------------------------'
    print 
    print chebi.client.service.getCompleteEntity('CHEBI:40721')
    print '------------------------------------------'
    print 
    print chebi.client.service.getCompleteEntity('CHEBI:37096')

    print '\n\n\n\n\n\n'

    print chebi.client.service.getCompleteEntity('CHEBI:14737')
    print '------------------------------------------'
    print 
    print chebi.client.service.getCompleteEntity('CHEBI:15980')
    print '------------------------------------------'
    print 
    print chebi.client.service.getCompleteEntity('CHEBI:18697')
    print '------------------------------------------'
    print 
    print chebi.client.service.getCompleteEntity('CHEBI:35972')