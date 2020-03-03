#!/usr/bin/env python

import hashlib
from SOAPpy import WSDL, SOAPProxy, Errors
from socket import error as socket_error
from xml.sax import SAXParseException

class brenda_annotation():
    '''
    Use Brenda database to annotate the entities found by the 
    NER algorithm
    '''
    def __init__(self,user='',pswrd=''):
            self.user = user
            self.pswrd = pswrd

    def access_protocol(self):
        '''
        Define the accession protocol information 
        '''
        self.endpointURL = "https://www.brenda-enzymes.org/soap/brenda_server.php"
        self.password = hashlib.sha256(self.pswrd).hexdigest()
        self.client = SOAPProxy(self.endpointURL)
        return None

    def run_function(self, function, ecNumber=None, organism=None, recommendedName=None):
        '''
        Run the brenda function of interest with the required parameters
        '''
        search = ","
        resultString = ''
        if ecNumber != None:
                search += ('ecNumber*' + ecNumber + '#')
        if recommendedName != None:
                search += ('recommendedName*' + recommendedName + '#')
        if organism != None:
                search += ('organism*' + organism + '#')
        parameters = "%s,"%(self.user)+self.password+search
        for n in range(10):
                try:
                        resultString = getattr(self.client,function)(parameters)
                        break
                except (socket_error, SAXParseException, Errors.HTTPError):
                        pass
        return resultString
