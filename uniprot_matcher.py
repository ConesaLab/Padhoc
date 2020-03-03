#!/usr/bin/env python

import urllib, urllib2

class uniprot_matcher():
    '''
    Match a protein to the UniProt database to check if it
    is a "real" protein, and to obtain the UniProt identifier
    '''

    def __init__(self, organisms_file):

        self.organisms_file = organisms_file
        self.ec2Uniprot_dict = {}
        self.uniprot_dictio = {}

    def uniprot_accession(self, inp, out, query):
        '''
        Protocol to access the uniprot database.
        '''
        if query not in self.uniprot_dictio.keys():
            url = 'http://www.uniprot.org/uploadlists/'
            params = {
                    'from': inp,
                    'to': out,
                    'format':'tab',
                    'query': query
            }
            data = urllib.urlencode(params)
            request = urllib2.Request(url, data)
            response = urllib2.urlopen(request)
            page = response.read(200000)
            self.uniprot_dictio[query] = page
        else:
            page = self.uniprot_dictio[query]
        return page


    def uniprotOrganism(self, taxa):
        '''
        Convert from organism to UniProt organism using the
        conversion file. For this function we need the 
        organism taxonomic number
        '''
        uni_name = False
        for line in self.organisms_file.split('\n'):
                line = line.split()
                try:
                        if line[2][:-1]== taxa:
                                uni_name = line[0]
                except IndexError:
                        pass
        if uni_name == False:
            uni_name = 'NoOrganismMatch'
        return uni_name


    def fill_ec_dictio(self, ec, uniprotID, uniprotGene):
        '''
        Fill the dictionary of ECs to uniprot identifiers
        '''
        if ec in self.ec2Uniprot_dict.keys():
                self.ec2Uniprot_dict[ec].append((uniprotID, uniprotGene))
        else:
                self.ec2Uniprot_dict[ec] = [(uniprotID, uniprotGene)]
        return None

    
    def EC2Uniprot(self, converter, taxa):
        '''
        Use the expasy file to convert from EC to UniProt.
        Filter by organism using the uniprot organisms file
        '''
        organism = self.uniprotOrganism(taxa)
        for line in converter.split('//'):
                if line.strip().startswith('ID'):
                        line = line.split('\n')
                        ec_id = line[1].split()[1]
                        for l in line:
                                if l.startswith('DR'):
                                        l = l.replace('DR', '').replace(' ', '').split(';')
                                        for uniprot in l[:-1]:
                                                uniprotID = uniprot.split(',')[0]
                                                uniprotGene = uniprot.split(',')[1]
                                                if uniprotGene.split('_')[-1] == organism:
                                                        self.fill_ec_dictio(ec_id, uniprotID, uniprotGene)
        return None


if __name__ == '__main__':

    uniprot = uniprot_matcher('')
    a = uniprot.uniprot_accession('GENENAME', 'ACC', 'T11')
    print a