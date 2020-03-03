#!/usr/bin/env python

from neo4j.v1 import GraphDatabase, basic_auth
from chebi_from_string import chebi_from_string
from uniprot_queries import uniprot_queries

import re, sys


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


class graph_database():
    '''
    Create a Neo4j graph database, create the tools to fill it
    July, 5, 2017: add indexes
    '''
    def __init__(self, chebi, uniprot, brenda='',organism='',user='',pasw=''):
        self.organism = organism
        self.password = pasw
        self.user = user
        self.brenda = brenda
        self.chebi = chebi
        self.uniprot = uniprot
        self.EC2textname = {}

    def connect(self):
        '''
        Connect to the graph db using neo4j credentials
        '''
        self.driver = GraphDatabase.driver("bolt://localhost:7687",auth=basic_auth(self.user,self.password))
        self.session = self.driver.session()
        self.session.run("CREATE CONSTRAINT ON (n:Compound) ASSERT n.id IS UNIQUE")
        self.session.run("CREATE CONSTRAINT ON (n:Protein) ASSERT n.id IS UNIQUE")
        return None
    
    def extractsynonyms(self):
        dico_synonyms_prot={}
        syn=self.session.run("MATCH (n) RETURN n.synonyms as s,n.id as i")
        for r in syn:
            if not r['s']==None:
                for i in r['s']:
                    if len(i) > 3:
                        dico_synonyms_prot[clean(i)]=r['i']
        return dico_synonyms_prot

    
    # def create_protein(self, specie, uniID, uniName): #add by Cecile 08/08/2017
    #     '''
    #     Add a protein to the database
    #     '''
    #     torun="CREATE (a:Protein {specie: '%s', uniprotEntryName: '%s', id: '%s'})" % (specie, uniName, uniID)
    #     self.session.run(torun)
    #     return None
    
    def check_relationship_link(self, nodeA, nodeB, relationship): #add by Cecile 08/08/2017
        '''
        Check if a relationship already exists
        Do not check for the >
        '''
        result = self.session.run('MATCH (a)-[r:%s]-(b) '
                'WHERE a.id="%s" AND b.id="%s" RETURN a,r,b'%(relationship, nodeA, nodeB))
        if len([i for i in result]) > 0:
                is_rel = True
        else:
                is_rel = False
        return is_rel


    def check_species_relationship(self, species1, species2, relationship):
        '''
        Check relationship between species
        '''
        result = self.session.run('MATCH (a)-[r:%s]-(b) '
                'WHERE a.specie="%s" AND b.specie="%s" RETURN a,r,b'%(relationship, species1, species2))
        if len([i for i in result]) > 0:
                is_rel = True
        else:
                is_rel = False
        return is_rel
    

    def add_orthology_relationship(self, protein1, protein2, cluster):
        '''
        Add a orthology relationship between two db nodes.
        The relationship is added only if the two nodes are present.
        The specie is not necessary because uniprot IDs are species specific
        '''
        protein1_exists = self.check_protein(protein1)
        protein2_exists = self.check_protein(protein2)
        if (protein1_exists == True) and (protein2_exists == True):
            is_rel, r = self.check_relationship(protein1, protein2, 'Orthology_relationship')
            if is_rel == False:
                add_rel = 'MATCH (n), (y) WHERE n.id = "%s" and y.id = "%s"\
                MERGE (n)-[:Orthology_relationship {orthology_group: ["%s"]}]->(y)'%(protein1, protein2, cluster)
                self.session.run(add_rel)
            elif is_rel == True:
                updated_property = [cluster]
                for rel in r:
                    updated_property += rel['r']['orthology_group']
                add_rel = 'MATCH (n)-[r:Orthology_relationship]-(y) WHERE n.id = "%s" and y.id = "%s"\
                SET r.orthology_group = "%s"'%(protein1, protein2, list(set(updated_property)))
        return None


    # def fill_graph_db_ortho(self,orthologs): #add cecile 08/09/2017
    #     '''
    #     Fill the graph database with the ortholog relations stored in the orthologs dictionary
    #     Orthologs is a dictionary containing: ortho group; seq id; specie
    #     '''
    #     #print orthologs
    #     for og in orthologs:#og is the ortholog group
    #         print "OG "+og
    #         for seq in orthologs[og]: #seq is the sequence format tr|A9RJQ8|A9RJQ8_PHYPA\n
    #             print "seq "+seq
    #             print orthologs[og][seq]
    #             s1=seq.rstrip()
    #             tabseq=s1.split('|')
    #             specie1=orthologs[og][seq]
    #             specie1=re.sub('.fasta$','',specie1)
    #             specie1=re.sub('_',' ',specie1)
    #             for seq2 in orthologs[og]:
    #                 s2=seq2.rstrip()
    #                 tabseq2=s2.split('|')
    #                 specie2=orthologs[og][seq2]
    #                 specie2=re.sub('.fasta$','',specie2)
    #                 specie2=re.sub('_',' ',specie2)
    #                 if not seq == seq2: #not the two same proteins
    #                     ispres=self.check_protein(tabseq[1])
    #                     if not ispres:
    #                         #add the protein to the db
    #                         print "specie %s, ts1 %s, ts2 %s" %(specie1,tabseq[1],tabseq[2])
    #                         self.create_protein(specie1,tabseq[1],tabseq[2])
    #                     ispres2=self.check_protein(tabseq2[1])
    #                     if not ispres2:
    #                         print "specie2 %s, ts1 %s, ts2 %s" %(specie2,tabseq2[1],tabseq2[2])
    #                         self.create_protein(specie2,tabseq2[1],tabseq2[2])
    #                     if not self.check_relationship_link(tabseq[1],tabseq2[1],'Ortho_Inparanoid'):
    #                         print "OrthoRel ts1 %s, ts21 %s" %(tabseq[1],tabseq2[1])
    #                         self.create_orthorelationship(tabseq[1], tabseq2[1])
    #     return 'NA'
    

    def create_blastp_orthology_relationship(self, nodeA, nodeB):
        '''
        Create a Ortho_Inparanoid relationship between two existing nodes in the database
        Relationship in neo4j can only have one type!
        '''
        self.session.run("MATCH (a), (b) WHERE a.id={nodeA} AND b.id={nodeB} \
                CREATE (a)-[:Orthology_relationship]\
                ->(b)",{"nodeA":nodeA, "nodeB":nodeB})
        return None


    def check_specie(self):
        '''
        Check if the organism of interest is in the database
        '''
        #CECILE COMMENT: Do we want to have the taxid and the specie name in the DB?
        specie = self.session.run("MATCH (n) WHERE n.specie = {specie} RETURN n.specie as specie",
            {"specie":self.organism})
        species = [i['specie'].encode('utf-8') for i in specie]
        if len(species) > 0:
            if self.organism in species:
                is_spec = True
            else:
                is_spec = False
        else:
            is_spec = False
        return is_spec


    def obtain_uniIDs_from_specie(self, specie):
        '''
        obtain uniIDs from all proteins of a given specie
        '''
        uniIDs = []
        query = "MATCH (n) WHERE n.specie = '%s' RETURN DISTINCT n.uniprotID as uniID"%(specie)
        result = self.session.run(query)
        for elem in result:
            uniIDs.append(elem['uniID'])
        return uniIDs


    def create_enzyme(self, uniID, uniName, syns_list, ec):
        '''
        Add an enzyme entity to the database for the curent organism (in self.organism)
        '''
        self.session.run("CREATE (a:Enzyme:Protein {uniprotID: {uniID}, synonyms: {syns}, "
                "ECs: {ec}, id: {uniID}, uniprotEntryName: {uniEntry}, specie: {specie}})",
                {"uniID": uniID, "syns": syns_list, "ec": [ec], "uniEntry": uniName, "specie":self.organism})
        return None


    def set_enzyme_properties(self, uniID, syns_list, ecs):
        '''
        Use an existing node and change its variable properties
        '''
        self.session.run("MATCH (a) WHERE a.id = {uniID} SET a.synonyms = {synonyms}, a.ECs = {ECs}",
                {"uniID": uniID, "synonyms": syns_list, "ECs":ecs})
        return None


    def extract_species(self):
        '''
        '''
        result = self.session.run("MATCH (n) RETURN DISTINCT n.specie")
        result = [r['n.specie'] for r in result]
        return result

    def add_prop_enzyme(self, uniID, syns_list, ec):
        '''
        Use an existing enzyme and add new properties.
        '''
        result = self.session.run("MATCH (a) WHERE a.id = {uniID} RETURN a.synonyms as "
                "synonyms, a.ECs as EC", {"uniID":uniID})
        for n in result:
                synonyms = set(n["synonyms"])
                EC = set(n["EC"])
                for syn in syns_list:
                        synonyms.add(syn)
                EC.add(ec)
                self.set_enzyme_properties(uniID, list(synonyms), list(EC))
        return None


    def create_second_level_enzyme(self, uniID, synonyms, sentences, ECs, cluster):
        '''
        Create an enzyme fron the clustering step.
        This enzyme includes the names and the sentences where it has been found
        We also create an edge that connects the cluster to the enzymes that
        make that clusterin the first level graph
        '''
        cluster = 'proteinCluster' + str(cluster)
        prot_exist = self.check_compressed_protein(uniID)
        if prot_exist == False:
            self.session.run("CREATE (a:Enzyme:Protein:Compressed {uniprotIDs: {uniID}, synonyms: {syns}, "
                    "sentences: {sentences}, ECs: {ECs}, id: {cluster}})",
                    {"uniID": [uniID], "syns": synonyms, "sentences": sentences, "ECs": ECs, "cluster": cluster})
            #Create link with lower level
            self.session.run("MATCH (n:Protein:Compressed), (y:Protein) WHERE {uniID}\
                IN n.uniprotIDs AND y.id = {uniID} CREATE (n)-[r:to_compressed]->(y)",
              {"uniID": uniID})
        else:
            # In case a node is repeated, or the script is run twice
            self.modify_second_level_enzyme(uniID, uniID, sentences, synonyms, ECs)
        return None


    def modify_second_level_enzyme(self, old_uniID, new_uniID, sentences, synonyms, ECs):
        '''
        Modify an existing second level enzyme, 
        adding new proteins to the 'cluster'
        '''
        result = self.session.run("MATCH (n:Protein:Enzyme:Compressed) WHERE {uniID} IN n.uniprotIDs\
         RETURN n.uniprotIDs AS uniIDs, n.synonyms AS nsyns, n.sentences AS nSents, n.ECs AS nECs",{"uniID": old_uniID})

        for r in result:
            uniIDs = r['uniIDs']; nsentences = r['nSents']; syns = r['nsyns']; nECs = r['nECs']
            if new_uniID not in uniIDs:
                uniIDs.append(new_uniID)
            for sentence in sentences:
                if sentence not in nsentences:
                    nsentences.append(sentence)
            for syn in synonyms:
                if syn not in syns:
                    syns.append(syn)
            for EC in nECs:
                if EC not in ECs:
                    ECs.append(EC)

        self.session.run("MATCH (n:Protein:Compressed) WHERE {uniID} IN n.uniprotIDs\
         SET n.synonyms = {synonyms}, n.uniprotIDs = {uniIDs}, n.sentences = {sentences}, n.ECs = {ECs}\
          WITH n MATCH (y:Protein) WHERE y.id = {new_uniID} MERGE (n)-[r:to_compressed]->(y)",
        {"uniID": old_uniID, "uniIDs": uniIDs, "synonyms": syns, "sentences": nsentences, "new_uniID": new_uniID, "ECs": ECs})
        return None


    def create_compound(self, chebiID, chebiName):
        '''
        Add a compound to the database
        '''
        self.session.run("CREATE (a:Compound {chebiID: {chebiID}, "
                "compoundName: {compoundName}, id: {chebiID}})",
                {"chebiID": chebiID, "compoundName": chebiName})
        return None


    def add_protein(self, uniprotID, uniprotName, uniGenes, uniProteins):
        '''
        Add the protein to the neo4j database as a protein instance
        '''
        is_node = self.check_protein(uniprotID)
        if is_node == False:
            protList = []; geneList = []
            for prot in uniProteins.split('('):
                protList.append(re.sub('\)', '', prot).encode('utf8'))
            for gene in uniGenes.split():
                geneList.append(re.sub(':', '', gene).encode('utf8'))
            query = "CREATE (n:Protein {uniprotID:'%s', id:'%s', uniProtEntryName: '%s',\
             uniprotGenesNames: %s, uniprotProteinNames: %s, specie: '%s'})"%\
             (uniprotID, uniprotID, uniprotName, geneList, protList, self.organism)
            self.session.run(query)
        return None


    def create_second_level_compound(self, chebiID, compNames, sentences, cluster):
        '''
        Create a compound for the clustering step
        These compounds are a set of metabolites that cluster together
        Create relationship with the first level graph corresponding compounds
        '''
        cluster = 'compoundCluster' + str(cluster)
        comp_exist = self.check_compressed_compound(chebiID)
        if comp_exist == False:
            self.session.run("CREATE (a:Compound:Compressed {chebiIDs: {chebiID},\
                compoundNames: {compNames}, sentences: {sentences}, id: {cluster}})",
                    {"chebiID": [chebiID], "compNames": compNames, "sentences": sentences, "cluster": cluster})

            self.session.run("MATCH (n:Compound:Compressed), (y:Compound) WHERE {chebiID}\
                IN n.chebiIDs AND y.id = {chebiID} CREATE (n)-[r:to_compressed]->(y)",
              {"chebiID": chebiID})
        else:
            # In case a node is repeated, or the script is run twice
            self.modify_second_level_compound(chebiID, chebiID, compNames, sentences)
        return None

    def modify_second_level_compound(self, chebiID_old, chebiID_new, compNames, sentences):
        '''
        Modify an existing second level compound, 
        adding new compounds to the 'cluster'
        '''
        result = self.session.run("MATCH (n:Compound:Compressed) WHERE {chebiID} IN n.chebiIDs\
         RETURN n.chebiIDs AS chebiIDs, n.compoundNames AS cNames, n.sentences AS nSents",{"chebiID": chebiID_old})

        for r in result:
            chebiIDs = r['chebiIDs']; nsentences = r['nSents']; cNames = r['cNames']
            if chebiID_new not in chebiIDs:
                chebiIDs.append(chebiID_new)
            for sentence in sentences:
                if sentence not in nsentences:
                    nsentences.append(sentence)
            for name in compNames:
                if name not in cNames:
                    cNames.append(name)

        self.session.run("MATCH (n:Compound:Compressed) WHERE {chebiID} IN n.chebiIDs\
         SET n.compoundNames = {compNames}, n.chebiIDs = {chebiIDs}, n.sentences = {sentences}\
          WITH n MATCH (y:Compound) WHERE y.id = {new_chebiID} MERGE (n)-[r:to_compressed]->(y)",
        {"chebiID": chebiID_old, "chebiIDs": chebiIDs, "compNames": cNames, "sentences": nsentences, "new_chebiID": chebiID_new})
        return None


    def create_compressed_relationship(self, nodeA, nodeB, nodeA_type, nodeB_type, sentences):
        '''
        Create a compressed relatinship between two nodes,
        these two nodes must be from the compressed level of the graph
        '''
        if nodeA_type == 'Compound':
            nodeA_prop = 'chebiIDs'
        elif nodeA_type == 'Protein':
            nodeA_prop = 'uniprotIDs'
        if nodeB_type == 'Compound':
            nodeB_prop = 'chebiIDs'
        elif nodeB_type == 'Protein':
            nodeB_prop = 'uniprotIDs'
        if nodeA_type == None or nodeB_type == None:
            return
        result = self.session.run('MATCH (n:Compressed)-[r]->(y:Compressed)\
         WHERE "%s" IN n.%s AND "%s" IN y.%s RETURN r'%(nodeA, nodeA_prop, nodeB, nodeB_prop))
        if len([r for r in result]) == 0:
            self.session.run('MATCH (n:Compressed), (y:Compressed) WHERE "%s" IN n.%s\
             AND "%s" IN y.%s MERGE (n)-[r:compressed_relationship {sentences: %s}]-(y)'
             %(nodeA, nodeA_prop, nodeB, nodeB_prop, list(sentences)))
        else:
            for r in result:
                for sent in r["sentences"]:
                    sentences.add(sent)
            self.session.run('MATCH (n:Compressed), (y:Compressed) WHERE "%s" IN n.%s AND "%s" IN y.%s\
             MERGE (n)-[r:compressed_relationship {sentences: %s}]-(y)'
             %(nodeA, nodeA_prop, nodeB, nodeB_prop, list(sentences)))
        return None


    def create_brenda_relationship(self, nodeA, nodeB, ec, reaction, specie):
        '''
        Create a relationship between to existing nodes in the database
        '''
        self.session.run('MATCH (a), (b) WHERE a.id={nodeA} AND b.id={nodeB} '
                "CREATE (a)-[:Brenda_relationship {ECs: {ec}, reactionsBrenda: {reaction}, species: {specie}}]"
                "->(b)",{"nodeA":nodeA, "nodeB":nodeB, "ec":[ec], "reaction":[reaction], "specie":[specie]})
        return None


    def check_relationship(self, nodeA, nodeB, relationship):
        '''
        Check if a relationship already exists
        '''
        result = self.session.run('MATCH (a)-[r:%s]->(b) '
                'WHERE a.id="%s" AND b.id="%s" RETURN a,r,b'%(relationship,nodeA, nodeB))
        result = [i for i in result]
        if len(result) > 0:
                is_rel = True
        else:
                is_rel = False
        return is_rel, result


    def update_brenda_relationship(self, nodeA, nodeB, ec, reaction, specie):
        '''
        Update the information from ane xisting relationship
        ec, reaction and specie must be a list
        '''
        result = self.session.run('MATCH (a)-[r:Brenda_relationship]->(b) '
                'WHERE a.id={nodeA} AND b.id={nodeB} '
                'RETURN r.ECs as ecs, r.reactionsBrenda as reac, r.species as species',
                {"nodeA":nodeA, "nodeB":nodeB})
        for n in result:
                ecs = (set(n["ecs"]+ec))
                reactions = (set(n["reac"])|set(reaction))
                species = (set(n["species"]+specie))
                self.set_reaction_properties(nodeA, nodeB, list(ecs), list(reactions), list(species))
        return None


    def set_reaction_properties(self, nodeA, nodeB, ec, reaction, specie):
        '''
        Create a relationship between to existing nodes in the database
        all properties musy be list type
        '''
        self.session.run('MATCH (a)-[r:Brenda_relationship]->(b) WHERE a.id={nodeA} AND b.id={nodeB} '
                "SET r.ECs = {ec}, r.reactionsBrenda = {reaction}, r.species = {specie}",
                {"nodeA":nodeA, "nodeB":nodeB, "ec":ec, "reaction":reaction, "specie":specie})
        return None


    def create_prop_relationship(self, nodeA, nodeB, relationship, sentence):
        '''
        Create the relationship and add the sentence information
        '''
        self.session.run('MATCH (a), (b) WHERE a.name="%s" AND b.name="%s" '
                "CREATE (a)-[:%s {Sentence: %s}]->(b)"%(nodeA, nodeB, relationship, sentence))
        return None


    def check_protein(self, ID):
        '''
        Check if a node exists. Will be used before adding a node
        '''
        result = self.session.run("MATCH (a:Protein) WHERE a.id={id} "
                "RETURN a.id AS id", 
                {"id": ID})
        if len([i for i in result]) > 0:
                is_node = True
        else:
                is_node = False
        return is_node


    def check_compound(self, ID):
        '''
        Check if a node exists. Will be used before adding a node
        '''
        result = self.session.run("MATCH (a:Compound) WHERE a.id={id} "
                "RETURN a.id AS id", 
                {"id": ID})
        if len([i for i in result]) > 0:
                is_node = True
        else:
                is_node = False
        return is_node


    def check_compressed_protein(self, ID):
        '''
        Check if a protein already exists for the compressed graph
        '''
        result = self.session.run("MATCH (a:Compressed:Protein) WHERE {id} IN a.uniprotIDs\
                 RETURN a.uniprotIDs AS id", 
                {"id": ID})
        if len([i for i in result]) > 0:
            is_node = True
        else:
            is_node = False
        return is_node


    def check_compressed_compound(self, chebiID):
        '''
        Check if a compound already exists for the compressed graph
        '''
        result = self.session.run("MATCH (a:Compressed:Compound) WHERE {id} IN a.chebiIDs\
                 RETURN a.chebiIDs AS id", 
                {"id": chebiID})
        if len([i for i in result]) > 0:
            is_node = True
        else:
            is_node = False
        return is_node


    def is_enzyme(self, syn):
        '''
        Check if the entity is into the synonyms of any enzyme
        '''
        result = self.session.run("MATCH (a) WHERE '%s' IN a.Synonyms "
                "RETURN a.name, labels(a)"%(syn))
        return result


    def is_compound(self, compound):
        '''
        Check if the entity is included within any compound 
        (node) name
        '''
        result = self.session.run('MATCH (a) WHERE a.name = "%s" '
                'RETURN a.name, labels(a)'%(compound))
        return result


    def brenda_obtain_EC(self):
        '''
        Obtain the EC numbers belonging to the desired organism
        '''
        ecNumbers = self.brenda.run_function('getEcNumbersFromOrganism', organism=self.organism)
        ecNumbers = ecNumbers.split('!')
        return ecNumbers


    def parse_brenda_output(self, input_string, entity_type):
        '''
        Brenda output is a complex string, here the string is parsed
        and a list with all the relevant elements is outputed
        '''
        entity_set = set()
        if input_string != '':
                for elem in input_string.split('!'):
                        reaction = None
                        for entity in elem.split('#'):
                                if entity.startswith(entity_type):
                                        entry = entity.split('*')[1]
                                if entity.startswith('reactionPartners'):
                                        reaction = entity.split('*')[1]
                        if entity_type in ('substrate', 'product'):
                                entity_set.add((entry,reaction))
                        else:
                                entity_set.add(entry)
        return entity_set



    def fill_graph_db(self, ec, substrates, products, synonyms, ec2uniprot):
        '''
        Use the Brenda information to fill the graph database, create a node
        for each protein that correspond to the enzyme and establish the relationships 
        according to the connection by compound
        '''
        for prot in ec2uniprot[ec]:
                uniID = prot[0]; uniName = prot[1]
                is_node = self.check_protein(uniID)
                if is_node == True:
                        self.add_prop_enzyme(uniID, list(synonyms), ec)
                else:
                    self.create_enzyme(uniID, uniName, list(synonyms), ec)

                for subs in substrates:
                        substrate = subs[0]; reaction = subs[1]
                        chid, name = self.chebi.chebi(substrate)
                        if chid != None:
                                is_node = self.check_compound(chid)
                                if is_node == False:
                                        self.create_compound(chid, name)
                                is_rel, r = self.check_relationship(chid, uniID, 'Brenda_relationship')
                                if is_rel == False:
                                        self.create_brenda_relationship(chid, uniID, ec, reaction, self.organism)
                                else:
                                        self.update_brenda_relationship(chid, uniID, [ec], [reaction], [self.organism])

                for prod in products:
                        product = prod[0]; reaction = prod[1]
                        chid, name = self.chebi.chebi(product)
                        if chid != None:
                                is_node = self.check_compound(chid)
                                if is_node == False:
                                        self.create_compound(chid, name)
                                is_rel, r = self.check_relationship(uniID, chid, "Brenda_relationship")
                                if is_rel == False:
                                        self.create_brenda_relationship(uniID, chid, ec, reaction, self.organism)
                                else:
                                        self.update_brenda_relationship(uniID, chid, [ec], [reaction], [self.organism])
        return None


    def create_database(self, ec2uniprot):
        '''
        Just as it was done for the MySQL database, use the BRENDA database
        information to fill the Neo4j database. The main difference with
        the mySQL database is that here we add the relationships as well
        '''
        process = False
        self.connect()
        self.EC2textname_dictio('EC2entryName.txt')
        ecNumbers = self.brenda_obtain_EC()
        for EC in ecNumbers:
            print EC
#            if EC == '2.7.11.15':
#                process = True
#            result = self.session.run('MATCH (n) WHERE n.ECs = "%s" RETURN n'%(EC))
#            if len([r for r in result]) == 0:
#            if process == True:
            substrates = self.brenda.run_function('getSubstrate', organism=self.organism, ecNumber=EC)
            substrates = self.parse_brenda_output(substrates, 'substrate')
            products = self.brenda.run_function('getProduct', organism=self.organism, ecNumber=EC)
            products = self.parse_brenda_output(products, 'product')
            synonyms = self.brenda.run_function('getSynonyms', organism=self.organism, ecNumber=EC)
            synonyms = self.parse_brenda_output(synonyms, 'synonyms')
            recommendedName = self.brenda.run_function('getRecommendedName', ecNumber=EC)
            recommendedName = self.parse_brenda_output(recommendedName, 'recommendedName')
            synonyms = synonyms|recommendedName
            try:
                ec2uniprot[EC]
                self.fill_graph_db(EC, substrates, products, synonyms, ec2uniprot)
            except KeyError:
                try:
                    gene_name = self.EC2textname[EC]
                    uniIDs = self.uniprot.query_id(gene_name)
                    if uniIDs != ('','','',''):
                        ec2uniprot[EC] = [(uniIDs[0], uniIDs[1])]
                        self.fill_graph_db(EC, substrates, products, synonyms, ec2uniprot)
                except KeyError:
                    pass
        return None

    def EC2textname_dictio(self, EC2namefile):
        '''
        '''
        for line in open(EC2namefile).read().split('\n')[:-1]:
            EC_number = line.split('\t')[0].strip()
            description = line.split('\t')[1]
            self.EC2textname[EC_number] = description
        return None


    def check_db_length(self):
        '''
        Check the length of the graph db, this will help to take the
        decission of creating it or using the existing db
        '''
        nodes = None
        driver = GraphDatabase.driver("bolt://localhost:7687", 
                auth=basic_auth(self.user,self.password))
        session = driver.session()
        result = session.run('MATCH (n) RETURN count(n)')
        for r in result:
                for n in r:
                        nodes = r[0]
        session.close()
        return nodes


    def check_TM_relationship(self, ent1, ent2):
        '''
        Check if a given relationship exists
        '''
        result = self.session.run('MATCH (n)-[r]-(y)' 
                ' WHERE n.name = "%s" AND y.name = "%s" RETURN n.name, y.name'%(ent1, ent2))
        res = [i for i in result]
        return res


    def add_TM_check(self, key, value, graph_entity):
        '''
        Iterate over the list of keys and values, check if each of 
        these relationships exist, if the don't add them to the
        Neo4j graph DB
        '''
        if key in graph_entity.keys():
                myKeys = graph_entity[key]
                for v in value.keys():
                        sentence = set(value[v]['sentence'])
                        if v in graph_entity.keys():
                                myVals = graph_entity[v]
                                for element in myKeys:
                                        for target in myVals:
                                                rel = self.check_TM_relationship(element, target)
                                                if len(rel) == 0:
                                                        #Modify this statement to include properties to the
                                                        #relationship, this will allow the addition of the text
                                                        #Also need to add the type of relationship (binding, negative, ...)
                                                        self.create_prop_relationship(element, target, 'TM_relationship', sentence)
                #Maybe we need to include a property that states that
                #for the same relationship a lot of different entities 
                #(synonyms) were used
        return None



    def add_TM_relationship(self, cause_dictio, non_cause_dictio, graph_entity):
        '''
        Add the relationships gathered from the TEES text mining
        into the Neo4j graph DB.
        '''
        for key, value in cause_dictio.iteritems():
                self.add_TM_check(key, value, graph_entity)
        for key, value in non_cause_dictio.iteritems():
                self.add_TM_check(key, value, graph_entity)
        return None


    def extract_pattern(self, cypher_query):
        '''
        Extract a cypher query from the neo4j database
        '''
        result = self.session.run(cypher_query)
        result = [i for i in result]
        return result

if __name__ == '__main__':

    gd = graph_database('', '', '', 'Homo sapiens', 'neo4j', 'neo4j')
    gd.connect()

    gd.extractsynonyms()
    sys.exit()

    gd.obtain_uniIDs_from_specie('Citrus Clementina')

    #print dico_syns