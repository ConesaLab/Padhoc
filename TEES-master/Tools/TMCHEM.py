from __future__ import print_function
parse__version__ = "$Revision: 1.3 $"

import sys,os
import time, datetime
import sys
import shutil
import subprocess
import tempfile
import codecs

try:
    import xml.etree.cElementTree as ET
except ImportError:
    import cElementTree as ET

sys.path.append(os.path.dirname(os.path.abspath(__file__))+"/..")
import Utils.ElementTreeUtils as ETUtils
import Utils.Settings as Settings
import Utils.Download as Download
from Tools import Tool

#from __future__ import print_function
from xml.dom.minidom import parse
from optparse import OptionParser
import xml.dom.minidom
from subprocess import call
import glob

def test(progDir):
    return True

def run(input,output=None,rep="/home/salva/tmChem/"):
    #input is an elementTree object!
#    print("Loading corpus "+input,file=sys.stderr)

    outtmp=tempfile.mkdtemp()

    print("tmchem work directory at "+outtmp,file=sys.stderr)

 #   os.system("cp "+input+" "+input+"_beforetmchem")#save the previous version of the input file
 #   if input.endswith(".tar.gz"): #uncompress the input file if needed
 #       os.system("tar -xvzf "+input)

    # 1) xml to pubtator format

    # Open XML document using minidom parser
    #DOMTree = xml.dom.minidom.parse(input)
    #collection = DOMTree.documentElement
    # Get all the documents in the collection
    #abstracts = collection.getElementsByTagName("document")

    # pubtator output file: open the file

    root=input.getroot()
    outpubtator=outtmp+'/sentences.pubtator'
    print("sentencepubtator file :"+outpubtator,file=sys.stderr)
    f = open(outpubtator,'w')

    # foreach document, print each sentence
#    for abstract in abstracts:
    for abstract in root.findall("document"):
        sentences = abstract.findall('sentence')
        for sent in sentences:
#		print("sent in sentences")
#		print(sent.attrib['id'],file=sys.stderr)
            f.write(sent.attrib["id"]+"|t|"+sent.attrib['text']+"\n"+sent.attrib["id"]+"|a|\n\n")
#	    print(sent.attrib["id"]+"|t|"+sent.attrib['text']+"\n"+sent.attrib["id"]+"|a|\n",file=sys.stderr)
    f.close()
   
    cwd = os.getcwd()
    print(cwd)

    # 2) run tmChem
    outtmptmchem=outtmp
    prevdir=os.getcwd()#save previous directory
    os.chdir(rep)#go to tmchem directory
    print("perl tmChem.pl -i "+outtmp+" -o "+outtmptmchem,file=sys.stderr)
    os.system("perl tmChem.pl -i "+outtmp+" -o "+outtmptmchem)

    #TOTEST cecile
    #head, tail = os.path.split(output) #separe repository and name of the file
    #print(head)
    #print(outtmptmchem+"/sentences.pubtator.tmChem")
    shutil.copy2(outtmptmchem+"/sentences.pubtator.tmChem",cwd) #COPY THE OUTPUT IN THE OUTREPWITH THE REST OF THE FILES
    print(cwd + '\n\n')
    sys.exit()
    #end to test cecile

    os.chdir(prevdir)#go back to the previous directory
    
	
    # 3) convert pubtator to xml
    vale=100
    ftmchem=glob.glob(outtmptmchem+"/*tmChem")
    # output pubtator file ftmchem[0]
    print("output tmchem file in pubtator format "+ftmchem[0],file=sys.stderr)
    sentencemetabolite={} #dictionnaire

    #add the chemicals to the xml file
    f=open(ftmchem[0],'r')

    nb=0
    


    for line in f:
        #read line by line the pubtator file
        line=line.rstrip("\r\n")
        words = line.split('\t')
        if 't' in words[0]:
            idtmp=words[0].split('|')
            sentencemetabolite[idtmp[0]]={}
        else:
            if len(words)>1:
                #not the line with the sentence 
                #dictionnary with sentence id, position debut, position fin, metabolite name
                se=str(nb)+'-'+words[1]+'-'+words[2]
                sentencemetabolite[words[0]][se]=words[3]
                nb=nb+1
    f.close()

    #print(input)
    tree=input
    root=tree.getroot()

    e=vale
    for sentence in root.findall("./document/"):#foreach sentence
        if sentence.get('id') in sentencemetabolite:
            sent=sentencemetabolite[sentence.get('id')]
            e=vale
            for posi in sorted(sent):
                tmpposi=posi.split('-')
                sentp1=int(float(tmpposi[2]))+1
                sentp1=str(sentp1)
                cos=tmpposi[1]+'-'+sentp1
                ide=sentence.get('id')+'.e'+str(e)
                obo=tmpposi[1]+'-'+tmpposi[2]
                #add node to the good sentence
                newentity=xml.etree.ElementTree.SubElement(sentence,tag='entity',attrib={"charOffset":cos,"given":"True","id":ide,"origBANNEROffset":obo,"source":"TmChem","text":sent[posi],"type":"Metabolite"})
                e=e+1
    #Write the output
    tree.write(output,encoding="UTF-8",method="xml")


if __name__=="__main__":
    import sys

    from optparse import OptionParser, OptionGroup
    # Import Psyco if available
    try:
        import psyco
        psyco.full()
        print("Found Psyco, using",file=sys.stderr)
    except ImportError:
        print("Psyco not installed",file=sys.stderr)

    optparser = OptionParser(description="TmChem named entity recognizer wrapper")
    optparser.add_option("-i", "--input", default=None, dest="finput", help="The xml file produce after BANNER (gene recognition), corpus interaction XML format", metavar="FILE")
    optparser.add_option("-t", "--outtmp", default="tmpTMCHEM", dest="outtmp", help="The temporary repository for pubtator format (create from the xml file produce after BANNER)")
    optparser.add_option("-o", "--output", default=None, dest="foutput", help="Output file in Interaction XML format.")
    optparser.add_option("-r", "--reptmchem", default="/home/vagrant/tmChem/", dest="rep", help="repertory with the tmchem program")
#file=sys.stderr    optparser.add_option("-p", "--processElement", default="sentence", dest="processElement", help="input element tag (usually \"sentence\" or \"document\")")
#    optparser.add_option("-s", "--split", default=False, action="store_true", dest="splitNewlines", help="Split BANNER entities at newlines")
#    optparser.add_option("--debug", default=False, action="store_true", dest="debug", help="Preserve temporary working directory")


    (options, args)=optparser.parse_args()

    run(options.finput,options.foutput,options.rep)


