#make one py file to plot up a pathway
#modified from code use in NB projecct
#KLongnecker, 17 April 2017
import pandas as pd
import os
import re
import matplotlib.pyplot as plt
import palettable as pal
import glob

from Bio import SeqIO
from Bio.KEGG.REST import *
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
from IPython.display import Image, HTML
   
from IPython.core.debugger import Tracer 
    
def gatherDetails(usePathway,folderName,CO_fromTSQ,CO_values):
    #check if the directories exist, one for pathway files
    if not os.path.exists(folderName):
        os.makedirs(folderName)
    #else:
        #raise ValueError('Krista - be careful, this folder already exists')
                   
    #for item in trimPath: #searching within one pathway at a time
    #only one pathway at a time
    genes = getKfrom_ko(usePathway)
    compounds = getCfrom_ko(usePathway)

    #have to track genes and compounds differently for the biopython plotting later on 
    setG = set(genes)
    setC = set(compounds)
    setT = set(CO_fromTSQ)
    #figure out which compounds are in the TSQ data...
    intCompounds = setC.intersection(setT)
        
    ## plot the pathway map for this pathway, get details from KEGG for plotting
    ###ideally I will make the colors a function of amount
    useColors = pal.colorbrewer.qualitative.Set1_4.hex_colors
    size = 20 #turns out I can increase the size of the compounds in the plots

    #useColors.insert(0,'#f7f7f7') ## insert white at beginning
    pathway = KGML_parser.read(kegg_get(usePathway, "kgml"))
    for element in pathway.orthologs:
        #print element.name
        for graphic in element.graphics: #[only color if in Prochlorococcus]
            tg = element.name[3:9] #skip over the 'ko:'
            if (tg in setG):
                #in the pathway for Pro
                graphic.bgcolor = useColors[2] #
 
        # Change the colours of compounds (mostly same as genes
        for element in pathway.compounds:
            for graphic in element.graphics:
                tc = element.name[4:10] #skip over the 'cpd:'
                if (tc in intCompounds):
                    #in the pathway AND in the set for this particular K means group
                    graphic.bgcolor = useColors[3] #
                    graphic.width = size
                    graphic.height = size
                        
    canvas = KGMLCanvas(pathway, import_imagemap=True)
    pdfName = 'mapWithColors_' + str(usePathway) + '.pdf'
    #Tracer()()
    canvas.draw(folderName + '/' + pdfName)
    pdfName = None #empty it in case that is where I am having issues


#set up a function to get the list of K orthologues for a given pathway (must be defined as ko00140 NOT map00140)
def getKfrom_ko(ko_id):
    pathway_file = kegg_get(ko_id).read()  # query and read the pathway
    K_list = []

    current_section = None
    for line in pathway_file.rstrip().split("\n"):
        section = line[:12].strip()  # section names are within 12 columns
        if not section == "":
            current_section = section
        if current_section == "ORTHOLOGY":
            K_identifiers = line[12:].split("; ")
            t = K_identifiers[0]
            K_id = t[0:6]

            if not K_id in K_list:
                K_list.append(K_id)
    return K_list

#set up a function to get the list of compounds for a given pathway (must be defined as ko00140 NOT map00140)
def getCfrom_ko(ko_id):
    pathway_file = kegg_get(ko_id).read()  # query and read the pathway
    compound_list = []

    current_section = None
    for line in pathway_file.rstrip().split("\n"):
        section = line[:12].strip()  # section names are within 12 columns
        if not section == "":
            current_section = section
        if current_section == "COMPOUND":
            compound_identifiers = line[12:].split("; ")
            t = compound_identifiers[0]
            compound_id = t[0:6]

            if not compound_id in compound_list:
                compound_list.append(compound_id)
    return compound_list

# A bit of code that will help us display the PDF output
def PDF(filename):
    return HTML('<iframe src=%s width=700 height=350></iframe>' % filename)
