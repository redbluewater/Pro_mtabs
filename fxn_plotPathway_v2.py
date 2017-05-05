#make one py file to plot up a pathway
#modified from code used in NB project
#this version gets rid of the crazy scaling I setup ...use log10 data...(I hope)
#KLongnecker, 4 May 2017
import pandas as pd
import numpy as np
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
    
def gatherDetails(enterPathway,folderName,useCO,CO_values):
    #check if the directories exist, one for pathway files
    if not os.path.exists(folderName):
        os.makedirs(folderName)
                   
    #only one pathway at a time
    setKeep = 1
    try:
        kegg_get(enterPathway).read()
    except:
        #use the ko map if there is nothing species specific...this can also fail...
        usePathway = 'ko' + enterPathway[3:8]
        setKeep = 0
        try:
            kegg_get(usePathway).read()
        except:
            pass
        
    if setKeep:
        usePathway = enterPathway

    #get the compounds and genes for this pathway     
    genes = getKfrom_ko(usePathway)
    compounds = getCfrom_ko(usePathway)
    
    #figure out which ones I have data for...
    setC = set(compounds)
    setT = set(useCO)
    intCompounds = setC.intersection(setT)
        
    ## plot the pathway map for this pathway, get details from KEGG for plotting (%must be at least 4 colors)
    useColors = pal.colorbrewer.diverging.PuOr_4.hex_colors
    #useColors = pal.colorbrewer.diverging.RdYlBu_11.hex_colors
    
    #set the color of the mtab based on its value, only scale the values from this particular pathway
    useCOsubset = CO_values.loc[intCompounds]
    cmin = useCOsubset[useCOsubset>0].min() #only consider values greater than 0
    cmax = useCOsubset.replace([np.inf],np.nan).dropna(how = 'all').max()
    
    size = 20 #increase the size of the compounds in the plots
        
    #can have all zeros...
    if sum(useCOsubset.dropna())==0:
        pass
        #print('No measured metabolites in pathway ' + usePathway)len(
    elif len(useCOsubset.value_counts())==1:
        #only two color options: yes/no
        dummy = useCOsubset.copy(deep = True)
        dummy.replace([np.inf],np.nan,inplace = True)
        for idx,item in enumerate(useCOsubset):
            if np.isnan(item):
                useCOsubset.iloc[idx] = int(0)
            else:
                useCOsubset.iloc[idx] = int(1) 
        
        #go get the pathway information and customize the plot
        pathway = KGML_parser.read(kegg_get(usePathway, "kgml")) #no choice in gene color: green

        # Change the colors of compounds
        for element in pathway.compounds:
            for graphic in element.graphics:
                tc = element.name[4:10] #skip over the 'cpd:'
                if (tc in intCompounds):
                    #in the pathway, set the color
                    tempColor = useCOsubset.loc[tc]
                    graphic.bgcolor = useColors[int(tempColor)] 
                    graphic.width = size
                    graphic.height = size

        canvas = KGMLCanvas(pathway, import_imagemap=True)
        pdfName = 'map' + str(usePathway[3:8]) + '_strain_' + usePathway[0:3] + '.pdf'
        canvas.draw(folderName + '/' + pdfName)
        pdfName = None #empty it in case that is where I am having issues         
        
    else:
        #dummy = useCOsubset.copy(deep = True)
        #dummy.replace([np.inf],np.nan,inplace = True)
        #for idx,item in enumerate(useCOsubset):
        #   if np.isnan(item):
        #      useCOsubset.iloc[idx] = 0
        #    elif np.isinf(item):
        #        useCOsubset.iloc[idx] = 10*cmax #make inf 10x the biggest value

        #now, find cmax again...use that downstream
        #cmax = useCOsubset.replace([np.inf],np.nan).dropna(how = 'all').max()

        #use histogram to make the bins (complete hack)
        #Tracer()()
        #now...put zero at beginning and inf at end
        #BUT - can actually have cases with values for all metabolites (novel concept)

        try:
            nz = useCOsubset.value_counts()[0] #count the number of zeros
            useColors.insert(0,'#31a354') #insert green at beginning only if there are zeros
            if np.array_equal(cmin,cmax):
                a,bin_edges = np.histogram(useCOsubset,bins = len(useColors)-2,range = (0,cmax))
            else:
                a,bin_edges = np.histogram(useCOsubset,bins = len(useColors)-2,range = (cmin,cmax))
                a = np.insert(a,0,nz)
                bin_edges = np.insert(bin_edges,0,0)
        except: #no zeros in the pathway (happens!)
            a,bin_edges = np.histogram(useCOsubset,bins = len(useColors)-2,range = (cmin,cmax))

        #then find the index for each number...this will be the index into useColors
        #Tracer()()
        useIdx = np.digitize(useCOsubset,bin_edges)-1
            
        color_df = pd.DataFrame({'mtab': useCOsubset,'idx':useIdx})
        #Tracer()()
        #go get the pathway information and customize the plot
        pathway = KGML_parser.read(kegg_get(usePathway, "kgml")) #no choice in gene color: green

        # Change the colors of compounds
        for element in pathway.compounds:
            for graphic in element.graphics:
                tc = element.name[4:10] #skip over the 'cpd:'
                if (tc in intCompounds):
                    #in the pathway, set the color
                    tempColor = color_df.loc[tc,'idx']
                    graphic.bgcolor = useColors[int(tempColor)] 
                    graphic.width = size
                    graphic.height = size

        canvas = KGMLCanvas(pathway, import_imagemap=True)
        pdfName = 'map' + str(usePathway[3:8]) + '_strain_' + usePathway[0:3] + '.pdf'
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
