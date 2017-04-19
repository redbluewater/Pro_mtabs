#make one py file to plot up a figure...
#KLongnecker, 18 April 2017
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import palettable as pal
from Bio import SeqIO
from Bio.KEGG.REST import *
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
from IPython.display import Image, HTML
import os

from IPython.core.debugger import Tracer
#os._exit(0)

def plotGroup(oneGroup,prunedBRITE,useCO,mtabPruned,oneStrain):
    #oneStrain = 'pmg'
    shortList = prunedBRITE.loc[(prunedBRITE['B']==oneGroup)] 
    onePath = shortList.loc[:,'map']
    onePath_ann = []
    for item in onePath:
        onePath_ann.append(oneStrain + item)

    gatherGroup = pd.DataFrame()
    for one in onePath_ann:
        #Tracer()()
        mCpds = set(getCfrom_ko(one))
        ProData= set(useCO)
        handh = mCpds.intersection(ProData)
        for cpd in handh:
            #print(cpd)
            tm = mtabPruned.loc[cpd,:]
            gatherGroup = gatherGroup.append(tm)

    hfont = {'fontname':'Palatino'}
    plt.title(oneGroup)
    plt.xlabel('xlabel', **hfont)
    plt.pcolor(gatherGroup,cmap = 'PRGn')
    plt.yticks(np.arange(0.5, len(gatherGroup.index), 1), gatherGroup.index,fontsize = 5)
    #plt.xticks(np.arange(0.5, len(gatherGroup.columns), 1), gatherGroup.columns,rotation = 45)
    plt.xticks(np.arange(0.5,(len(list(mtabPruned)) + 0.1),1), gatherGroup.columns,rotation = 90)

    plt.show()
    
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