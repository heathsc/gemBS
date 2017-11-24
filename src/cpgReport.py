# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 12:02:34 2017

@author: marcos
"""

import os
import json

import matplotlib.pyplot as plt
import matplotlib.pylab as pylab

import numpy as np

from variantsReport import HtmlVariantsReport
from reportStats import RunBasicStats
from report import BasicHtml
from sphinx import BasicSphinx
from sphinx import ConfigSphinx

class CpgStats(object):
    """ Basic definition of cpg statistics """
    
    def __init__(self):
        """ Initializates a dictionary of basic values """

        self.data = {}
        self.methylation = {}
        
    def update(self,dictionary,methDict,infoReadsDict):
        """ Adds to previous data new set of values
        
            dictionary - Dictionary of values to be added
        """        
        self.data = dictionary
        self.methylation = methDict
        self.information_reads = infoReadsDict

        self.vTotalInfoReadsCpgs = []
        self.vReferenceCpGsInfoReadsCpgs = []
        self.vNonReferenceCpGsInfoReadsCpgs = []

        self.vTotalInfoReadsCpgsQ20 = []
        self.vReferenceCpGsInfoReadsCpgsQ20 = []
        self.vNonReferenceCpGsInfoReadsCpgsQ20 = []
        
        self.vAxisInfoReads = []
        
        
                        
    def getPercentage(self,total,subtotal):
        """Calculates the percentage for a given value"""
        if total != 0.0:
            return (float(subtotal)/float(total))*100
        else:
            return 0.0
            
    def getTableCpGs(self):
        """Returns a vector of CpGs Values"""
        vector = []
        
        nTotal = self.data["TotalDinucleotides"]
        nTotalQ20 = self.data["TotalQuality20"]
        nReferenceCG = self.data["TotalReferenceCGs"]
        nReferenceCGQ20 = self.data["TotalReferenceCGsQuality20"]
        nNonReferenceCG = self.data["TotalNonReferenceCGs"]
        nNonReferenceCGQ20 = self.data["TotalNonReferenceCGsQuality20"]
        
        vector.append(["Dinucleotide CG*","%i" %(nTotal),"%i" %(nTotalQ20),"%.2f %%" %(self.getPercentage(nTotal,nTotalQ20))])        
        vector.append(["Reference CG","%i" %(nReferenceCG),"%i" %(nReferenceCGQ20),"%.2f %%" %(self.getPercentage(nReferenceCG,nReferenceCGQ20))])
        vector.append(["Non Reference CG","%i" %(nNonReferenceCG),"%i" %(nNonReferenceCGQ20),"%.2f %%" %(self.getPercentage(nNonReferenceCG,nNonReferenceCGQ20))])
                        
        return vector

    def getTableSummary(self):
        """ Returns a vector of Summary Values"""
        vector = []
        
        #1. ALL CGs (Methylated,Intermediate Methylated,UnMethylated)        
        nTotalMethylation = [self.data["TotalMethylated"],self.data["TotalMethylatedQuality20"]]
        nTotalIntermediate = [self.data["TotalIntermediateMethylated"],self.data["TotalIntermediateMethylatedQuality20"]]
        nTotalUnMethylation = [self.data["TotalUnMethylated"],self.data["TotalUnmethylatedQuality20"]]
        
        #2. REFERENCE CGs (Methylated,Intermediate Methylated,UnMethylated)
        nTotalRefCGsMethylation = [self.data["TotalReferenceCGsMethylated"],self.data["TotalReferenceCGsMethylatedQuality20"]]
        nTotalRefCGsIntermediate = [self.data["TotalReferenceCGsIntermediateMethylated"],self.data["TotalReferenceCGsIntermediateMethylatedQuality20"]]
        nTotalRefCGsUnMethylation = [self.data["TotalReferenceCGsUnMethylated"],self.data["TotalReferenceCGsUnMethylatedQuality20"]]
        
        #3. NON REFERENCE CGs (Methylated,Intermediate Methylated,UnMethylated)
        nTotalNonRefCGsMethylation = [self.data["TotalNonReferenceCGsMethylated"],self.data["TotalNonReferenceCGsMethylatedQuality20"]]
        nTotalNonRefCGsIntermediate = [self.data["TotalNonReferenceCGsIntermediateMethylated"],self.data["TotalNonReferenceCGsIntermediateMethylatedQuality20"]]
        nTotalNonRefCGsUnMethylation = [self.data["TotalNonReferenceCGsUnMethylated"],self.data["TotalNonReferenceCGsUnMethylatedQuality20"]]
    
        nTotal = self.data["TotalDinucleotides"]
        nTotalQ20 = self.data["TotalQuality20"]
        vector.append(["All CGs","Methylated","%i" %(nTotalMethylation[0]),"%.2f %%" %(self.getPercentage(nTotal,nTotalMethylation[0])),"%i"%(nTotalMethylation[1]),"%.2f %%" %(self.getPercentage(nTotalQ20,nTotalMethylation[1]))])
        vector.append(["All CGs","Intermediate Methylated","%i" %(nTotalIntermediate[0]),"%.2f %%" %(self.getPercentage(nTotal,nTotalIntermediate[0])),"%i"%(nTotalIntermediate[1]),"%.2f %%" %(self.getPercentage(nTotalQ20,nTotalIntermediate[1]))])
        vector.append(["All CGs","UnMethylated","%i" %(nTotalUnMethylation[0]),"%.2f %%" %(self.getPercentage(nTotal,nTotalUnMethylation[0])),"%i"%(nTotalUnMethylation[1]),"%.2f %%" %(self.getPercentage(nTotalQ20,nTotalUnMethylation[1]))])
        
        nReferenceCG = self.data["TotalReferenceCGs"]
        nReferenceCGQ20 = self.data["TotalReferenceCGsQuality20"]
        vector.append(["Reference CGs","Methylated","%i" %(nTotalRefCGsMethylation[0]),"%.2f %%" %(self.getPercentage(nReferenceCG,nTotalRefCGsMethylation[0])),"%i"%(nTotalRefCGsMethylation[1]),"%.2f %%" %(self.getPercentage(nReferenceCGQ20,nTotalRefCGsMethylation[1]))])
        vector.append(["Reference CGs","Intermediate Methylated","%i" %(nTotalRefCGsIntermediate[0]),"%.2f %%" %(self.getPercentage(nReferenceCG,nTotalRefCGsIntermediate[0])),"%i"%(nTotalRefCGsIntermediate[1]),"%.2f %%" %(self.getPercentage(nReferenceCGQ20,nTotalRefCGsIntermediate[1]))])
        vector.append(["Reference CGs","UnMethylated","%i" %(nTotalRefCGsUnMethylation[0]),"%.2f %%" %(self.getPercentage(nReferenceCG,nTotalRefCGsUnMethylation[0])),"%i"%(nTotalRefCGsUnMethylation[1]),"%.2f %%" %(self.getPercentage(nReferenceCGQ20,nTotalRefCGsUnMethylation[1]))])
        
        nNonReferenceCG = self.data["TotalNonReferenceCGs"]
        nNonReferenceCGQ20 = self.data["TotalNonReferenceCGsQuality20"]
        vector.append(["Non Reference CGs","Methylated","%i" %(nTotalNonRefCGsMethylation[0]),"%.2f %%" %(self.getPercentage(nNonReferenceCG,nTotalNonRefCGsMethylation[0])),"%i"%(nTotalNonRefCGsMethylation[1]),"%.2f %%" %(self.getPercentage(nNonReferenceCGQ20,nTotalNonRefCGsMethylation[1]))])
        vector.append(["Non Reference CGs","Intermediate Methylated","%i" %(nTotalNonRefCGsIntermediate[0]),"%.2f %%" %(self.getPercentage(nNonReferenceCG,nTotalNonRefCGsIntermediate[0])),"%i"%(nTotalNonRefCGsIntermediate[1]),"%.2f %%" %(self.getPercentage(nNonReferenceCGQ20,nTotalNonRefCGsIntermediate[1]))])
        vector.append(["Non Reference CGs","UnMethylated","%i" %(nTotalNonRefCGsUnMethylation[0]),"%.2f %%" %(self.getPercentage(nNonReferenceCG,nTotalNonRefCGsUnMethylation[0])),"%i"%(nTotalNonRefCGsUnMethylation[1]),"%.2f %%" %(self.getPercentage(nNonReferenceCGQ20,nTotalNonRefCGsUnMethylation[1]))])
        
        return vector

    def setCpgInInformationReads(self,dictionaryInfoReads={},vectorInformationReads=[],vectorCpGs=[]):
        """Sets CpG values in a vector to a given position relative to the number of information reads

           dictionaryInfoReads - Dictionary of information reads { infoReads:cpgs ,..., infoReads:cpgs }
           vectorInformationReads - Vector of information reads to be expressed in x plot axis
           vectorCpGs - Vector of CpGs, would be same size like vectorInformationReads
        """        
        #Define vector for plotting Y bars
        for infoReads in sorted(dictionaryInfoReads):
            previous_ireads = 0
            current_ireads = 0
            assigned_cpgs = 0
            for index in range(0,len(vectorInformationReads)):
                current_ireads = vectorInformationReads[index]
                
                if infoReads >= previous_ireads and infoReads <= current_ireads:
                    vectorCpGs[index] = dictionaryInfoReads[infoReads]
                    assigned_cpgs = 1
                    break
                previous_ireads = current_ireads
            if assigned_cpgs == 0:
                vectorCpGs[-1] =  dictionaryInfoReads[infoReads]

    
    def setInformationReadsValues(self):
        """ Configures Information Reads Values to plot stacked Bar plot """
        dictInformationReadsReferenceCpgs = {}
        dictInformationReadsReferenceCpgsQ20 = {}
        dictInformationReadsNonReferenceCpgs = {}
        dictInformationReadsNonReferenceCpgsQ20 = {}
        dictInformationReadsTotalCpgsTotal = {}
        dictInformationReadsTotalCpgsQ20 = {}
        
        #Create integer keys
        for infoReads,cpgs in self.information_reads["informationReadsReferenceCGs"].iteritems():
            dictInformationReadsReferenceCpgs[int(infoReads)] = cpgs
        for infoReads,cpgs in self.information_reads["informationReadsNonReferenceCGs"].iteritems():
            dictInformationReadsNonReferenceCpgs[int(infoReads)] = cpgs
        for infoReads,cpgs in self.information_reads["informationReadsReferenceCGsQ20"].iteritems():
            dictInformationReadsReferenceCpgsQ20[int(infoReads)] = cpgs
        for infoReads,cpgs in self.information_reads["informationReadsNonReferenceCGsQ20"].iteritems():
            dictInformationReadsNonReferenceCpgsQ20[int(infoReads)] = cpgs
            
        #Create Total Values
        for infoReads,cpgs in dictInformationReadsReferenceCpgs.iteritems():
            dictInformationReadsTotalCpgsTotal[infoReads] = cpgs
        for infoReads,cpgs in dictInformationReadsNonReferenceCpgs.iteritems():
            if infoReads in dictInformationReadsTotalCpgsTotal:
                dictInformationReadsTotalCpgsTotal[infoReads] += cpgs
            else:
                dictInformationReadsTotalCpgsTotal[infoReads] = cpgs
        #Create Total Values Q>20
        for infoReads,cpgs in dictInformationReadsReferenceCpgsQ20.iteritems():
            dictInformationReadsTotalCpgsQ20[infoReads] = cpgs
        for infoReads,cpgs in dictInformationReadsNonReferenceCpgsQ20.iteritems():
            if infoReads in dictInformationReadsTotalCpgsQ20:
                dictInformationReadsTotalCpgsQ20[infoReads] += cpgs
            else:
                dictInformationReadsTotalCpgsQ20[infoReads] = cpgs
                
        #Define X axis for plotting information reads and CpGs
        self.vAxisInfoReads = [1,2,3,4,5,6,7,8,9,10,15,20,30,50,75,100,150,200,250,500,1000,1250,1500,2000,3000,5000,10000]                
                
        #Initialization CpGs Vectors
        for i in range(0,len(self.vAxisInfoReads)):
            self.vTotalInfoReadsCpgs.append(0)
            self.vTotalInfoReadsCpgsQ20.append(0)
            self.vReferenceCpGsInfoReadsCpgs.append(0)
            self.vReferenceCpGsInfoReadsCpgsQ20.append(0)
            self.vNonReferenceCpGsInfoReadsCpgs.append(0)
            self.vNonReferenceCpGsInfoReadsCpgsQ20.append(0)
                
                
        #Define vector for plotting
        self.setCpgInInformationReads(dictionaryInfoReads=dictInformationReadsTotalCpgsTotal,vectorInformationReads=self.vAxisInfoReads,vectorCpGs=self.vTotalInfoReadsCpgs)
        self.setCpgInInformationReads(dictionaryInfoReads=dictInformationReadsTotalCpgsQ20,vectorInformationReads=self.vAxisInfoReads,vectorCpGs=self.vTotalInfoReadsCpgsQ20)

        self.setCpgInInformationReads(dictionaryInfoReads=dictInformationReadsReferenceCpgs,vectorInformationReads=self.vAxisInfoReads,vectorCpGs=self.vReferenceCpGsInfoReadsCpgs)
        self.setCpgInInformationReads(dictionaryInfoReads=dictInformationReadsReferenceCpgsQ20,vectorInformationReads=self.vAxisInfoReads,vectorCpGs=self.vReferenceCpGsInfoReadsCpgsQ20)
        
        self.setCpgInInformationReads(dictionaryInfoReads=dictInformationReadsNonReferenceCpgs,vectorInformationReads=self.vAxisInfoReads,vectorCpGs=self.vNonReferenceCpGsInfoReadsCpgs)
        self.setCpgInInformationReads(dictionaryInfoReads=dictInformationReadsNonReferenceCpgsQ20,vectorInformationReads=self.vAxisInfoReads,vectorCpGs=self.vNonReferenceCpGsInfoReadsCpgsQ20)        
    
    def drawStackedCoveragePlot(self,pngFile,vectorAxis=[],vectorTotalCGs=[],vectorReferenceCGs=[],vectorNonReferenceCGs=[],
                                title='Coverage CpGs'):
        """ From matplot lib plots Informative reads for three types of CpGs [Total CpGs,Reference CpGs,Non Reference CpGs]
        
            pngFile - PNG output file for Informative Reads
            vectorAxis - Vector of Axis Legends. Each value is a given number of informative reads
            vectorTotalCGs - Vector of total CGs, coordinated with vectorAxis, each value is a given number of CpGs
            vectorReferenceCGs - Vector of Reference CGs, coordinated with vectorAxis, each value is a given number of CpGs
            vectorNonReferenceCGs - Vector of Non Reference CGs, coordinated with vectorAxis, each value is a given number of CpGs 
            title - Title plot
        """  
        figure, ax = plt.subplots()
        
        axisX = np.arange(len(vectorAxis))
        width = 1

        plotTotalCGs = plt.bar(axisX,vectorTotalCGs,width,color='#2f27d6')
        plotRefCGs = plt.bar(axisX,vectorReferenceCGs,width,color='#65d627')
        plotNoRefCGs = plt.bar(axisX,vectorNonReferenceCGs,width,color='#d62728')
        
        plt.ylabel('#CpGs')
        plt.xlabel('Informative Reads')
        plt.title(title)
        plt.xticks(axisX,vectorAxis,fontsize = 8,rotation='vertical')
        
        plt.legend((plotTotalCGs[0],plotRefCGs[0],plotNoRefCGs[0]),('Total CGs*','Reference CGs','Non Reference CGs'),fontsize = 9,loc='best')
             
        pylab.savefig(pngFile)
        
        plt.close(figure)     
     
    def boxplotMethylationStatus(self,pngFile=None,vectorMethylationValues=[],title='Methylation Levels CGs'):
        """ Prints a Boxplot in a png file 
            
            pngFile - PNG file to store the data
            vectorMethylationValues - Vector of Methylation Values to be printed as a bar plot
            title - Title plot
        """
        #1. Plot Configuration
        figure, ax = plt.subplots(1,1,figsize=(5.7,5))
        
        ax.set_title(title)

        #2. Fitting Line for the histogram distribution
        #2.0 Vector of value
        x = vectorMethylationValues
        #2.3 Number of bins
        num_bins = 50
        #2.4 Data histogram
        n, bins,patches = ax.hist(x,num_bins,normed=1,facecolor='green',alpha=0.5)
        ax.set_xlim(0.0, 1.0)
        
        axes = plt.gca()
        ylim = axes.get_ylim()
                       
        #3. Axis Configuration
        ax.set_xlabel("Methylation Values")
        ax.set_ylabel("Density")
        
        #2.5 Horizontal Boxplot
        ax.boxplot(x,0,'b+',0, widths = ylim[-1]/8,positions=[ylim[-1]/2],manage_xticks=False) 
        ax.axes.set_ylim(ylim)
        

        pylab.savefig(pngFile)        
        plt.close(figure)
                

class HtmlCpgReport(HtmlVariantsReport):
    """ Builds Html CpG Report """
     
    def __init__(self):
        """ Html CpG Class Constructor """
        HtmlVariantsReport.__init__(self)
        
        
    def addRow(self,vTableHtml=None,fields=[],isOdd=True):
        """ Add Row html tags
            
            vectorHtml - Vector of tags to add
            fields - Vector of fields to be printed in a row
            isOdd - Boolean value
        """
        if isOdd == True: 
            vTableHtml.append("   <TR class=\"odd\">\n")
        else:
            vTableHtml.append("   <TR>\n")

        for field in fields:
            vTableHtml.append("     <TD> %s </TD> \n" %(field))
        
        vTableHtml.append("   </TR>\n")        


    def createStatsTable(self,color=None,cpgValues=None,tableTitle=None,headers=None):
        """Create Statistics Table
        
           color - Table color, could be green or blue  
           cpgValues - Dictionary of values to be printed
           tableTitle - Table Title
           headers - list of values to be printed as headers
           returns vector of HTML tags
        """
        vTableHtml = []   
        
        #0.Table Title
        vTableHtml.append('<H3 id="section"> %s </H3>\n' %(tableTitle))
     
        #1.Table Definition
        if color == "green":
            vTableHtml.append("  <TABLE id=\"green\">\n")
        else:
            vTableHtml.append("  <TABLE id=\"hor-zebra\">\n")

        #2.Table header
        vTableHtml.append("   <TR>\n")
        for conceptHeader in headers:        
            vTableHtml.append("    <TH scope=\"col\">%s</TH>\n"%(conceptHeader))
            
        vTableHtml.append("   </TR>\n")
        
        #3.Fullfill table
        isOdd = True
        for fields in cpgValues:    
            self.addRow(vTableHtml=vTableHtml,fields=fields,isOdd=isOdd)
            isOdd = not isOdd
        
        vTableHtml.append(" </TABLE>\n")
        return vTableHtml
    
    
    def createSummaryTable(self,color=None,cpgValues=None,tableTitle=None):
        """Create Summary Table
        
           color - Table color, could be green or blue  
           cpgValues - Vector of values to be printed
           tableTitle - Table Title
           headers - list of values to be printed as headers
           returns vector of HTML tags
        """
        vTableHtml = []   
        
        
        #Based On: https://www.w3.org/WAI/tutorials/tables/irregular/   
        #0.Table Title
        vTableHtml.append('<H3 id="section"> %s </H3>\n' %(tableTitle))
        
        #1.Table Definition
        if color == "green":
            vTableHtml.append("  <TABLE id=\"green\">\n")
        else:
            vTableHtml.append("  <TABLE id=\"hor-zebra\">\n")
            
        #2.Column group span
        vTableHtml.append("      <COL>\n")
        vTableHtml.append("          <COLGROUP span=\"2\"> </COLGROUP>\n")
        vTableHtml.append("          <COLGROUP span=\"2\"> </COLGROUP>\n")
        vTableHtml.append("      <TR>")
        vTableHtml.append("          <TD rowspan=\"2\"></TD>")
        vTableHtml.append("          <TH colspan=\"2\" scope=\"colgroup\">Total</TH>")
        vTableHtml.append("          <TH colspan=\"2\" scope=\"colgroup\">Q>20</TH>")
        vTableHtml.append("      </TR>")
        vTableHtml.append("      <TR>")
        vTableHtml.append("         <TH scope=\"col\">#</TH>")
        vTableHtml.append("         <TH scope=\"col\">%</TH>") 
        vTableHtml.append("         <TH scope=\"col\">#</TH>")
        vTableHtml.append("         <TH scope=\"col\">%</TH>")
        vTableHtml.append("      </TR>")
        
        vTableHtml.append("      <TR class=\"odd\">")
        vTableHtml.append("         <TH scope=\"row\"> %s %s </TH>" %(cpgValues[0][0],cpgValues[0][1]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[0][2]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[0][3]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[0][4]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[0][5]))
        vTableHtml.append("      </TR>")
        vTableHtml.append("      <TR class=\"odd\">")
        vTableHtml.append("         <TH scope=\"row\"> %s %s </TH>" %(cpgValues[1][0],cpgValues[1][1]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[1][2]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[1][3]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[1][4]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[1][5]))
        vTableHtml.append("      </TR>")
        vTableHtml.append("      <TR class=\"odd\">")
        vTableHtml.append("         <TH scope=\"row\"> %s %s </TH>" %(cpgValues[2][0],cpgValues[2][1]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[2][2]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[2][3]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[2][4]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[2][5]))
        vTableHtml.append("      </TR>")
        
        vTableHtml.append("      <TR>")
        vTableHtml.append("         <TH scope=\"row\"> %s %s </TH>" %(cpgValues[3][0],cpgValues[3][1]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[3][2]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[3][3]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[3][4]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[3][5]))
        vTableHtml.append("      </TR>")
        vTableHtml.append("      <TR>")
        vTableHtml.append("         <TH scope=\"row\"> %s %s </TH>" %(cpgValues[4][0],cpgValues[4][1]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[4][2]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[4][3]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[4][4]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[4][5]))
        vTableHtml.append("      </TR>")
        vTableHtml.append("      <TR>")
        vTableHtml.append("         <TH scope=\"row\"> %s %s </TH>" %(cpgValues[5][0],cpgValues[5][1]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[5][2]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[5][3]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[5][4]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[5][5]))
        vTableHtml.append("      </TR>")
               
        vTableHtml.append("      <TR class=\"odd\">")
        vTableHtml.append("         <TH scope=\"row\"> %s %s </TH>" %(cpgValues[6][0],cpgValues[6][1]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[6][2]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[6][3]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[6][4]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[6][5]))
        vTableHtml.append("      </TR>")
        vTableHtml.append("      <TR class=\"odd\">")
        vTableHtml.append("         <TH scope=\"row\"> %s %s </TH>" %(cpgValues[7][0],cpgValues[7][1]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[7][2]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[7][3]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[7][4]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[7][5]))
        vTableHtml.append("      </TR>")
        vTableHtml.append("      <TR class=\"odd\">")
        vTableHtml.append("         <TH scope=\"row\"> %s %s </TH>" %(cpgValues[8][0],cpgValues[8][1]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[8][2]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[8][3]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[8][4]))
        vTableHtml.append("         <TD> %s </TD>" %(cpgValues[8][5]))
        vTableHtml.append("      </TR>")
        
        vTableHtml.append(" </TABLE>\n")
               
        return vTableHtml
        
        
    def run(self,vectorHtml=None,cpgStats=None,output_dir=None,name=None,parentName=None,parentDocument=None):
        """ Run Html CpG Report 
        
            vectorHtml - Vector of html tags
            cpgStats - CpG Stats instance 
            ouput_dir - Ouput Directory to store the report
            name - Basic name
            parentName - Parent Name
            parentDocument - Parent Document HTML
        """    
        
        #1. Add HTML Report Header
        self.addHtmlStartReport(vectorHtml=vectorHtml,parentName=parentName,currentName=name,parentDocument=parentDocument)
        
        #2. Add Total Dinucleotides Tables
        self.addSpaceSection(vectorHtml=vectorHtml) 
        vectorHtml.extend(self.createStatsTable(color='blue',cpgValues=cpgStats.getTableCpGs(),tableTitle="Dinucleotides Statistics",headers=["Type","Total","Q>20","%"]))
        
        #3.Informative Reads Plots
        self.addSpaceSection(vectorHtml=vectorHtml) 
        infoReadsPlot = "%s/IMG/%s.information.reads.png" %(output_dir,name)
        infoReadsQ20Plot = "%s/IMG/%s.information.reads.q20.png" %(output_dir,name)
        cpgStats.setInformationReadsValues()
        
        cpgStats.drawStackedCoveragePlot(infoReadsPlot,vectorAxis=cpgStats.vAxisInfoReads,vectorTotalCGs=cpgStats.vTotalInfoReadsCpgs,
                                         vectorReferenceCGs=cpgStats.vReferenceCpGsInfoReadsCpgs,vectorNonReferenceCGs=cpgStats.vNonReferenceCpGsInfoReadsCpgs,title='Coverage CpGs')
                                    
        cpgStats.drawStackedCoveragePlot(infoReadsQ20Plot,vectorAxis=cpgStats.vAxisInfoReads,vectorTotalCGs=cpgStats.vTotalInfoReadsCpgsQ20,
                                         vectorReferenceCGs=cpgStats.vReferenceCpGsInfoReadsCpgsQ20,vectorNonReferenceCGs=cpgStats.vNonReferenceCpGsInfoReadsCpgsQ20,title='Coverage CpGs Q>20')
        
        vectorHtml.extend(self.buildTwoPlots(color="green",tableTitleFirst="All CpGs Informative Reads",tableTitleSecond="Q>20 CpGs Informative Reads",
                    pathFirstImage="IMG/%s"%(os.path.basename(infoReadsPlot)),pathSecondImage="IMG/%s"%(os.path.basename(infoReadsQ20Plot)))) 
        
        #4. Plot for Methylation Level Status
        self.addSpaceSection(vectorHtml=vectorHtml)
        methylationLevelAllCGsPlot = "%s/IMG/%s.all.cgs.png" %(output_dir,name)
        methylationLevelReferenceCGsPlot = "%s/IMG/%s.reference.cgs.png" %(output_dir,name)
        methylationLevelNonReferenceCGsPlot = "%s/IMG/%s.nonreference.cgs.png" %(output_dir,name)
        
        totalList = cpgStats.methylation["referenceCGsMethValuesQ20"]
        totalList.extend(cpgStats.methylation["nonReferenceCGsMethValuesQ20"])
        
        cpgStats.boxplotMethylationStatus(pngFile=methylationLevelAllCGsPlot,vectorMethylationValues=totalList,title='Methylation Levels All CGs Q>20')
        cpgStats.boxplotMethylationStatus(pngFile=methylationLevelReferenceCGsPlot,vectorMethylationValues=cpgStats.methylation["referenceCGsMethValuesQ20"],title='Methylation Levels Reference CGs Q>20')
        cpgStats.boxplotMethylationStatus(pngFile=methylationLevelNonReferenceCGsPlot,vectorMethylationValues=cpgStats.methylation["nonReferenceCGsMethValuesQ20"],title='Methylation Levels Non Reference CGs Q>20')
        
        vectorHtml.extend(self.buildThreePlots(color="green",tableTitleFirst="Methylation levels all CGs Q>20",
                              tableTitleSecond="Methylation levels Reference CGs Q>20",tableTitleThird="Methylation levels Non Reference CGs Q>20",
                              pathFirstImage="IMG/%s"%(os.path.basename(methylationLevelAllCGsPlot)),
                              pathSecondImage="IMG/%s"%(os.path.basename(methylationLevelReferenceCGsPlot)),
                              pathThirdImage="IMG/%s"%(os.path.basename(methylationLevelNonReferenceCGsPlot))))
                              
        #5.Final Summary Table                               
        self.addSpaceSection(vectorHtml=vectorHtml) 
        vectorHtml.extend(self.createSummaryTable(color='blue',cpgValues=cpgStats.getTableSummary(),tableTitle="Summary"))
       

        

class IndexCpgHtml(HtmlCpgReport):
    """ Class which defines Index Sample Report """
    
    def __init__(self,output_dir=None,name_project=None,dictionary_samples=None):
        """  Class constructor
        
             output_dir -- Output directory to store HTML reports
             name_project -- Project name
             dictionary_samples -- Dictionary of samples and list of stats objects. {"sample":"cpgStats"}
        """
        #Call Parent class constructor
        HtmlVariantsReport.__init__(self)        
        
        self.output_dir = output_dir
        self.name_project = name_project
        self.cpg_html_document = "%s/%s.html" %(self.output_dir,self.name_project)
        self.dictionary_samples = dictionary_samples
        
    
    def createLinksTables(self,vector_html=None,vector_links=None,title=None,color=None):
        """ Create a set of link tables 
        
            vector_html - vector of html tags
            vector_links - vector of html documents for build acces table [[htmlN,sample_name_N] , [htmlN+1,sample_name_N+1] ]
            title - Table header
            color - color envioronment
        """
        
        if color == "green":
            vector_html.append("  <TABLE id=\"green\">\n")
        else:
            vector_html.append("  <TABLE id=\"hor-zebra\">\n")
            
        vector_html.append("   <TR> <TH scope=\"col\"> %s </TH> </TR> \n" %(title))

        linksClass = ''	
        if color == "green":
            linksClass += '"linkGreen"'
        else:
            linksClass += '"link"'
  
        isOdd = True
        for sample_link in vector_links:
            if isOdd == True: 
                vector_html.append('   <TR class="odd"> <TD> <a class=%s href="%s"> %s </TD> </TR> \n' %(linksClass,sample_link[0],sample_link[1]))
            else:
                vector_html.append('   <TR> <TD> <a class=%s href="%s"> %s </TD> </TR> \n' %(linksClass,sample_link[0],sample_link[1]))

            isOdd = not isOdd
            
        vector_html.append(" </TABLE>\n")
        
        
    def run(self,vectorHtml=None):
        """ Run Lane HTML Documentation Building
        
            vectorHtml - vector of html tags
        """ 
        #1. Add HTML Report Header
        self.addHtmlReportHeader(vectorHtml)
        
        vectorHtml.append('<H1 id="title"> <U> CpG Report from Methylation Pipeline Project %s </U> </H1>' %(self.name_project))
        
        vector_sample_links = []
        
        for sample, stats in sorted(self.dictionary_samples.iteritems()):            
            sampleHtml = "%s/%s_cpg.html" %(self.output_dir,sample)
            #Create Build HTML
            vSampleHtml = []    
            
            sampleCpgReport = HtmlCpgReport()
            sampleCpgReport.run(vectorHtml=vSampleHtml,cpgStats=stats,output_dir=self.output_dir,name=sample,parentName=self.name_project,parentDocument=self.cpg_html_document)    
            self.closeHtmlReport(vectorHtml=vSampleHtml)
            RunBasicStats.saveDocument(file_name=sampleHtml,vectorContent=vSampleHtml)            
            
            vector_sample_links.append([os.path.basename(sampleHtml),sample])
            
            
        #2.Links table
        self.createLinksTables(vector_html=vectorHtml,vector_links=vector_sample_links,title="Sample Reports",color="blue")
        
        #3.Close Html
        self.closeHtmlReport(vectorHtml=vectorHtml)


class SphinxCpgsReport(BasicSphinx):
    """ Class responsable for printing Sphinx report per CpG """
    
    def __init__(self):
        """Sphinx Variants Report building"""
        
        #Call Parent class constructor
        BasicSphinx.__init__(self,mapping_stats=None,png_mapq_histogram="",png_insert_size_histogram="")   
        
        
    def addRow(self,ident=20,lenCell=30,vectorSphinx=None,fields=None):
        """ Add row sphinx
            
            ident - Characters to ident
            lenCell - Length of the Cell in each row
            vectorSphinx - Vector of tags to add
            fields - Vector of fields to be printed in a row
        """
        cells = 0
        cellContents = ""
        #Add ident chars
        for i in range(ident):
            cellContents += " "        
        
        lenFields = len(fields)

        for field in fields:
            isLeft = False
            isRight = False
            if cells == 0:
                isLeft = True
            elif cells == (lenFields-1):
                isRight = True

            cellContents += self.addCell(text="%s" %(field),lenCell=lenCell,isLeft=isLeft,isRight=isRight)
            cells += 1
        
        vectorSphinx.append(cellContents) 
            
        #Underline
        self.addSimpleLine(ident=ident,vectorSphinx=vectorSphinx,cells=cells,lenCell=lenCell)        
        

    def createStatsTable(self,ident=20,lenCell=30,vectorValues=None,headers=None):
        """Create Statistics Table
        
           ident - Characters to ident 
           lenCell - Cell length
           vectorValues - Vector of values
           headers - List Of Headers
           returns vector of Sphinx tags
        """
        
        vTableSphinx = []   
        cellsHeader = 0
        
        #1.Table header
        cellsHeader = len(headers)
        #2. Write overline
        self.addSimpleLine(ident=ident,vectorSphinx=vTableSphinx,cells=cellsHeader,lenCell=lenCell)
  
        #3.Stats table
        cellContents = ""
        
        #Add ident chars
        for i in range(ident):
            cellContents += " "  
       
        i = 0
        total = len(headers)
        
        for conceptHeader in headers:
            isLeft = False
            isRight = False
            if i == 0:
                isLeft = True
            elif i == total-1:
                isRight = True
                
            cellContents += self.addCell(text=conceptHeader,lenCell=lenCell,isLeft=isLeft,isRight=isRight)
            i = i + 1
            
        vTableSphinx.append(cellContents)
        
        #Header Line 
        self.addHeaderLine(ident=ident,vectorSphinx=vTableSphinx,cells=cellsHeader,lenCell=lenCell)
                
        #2. Table Contents
        for rowValue in vectorValues:       
            self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,fields=rowValue) 
 
            
        return vTableSphinx
        
        
    def buildImage(self,ident=20,pathImage=None):
        """ Build image
  
            pathImage - PNG file to be shown at sphinx document      
        """
        #1. SPHINX TAG BUILDING 
        vectorSphinx = []
   
        imageLine = "    .. image:: %s\n" %(pathImage)
        vectorSphinx.append(imageLine)

        return vectorSphinx
        
        
    def run(self,vectorSphinx=None,ident=20,lenCell=30,cpgStats=None,output_dir=None,name=None):
        """ Run Lane SPHINX Documentation Building
        
            vectorSphinx - vector of Sphinx tags
            ident - Characters to ident
            lenCell - Length of table cells in characters
            cpgStats - CpG Stats Instance
            ouput_dir - Ouput Directory to store the report
            name - Basic name
        """  
        
        #0. Add Chapter Line       
        self.addTopSection(ident=ident,vectorSphinx=vectorSphinx,title="CpGs Statistics for %s" %(name))   
        vectorSphinx.append("\n")
        
        #1. Add Total Dinucleotides Tables
        self.addSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Dinucleotides Statistics")
        vectorSphinx.append("\n")
        vectorSphinx.extend(self.createStatsTable(ident=ident,lenCell=lenCell,vectorValues=cpgStats.getTableCpGs(),headers=["Type","Total","Q>20","%"]))
        vectorSphinx.append("\n")  
          
        #3. Information Reads Stacked Plot
        vectorSphinx.append("\n")
        cpgStats.setInformationReadsValues()
        infoReadsPlot = "%s/%s_information_reads.png" %(output_dir,name)
        infoReadsQ20Plot = "%s/%s_information_reads_q20.png" %(output_dir,name)
        cpgStats.drawStackedCoveragePlot(infoReadsPlot,vectorAxis=cpgStats.vAxisInfoReads,vectorTotalCGs=cpgStats.vTotalInfoReadsCpgs,
                                         vectorReferenceCGs=cpgStats.vReferenceCpGsInfoReadsCpgs,vectorNonReferenceCGs=cpgStats.vNonReferenceCpGsInfoReadsCpgs,title='Coverage CpGs')
                                    
        cpgStats.drawStackedCoveragePlot(infoReadsQ20Plot,vectorAxis=cpgStats.vAxisInfoReads,vectorTotalCGs=cpgStats.vTotalInfoReadsCpgsQ20,
                                         vectorReferenceCGs=cpgStats.vReferenceCpGsInfoReadsCpgsQ20,vectorNonReferenceCGs=cpgStats.vNonReferenceCpGsInfoReadsCpgsQ20,title='Coverage CpGs Q>20')
        
        self.addSubSubSection(ident=ident,vectorSphinx=vectorSphinx,title="All CpGs Information Reads")
        vectorSphinx.append("\n")
        vectorSphinx.extend(self.buildImage(ident=0,pathImage=os.path.basename(infoReadsPlot)))
        vectorSphinx.append("\n")
        self.addSubSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Q>20 CpGs Information Reads")
        vectorSphinx.append("\n")
        vectorSphinx.extend(self.buildImage(ident=0,pathImage=os.path.basename(infoReadsQ20Plot)))
        vectorSphinx.append("\n")
        
        #4. Plot for Methylation Level Status
        vectorSphinx.append("\n")
        methylationLevelAllCGsPlot = "%s/%s_all_cgs.png" %(output_dir,name)
        methylationLevelReferenceCGsPlot = "%s/%s_reference_cgs.png" %(output_dir,name)
        methylationLevelNonReferenceCGsPlot = "%s/%s_nonreference_cgs.png" %(output_dir,name)
        
        totalList = cpgStats.methylation["referenceCGsMethValuesQ20"]
        totalList.extend(cpgStats.methylation["nonReferenceCGsMethValuesQ20"])
        
        cpgStats.boxplotMethylationStatus(pngFile=methylationLevelAllCGsPlot,vectorMethylationValues=totalList,title='Methylation Levels All CGs Q>20')
        cpgStats.boxplotMethylationStatus(pngFile=methylationLevelReferenceCGsPlot,vectorMethylationValues=cpgStats.methylation["referenceCGsMethValuesQ20"],title='Methylation Levels Reference CGs Q>20')
        cpgStats.boxplotMethylationStatus(pngFile=methylationLevelNonReferenceCGsPlot,vectorMethylationValues=cpgStats.methylation["nonReferenceCGsMethValuesQ20"],title='Methylation Levels Non Reference CGs Q>20')
        
        self.addSubSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Methylation levels all CGs Q>20")
        vectorSphinx.append("\n")
        vectorSphinx.extend(self.buildImage(ident=0,pathImage=os.path.basename(methylationLevelAllCGsPlot)))
        vectorSphinx.append("\n")
        self.addSubSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Methylation levels Reference CGs Q>20")
        vectorSphinx.append("\n")
        vectorSphinx.extend(self.buildImage(ident=0,pathImage=os.path.basename(methylationLevelReferenceCGsPlot)))
        vectorSphinx.append("\n")
        self.addSubSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Methylation levels Non Reference CGs Q>20")
        vectorSphinx.append("\n")
        vectorSphinx.extend(self.buildImage(ident=0,pathImage=os.path.basename(methylationLevelNonReferenceCGsPlot)))
        vectorSphinx.append("\n")
        
        #5. Summary Table
        self.addSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Summary")
        vectorSphinx.append("\n")
        vectorSphinx.extend(self.createStatsTable(ident=ident,lenCell=lenCell,vectorValues=cpgStats.getTableSummary(),headers=["Type","Status","Total #","Total %","Q>20 #","Q>20 %"]))
        vectorSphinx.append("\n")


class SumupCpgSphinx(BasicSphinx):
    """ Class which defines Index Sample Report """

    def __init__(self,output_dir=None,name_project=None,dictionary_samples=None):
        """  Class constructor
        
             output_dir -- Output directory to store HTML reports
             name_project -- Project name
             vector_sample -- Vector of samples
        """
        #Call Parent class constructor
        BasicSphinx.__init__(self)
        
        self.output_dir = output_dir
        self.name_project = name_project
        self.cpg_sphinx_document = "%s/%s.rst" %(self.output_dir,self.name_project)
        self.dictionary_samples = dictionary_samples
                
    def run(self,vectorSphinx=None):
        """ Run Sumup Sphinx Documentation Building
        
            vectorSphinx - vector of html tags
        """ 
        #1. Add SPHINX Report Header
        #Scaped name
        scaped_name = ""
        for character in self.name_project:
            if character == "_":
                scaped_name += "\\%s" %(character)
            else:
                scaped_name += "%s" %(character) 
                
        self.addPartLine(ident=0,vectorSphinx=vectorSphinx,title="CpG Report from Methylation Pipeline Project %s"%(scaped_name))

        vectorDocuments = []
          
        #2. Create Sample and lanes Report 
        for sample, stats in sorted(self.dictionary_samples.iteritems()):            
            sampleSphinx = "%s/%s.rst" %(self.output_dir,sample)
            vectorDocuments.append("%s" %(sample))
       
            #Create Build SPHINX
            vector_sphinx = []
            sphinxReport = SphinxCpgsReport()
            
            sphinxReport.run(vectorSphinx=vector_sphinx,ident=0,lenCell=30,cpgStats=stats,output_dir=self.output_dir,name=sample)

            #Store Sphinx file
            RunBasicStats.saveDocument(file_name=sampleSphinx,vectorContent=vector_sphinx)
             
        #3.Toc Table
        vectorSphinx.append("\n")
        self.addTocTree(vectorSphinx=vectorSphinx,vectorDocuments=vectorDocuments)
               

def buildCpgReport(inputs=None,output_dir=None,name=None):
    """ Build CpG report.
    
        inputs -- dictionary of samples and json stats {"sample":"stats_cpg.json,meth_values_cpg.json"}
        output_dir -- Output directory to store html documents.
        name --  Name basic to build output results.
    """
    
    #Check output directory
    if not os.path.exists("%s/IMG/" %(output_dir)):
        os.makedirs("%s/IMG/" %(output_dir))   
        
    #Dictionary Samples CpgStats {"sample":"stats_cpg.json"}
    dict_samples_cpg = {} 

    for sample,json_files in inputs.iteritems(): 
        #CpgStats object
        cpgStats = CpgStats()
        #Json Parsing
        with open(json_files[0], 'r') as file_json:
            with open(json_files[1], 'r') as meth_json:
                with open(json_files[2], 'r') as information_reads_json:
                    cpgStats.update(dictionary=json.load(file_json),methDict=json.load(meth_json),infoReadsDict=json.load(information_reads_json))
            
        #Update directory of samples
        dict_samples_cpg[sample] = cpgStats
        
    #IndexHtml object
    vector_index_html = []
    indexCpgHtml = IndexCpgHtml(output_dir=output_dir,name_project=name,dictionary_samples=dict_samples_cpg)
    indexCpgHtml.run(vectorHtml=vector_index_html)
    RunBasicStats.saveDocument(file_name=indexCpgHtml.cpg_html_document,vectorContent=vector_index_html)
    
    #CSS object
    cssBuilder = BasicHtml()
    vector_css = []
    cssBuilder.buildStyleSheet(vector_css)
    RunBasicStats.saveDocument(file_name="%s/style.css" %(output_dir),vectorContent=vector_css) 
        


def buildSphinxCpgReport(inputs=None,output_dir=None,name=None):
    """ Build Sphinx Variant report
    
        inputs -- dictionary of samples and json stats {"sample":"stats_cpg.json,meth_values_cpg.json,info_reads_cpg.json"}
        output_dir -- Output directory to store html documents.
        name --  Name basic to build output results.
    """
        
    #Check output directory
    if not os.path.exists("%s/IMG/" %(output_dir)):
        os.makedirs("%s/IMG/" %(output_dir)) 
        
    #Dictionary Samples CpgStats {"sample":"stats_cpg.json"}
    dict_samples_cpg = {} 

    for sample,json_files in inputs.iteritems(): 
        #CpgStats object
        cpgStats = CpgStats()
        #Json Parsing
        with open(json_files[0], 'r') as file_json:
            with open(json_files[1], 'r') as meth_json:
                with open(json_files[2], 'r') as information_reads_json:
                    cpgStats.update(dictionary=json.load(file_json),methDict=json.load(meth_json),infoReadsDict=json.load(information_reads_json)) 
            
        #Update directory of samples
        dict_samples_cpg[sample] = cpgStats        
        
    #SumupSphinx object
    vector_sumup_sphinx = []
    sumupCpgSphinx = SumupCpgSphinx(output_dir=output_dir,name_project=name,dictionary_samples=dict_samples_cpg)
    sumupCpgSphinx.run(vectorSphinx=vector_sumup_sphinx)
    RunBasicStats.saveDocument(file_name=sumupCpgSphinx.cpg_sphinx_document,vectorContent=vector_sumup_sphinx)
    
    #Config python file
    cfgFile = ConfigSphinx(path_config_file="%s/conf.py" %(output_dir),path_makefile_file="%s/Makefile" %(output_dir),master_file=name,project_name=name,main_title='CPGs REPORT')
    cfgFile.run()     
        
        