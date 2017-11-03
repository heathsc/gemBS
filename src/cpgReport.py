# -*- coding: utf-8 -*-
"""
Created on Thu Apr  6 12:02:34 2017

@author: marcos
"""

import os
import json
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import matplotlib.mlab as mlab
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
        
    def update(self,dictionary,methDict):
        """ Adds to previous data new set of values
        
            dictionary - Dictionary of values to be added
        """        
        self.data = dictionary
        self.methylation = methDict
                        
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
        
        vector.append(["Dinucleotides",self.data["TotalDinucleotides"],"100 %"])
        vector.append(["Q>20","%i" %(nTotalQ20),"%.2f %%" %(self.getPercentage(nTotal,nTotalQ20))])
                
        return vector

    def getTableZygosity(self):
        """ Returns a vector of Zygosity Values"""
        vector = []
        
        nTotal = self.data["TotalDinucleotides"]
        nTotalHomozygous = self.data["TotalHomozygous"]
        nTotalHomozygousQ20 = self.data["TotalHomozygousQuality20"]
        nTotalAlternativeCX = self.data["TotalAlternativeCX"]
        nTotalAlternativeCXQ20 = self.data["TotalAlternativeCXQuality20"]
        
        vector.append(["Homozygous","%i"%(nTotalHomozygous),"%.2f %%" %(self.getPercentage(nTotal,nTotalHomozygous))])
        vector.append(["Homozygous Q>20","%i"%(nTotalHomozygousQ20),"%.2f %%" %(self.getPercentage(nTotal,nTotalHomozygousQ20))])
        vector.append(["AlternativeCX","%i"%(nTotalAlternativeCX),"%.2f %%" %(self.getPercentage(nTotal,nTotalAlternativeCX))])
        vector.append(["AlternativeCX Q>20","%i"%(nTotalAlternativeCXQ20),"%.2f %%" %(self.getPercentage(nTotal,nTotalAlternativeCXQ20))])
        
        return vector


    def getTableMethylationStatus(self):
        """ Returns a vector of Methylation Value"""
        vector = []
        
        nTotal = self.data["TotalDinucleotides"]
        nTotalMethylated = self.data["TotalMethylated"]
        nTotalMethylatedQ20 = self.data["TotalMethylatedQuality20"]
        nTotalIntermediateMethylated = self.data["TotalIntermediateMethylated"]
        nTotalIntermediateMethylatedQ20 = self.data["TotalIntermediateMethylatedQuality20"]
        nTotalUnmethylated = self.data["TotalUnMethylated"]
        nTotalUnmethylatedQ20 = self.data["TotalUnmethylatedQuality20"]

        vector.append(["Methylated","%i"%(nTotalMethylated),"%.2f %%" %(self.getPercentage(nTotal,nTotalMethylated))])
        vector.append(["Methylated Q>20","%i"%(nTotalMethylatedQ20),"%.2f %%" %(self.getPercentage(nTotal,nTotalMethylatedQ20))])
        vector.append(["Intermediate Methylated","%i"%(nTotalIntermediateMethylated),"%.2f %%" %(self.getPercentage(nTotal,nTotalIntermediateMethylated))])
        vector.append(["Intermediate Methylated Q>20","%i"%(nTotalIntermediateMethylatedQ20),"%.2f %%" %(self.getPercentage(nTotal,nTotalIntermediateMethylatedQ20))])
        vector.append(["UnMethylated","%i"%(nTotalUnmethylated),"%.2f %%" %(self.getPercentage(nTotal,nTotalUnmethylated))])
        vector.append(["UnMethylated Q>20","%i"%(nTotalUnmethylatedQ20),"%.2f %%" %(self.getPercentage(nTotal,nTotalUnmethylatedQ20))])
        
        return vector
        
    def getTableCGs(self):
        """ Returns a vector of CGs value """
        vector = []
        
        nTotal = self.data["TotalDinucleotides"]
        nNonReferenceCGs = self.data["TotalNonReferenceCpGs"]
        nNonReferenceCGsQ20 = self.data["TotalNonReferenceCpGsQuality20"]
        
        nCGsSNPs = self.data["TotalSnpsReferenceCpGsCalled"]
        nCGsSNPsQ20 = self.data["TotalSnpsReferenceCpGsCalledQuality20"]
        
        vector.append(["Non Reference CpGs","%i"%(nNonReferenceCGs),"%.2f %%" %(self.getPercentage(nTotal,nNonReferenceCGs))])
        vector.append(["Non Reference CpGs Q>20","%i"%(nNonReferenceCGsQ20),"%.2f %%" %(self.getPercentage(nTotal,nNonReferenceCGsQ20))])
        
        vector.append(["SNPs (CX) at Reference CpGs","%i"%(nCGsSNPs),"%.2f %%" %(self.getPercentage(nTotal,nCGsSNPs))])
        vector.append(["SNPs (CX) at Reference CpGs Q>20","%i"%(nCGsSNPsQ20),"%.2f %%" %(self.getPercentage(nTotal,nCGsSNPsQ20))])
        
        return vector
        
    def getTableStatusHomozygosity(self):
        """ Returns Vector Status Homozygosity """
        vector = []
        
        vector.append(["HomozygousMethylated",self.data["TotalHomozygousMethylated"]])
        vector.append(["HomozygousMethylated Q>20",self.data["TotalHomozygousMethylatedQuality20"]])
        
        vector.append(["HomozygousIntermediateMethylated",self.data["TotalHomozygousIntermediateMethylated"]])
        vector.append(["HomozygousIntermediateMethylated Q>20",self.data["TotalHomozygousIntermediateMethylatedQuality20"]])
        
        vector.append(["HomozygousUnMethylated",self.data["TotalHomozygousUnMethylated"]])
        vector.append(["HomozygousUnMethylated Q>20",self.data["TotalHomozygousUnMethylatedQuality20"]])
        
        return vector
        
    def getTableStatusAlternativeCX(self):
        """ Returns Vector Status Heterozygosity """
        vector = []
        
        vector.append(["AlternativeCXMethylated",self.data["TotalAlternativeCXMethylated"]])
        vector.append(["AlternativeCXMethylated Q>20",self.data["TotalAlternativeCXMethylatedQuality20"]])
        
        vector.append(["AlternativeCXIntermediateMethylated",self.data["TotalAlternativeCXIntermediateMethylated"]])
        vector.append(["AlternativeCXIntermediateMethylated Q>20",self.data["TotalAlternativeCXIntermediateMethylatedQuality20"]])
        
        vector.append(["AlternativeCXUnMethylated",self.data["TotalAlternativeCXUnMethylated"]])
        vector.append(["AlternativeCXUnMethylated Q>20",self.data["TotalAlternativeCXUnMethylatedQuality20"]])
        
        return vector
        
    def getTableStatusNonReference(self):
        """ Returns Vector Status Non Reference CGs """
        vector = []
        
        vector.append(["NonReferenceCpGsMethylated",self.data["TotalNonReferenceCpGsMethylated"]])
        vector.append(["NonReferenceCpGsMethylated Q>20",self.data["TotalNonReferenceCpGsMethylatedQuality20"]])
        
        vector.append(["NonReferenceCpGsIntermediateMethylated",self.data["TotalNonReferenceCpGsIntermediateMethylated"]])
        vector.append(["NonReferenceCpGsIntermediateMethylated Q>20",self.data["TotalNonReferenceCpGsIntermediateMethylatedQuality20"]])
        
        vector.append(["NonReferenceCpGsUnMethylated",self.data["TotalNonReferenceCpGsUnMethylated"]])
        vector.append(["NonReferenceCpGsUnMethylated Q>20",self.data["TotalNonReferenceCpGsUnMethylatedQuality20"]])
                
        return vector


    def getTableStatusSNPsReferenceCpGs(self):
        """ Returns Vector Status CpGs SNPs """
        vector = []
        
        vector.append(["SnpsReferenceCpGsCalledMethylated",self.data["TotalSnpsReferenceCpGsCalledMethylated"]])
        vector.append(["SnpsReferenceCpGsCalledMethylated Q>20",self.data["TotalSnpsReferenceCpGsCalledMethylatedQuality20"]])

        vector.append(["SnpsReferenceCpGsCalledIntermediateMethylated",self.data["TotalSnpsReferenceCpGsCalledIntermediateMethylated"]])
        vector.append(["SnpsReferenceCpGsCalledIntermediateMethylated Q>20",self.data["TotalSnpsReferenceCpGsCalledIntermediateMethylatedQuality20"]])

        vector.append(["SnpsReferenceCpGsCalledUnMethylated",self.data["TotalSnpsReferenceCpGsCalledUnMethylated"]])
        vector.append(["SnpsReferenceCpGsCalledUnMethylated Q>20",self.data["TotalSnpsReferenceCpGsCalledUnMethylatedQuality20"]])
        
        return vector

    def getVectorMethylationValues(self):
        """ Returns a vector of Methylation Values """        
        return self.methylation["MethValues"]

    def barplotCluster2(self,pngFile=None,vectorValues=None,title='Alternative CX Status',yLabel='#Dinucleotides'):
        """ Prints Clustered plot in a png file 
            Conceptual cluester plot
            
            pngFile - PNG file to store the data
            vectorValues - Vector of values to be printed as a bar plot
            title - Title of the plot
            yLabel - Y Label
        """
            
        wholeValues = (vectorValues[0][1],vectorValues[2][1],vectorValues[4][1])
        q20Values = (vectorValues[1][1],vectorValues[3][1],vectorValues[5][1])
        
        xLabels = ('Methylated','Intermediate','UnMethylated')
                   
        # the x locations for the groups
        groupsRange = [0, 3, 7]
        
        # the width of the bars
        width = 1
        
        figure, ax = plt.subplots()
        
        rects1 = ax.bar(groupsRange, wholeValues, width, color='r')
        rects2 = ax.bar([x+width for x in groupsRange], q20Values, width, color='y')
        
        # add some text for labels, title and axes ticks
        ax.set_ylabel(yLabel)
        ax.set_title(title)
        ax.set_xticks([ (x+width) for x in groupsRange] )
        ax.set_xticklabels(xLabels)
        
        ax.legend((rects1[0], rects2[0]), ('Total','Q>20'),loc='best')

        pylab.savefig(pngFile)
        
        plt.close(figure)
        
        
    def boxplot(self,pngFile=None):
        """ Prints a Boxplot in a png file 
            
            pngFile - PNG file to store the data
            vectorValues - Vector of values to be printed as a bar plot
        """
        #1. Plot Configuration
        figure, ax = plt.subplots(1,1)
        ax.set_title("Methylation Levels CGs Q>20")

        #2. Fitting Line for the histogram distribution
        #2.0 Vector of value
        x = self.getVectorMethylationValues()
        #2.1 Mean of distribution
        mu = np.mean(x)
        #2.2 Standard Deviation of distribution
        sigma = np.std(x)
        #2.3 Number of bins
        num_bins = 50
        #2.4 Data histogram
        n, bins,patches = ax.hist(x,num_bins,normed=1,facecolor='green',alpha=0.5)
        y = mlab.normpdf(bins,mu,sigma)
        ax.plot(bins,y,'r-')
        axes = plt.gca()
        ylim = axes.get_ylim()
                       
        #3. Axis Configuration
        ax.set_xlabel("Methylation Values")
        ax.set_ylabel("Density")
        
        #2.5 Horizontal Boxplot
        ax.boxplot(x,0,'b+',0, widths = ylim[-1]/4,positions=[ylim[-1]/2],manage_xticks=False) 
        ax.axes.set_ylim(ylim)

        pylab.savefig(pngFile)        
        plt.close(figure)
                

class HtmlCpgReport(HtmlVariantsReport):
    """ Builds Html CpG Report """
     
    def __init__(self):
        """ Html CpG Class Constructor """
        HtmlVariantsReport.__init__(self)
        
        
    def addRow(self,vTableHtml=None,concept=None,value=None,percentage=None,isOdd=True):
        """ Add Row html tags
            
            vectorHtml - Vector of tags to add
            concept - Name
            value - Value
            percentage - % Percentage Value
            isOdd - Boolean value
        """
        if isOdd == True: 
            vTableHtml.append("   <TR class=\"odd\">\n")
        else:
            vTableHtml.append("   <TR>\n")
        
        vTableHtml.append("   <TD> %s </TD> <TD> %s </TD> <TD> %s </TD> \n" %(concept,value,percentage))
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
        for field in cpgValues:    
            self.addRow(vTableHtml=vTableHtml,concept=field[0],value=field[1],percentage=field[2],isOdd=isOdd)
            isOdd = not isOdd
        
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
        vectorHtml.extend(self.createStatsTable(color='blue',cpgValues=cpgStats.getTableCpGs(),tableTitle="Dinucleotides Statistics",headers=["Type","#","%"]))
        
        #3.Methylation Status Table
        self.addSpaceSection(vectorHtml=vectorHtml) 
        vectorHtml.extend(self.createStatsTable(color='green',cpgValues=cpgStats.getTableMethylationStatus(),tableTitle="Methylation Status Statistics",headers=["Type","#","%"]))
        
        #4.Zigosity Table
        self.addSpaceSection(vectorHtml=vectorHtml) 
        vectorHtml.extend(self.createStatsTable(color='blue',cpgValues=cpgStats.getTableZygosity(),tableTitle="Zygosity Statistics",headers=["Type","#","%"]))
        
        #5.CGs Q>20 Methylation Status
        self.addSpaceSection(vectorHtml=vectorHtml) 
        boxPlot = "%s/IMG/%s.boxplot.png" %(output_dir,name)
        cpgStats.boxplot(pngFile=boxPlot)  
        vectorHtml.extend(self.buildPlot(color="green",tableTitle="CGs Q>20 Methylation Values",pathImage="IMG/%s"%(os.path.basename(boxPlot))))
                        
        #6. BarPlot Homozygous Methylation Status  
        self.addSpaceSection(vectorHtml=vectorHtml) 
        homozygousPlot = "%s/IMG/%s.homozygous.png" %(output_dir,name)
        cpgStats.barplotCluster2(pngFile=homozygousPlot,vectorValues=cpgStats.getTableStatusHomozygosity(),title='Homozygous Status',yLabel='#Dinucleotides')
        vectorHtml.extend(self.buildPlot(color="blue",tableTitle="Homozygous Status",pathImage="IMG/%s"%(os.path.basename(homozygousPlot))))
        
        #7. BarPlot AlternativeCX Methylation Status  
        self.addSpaceSection(vectorHtml=vectorHtml) 
        alternativeCXPlot = "%s/IMG/%s.alternativeCXPlot.png" %(output_dir,name)
        cpgStats.barplotCluster2(pngFile=alternativeCXPlot,vectorValues=cpgStats.getTableStatusAlternativeCX(),title='Alternative (CX) Status',yLabel='#Dinucleotides')
        vectorHtml.extend(self.buildPlot(color="green",tableTitle="Alternative (CX) Status",pathImage="IMG/%s"%(os.path.basename(alternativeCXPlot))))

        #8. NonReference and SNPs (CX) at Reference CpGs Tables
        self.addSpaceSection(vectorHtml=vectorHtml) 
        vectorHtml.extend(self.createStatsTable(color='blue',cpgValues=cpgStats.getTableCGs(),tableTitle="CpGs Statistics",headers=["Type","#","%"]))
    
        #9. BarPlot CpG NonReference Status 
        self.addSpaceSection(vectorHtml=vectorHtml) 
        statusNonReferencePlot = "%s/IMG/%s.status.cpg.nonreference.png" %(output_dir,name)         
        cpgStats.barplotCluster2(pngFile=statusNonReferencePlot,vectorValues=cpgStats.getTableStatusNonReference(),title='Non Reference CpGs Methylation Status',yLabel='#CpGs') 
        vectorHtml.extend(self.buildPlot(color="green",tableTitle="NonReference CpGs Status",pathImage="IMG/%s"%(os.path.basename(statusNonReferencePlot))))
    
        #10. BarPlot SNPs (CX) at Reference CpGs
        self.addSpaceSection(vectorHtml=vectorHtml) 
        cpgSnpPlots =  "%s/IMG/%s.status.cpg.snp.png" %(output_dir,name)
        cpgStats.barplotCluster2(pngFile=cpgSnpPlots,vectorValues=cpgStats.getTableStatusSNPsReferenceCpGs(),title='SNPs (CX) at Reference CpGs Methylation Status',yLabel='#CpGs') 
        vectorHtml.extend(self.buildPlot(color="blue",tableTitle="SNPs (CX) at Reference CpGs Status",pathImage="IMG/%s"%(os.path.basename(cpgSnpPlots))))


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
        
        
    def addRow(self,ident=20,lenCell=30,vectorSphinx=None,name_concept=None,statValue=None,percentValue=None):
        """ Add row sphinx
            
            ident - Characters to ident
            lenCell - Length of the Cell in each row
            vectorSphinx - Vector of tags to add
            name_concept - Name of the concept
            statValue - Value to be built
            percentValue - Percentage of values to be built
        """
        cells = 0
        cellContents = ""
        #Add ident chars
        for i in range(ident):
            cellContents += " "        
        
        cellContents += self.addCell(text="%s" %(name_concept),lenCell=lenCell,isLeft=True,isRight=False)
        cells += 1
        cellContents += self.addCell(text="%s" %(statValue),lenCell=lenCell,isLeft=False,isRight=False)
        cells += 1
        cellContents += self.addCell(text="%s" %(percentValue),lenCell=lenCell,isLeft=False,isRight=True)
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
        cellsHeader = 3
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
            self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,name_concept=rowValue[0],statValue=rowValue[1],percentValue=rowValue[2]) 
 
            
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
        vectorSphinx.extend(self.createStatsTable(ident=ident,lenCell=lenCell,vectorValues=cpgStats.getTableCpGs(),headers=["Type","#","%"]))
        vectorSphinx.append("\n")  
        
        #2.Methylation Status Table
        self.addSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Methylation Status Statistics")
        vectorSphinx.append("\n")
        vectorSphinx.extend(self.createStatsTable(ident=ident,lenCell=lenCell,vectorValues=cpgStats.getTableMethylationStatus(),headers=["Type","#","%"]))
        vectorSphinx.append("\n")  
        
        #3.Zigosity Table
        self.addSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Zygosity Statistics")
        vectorSphinx.append("\n")
        vectorSphinx.extend(self.createStatsTable(ident=ident,lenCell=lenCell,vectorValues=cpgStats.getTableZygosity(),headers=["Type","#","%"]))
        vectorSphinx.append("\n")
          
        #4. BoxPlot CGs Q>20 Methylation Status 
        self.addSubSubSection(ident=ident,vectorSphinx=vectorSphinx,title="CGs Q>20 Methylation Status")
        vectorSphinx.append("\n")
        boxPlot = "%s/IMG/%s_boxplot.png" %(output_dir,name)
        cpgStats.boxplot(pngFile=boxPlot)
        vectorSphinx.extend(self.buildImage(ident=0,pathImage="./IMG/%s"%(os.path.basename(boxPlot))))
        vectorSphinx.append("\n")
          
        #5. BarPlot Homozygous Methylation Status    
        self.addSubSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Homozygous Status")
        vectorSphinx.append("\n")
        homozygousPlot = "%s/IMG/%s_homozygous.png" %(output_dir,name)
        cpgStats.barplotCluster2(pngFile=homozygousPlot,vectorValues=cpgStats.getTableStatusHomozygosity(),title='Homozygous Status',yLabel='#Dinucleotides')
        vectorSphinx.extend(self.buildImage(ident=0,pathImage="./IMG/%s"%(os.path.basename(homozygousPlot))))
        vectorSphinx.append("\n")
              
        #6. BarPlot AlternativeCX Methylation Status  
        self.addSubSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Alternative (CX) Status")
        vectorSphinx.append("\n")
        alternativeCXPlot = "%s/IMG/%s_alternativeCX.png" %(output_dir,name)
        cpgStats.barplotCluster2(pngFile=alternativeCXPlot,vectorValues=cpgStats.getTableStatusAlternativeCX(),title='alternativeCX Status',yLabel='#Dinucleotides')
        vectorSphinx.extend(self.buildImage(ident=0,pathImage="./IMG/%s"%(os.path.basename(alternativeCXPlot))))
        vectorSphinx.append("\n")
        
        #7. NonReference and SNPs (CX) at Reference CpGss Tables
        self.addSubSection(ident=ident,vectorSphinx=vectorSphinx,title="CpGs Statistics")
        vectorSphinx.append("\n")
        vectorSphinx.extend(self.createStatsTable(ident=ident,lenCell=lenCell,vectorValues=cpgStats.getTableCGs(),headers=["Type","#","%"]))
        vectorSphinx.append("\n")

        #8. BarPlot CpG NonReference Status 
        self.addSubSection(ident=ident,vectorSphinx=vectorSphinx,title="NonReference CpGs Status")
        vectorSphinx.append("\n")
        statusNonReferencePlot = "%s/IMG/%s_status_cpg_nonreference.png" %(output_dir,name)         
        cpgStats.barplotCluster2(pngFile=statusNonReferencePlot,vectorValues=cpgStats.getTableStatusNonReference(),title='Non Reference CpGs Methylation Status',yLabel='#CpGs') 
        vectorSphinx.extend(self.buildImage(ident=0,pathImage="./IMG/%s"%(os.path.basename(statusNonReferencePlot))))
        vectorSphinx.append("\n")
            
        #9. BarPlot SNPs (CX) at Reference CpGs
        self.addSubSection(ident=ident,vectorSphinx=vectorSphinx,title="SNPs (CX) at Reference CpGs Status")
        vectorSphinx.append("\n")
        cpgSnpPlots =  "%s/IMG/%s_status_cpg_snp.png" %(output_dir,name)
        cpgStats.barplotCluster2(pngFile=cpgSnpPlots,vectorValues=cpgStats.getTableStatusSNPsReferenceCpGs(),title='SNPs (CX) at Reference CpGs Methylation Status',yLabel='#CpGs')
        vectorSphinx.extend(self.buildImage(ident=0,pathImage="./IMG/%s"%(os.path.basename(cpgSnpPlots))))
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
        vectorDocuments.append("SUMUP")   

          
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
                cpgStats.update(dictionary=json.load(file_json),methDict=json.load(meth_json))
            
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
                cpgStats.update(dictionary=json.load(file_json),methDict=json.load(meth_json)) 
            
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
        
        