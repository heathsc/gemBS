# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 12:51:01 2017

@author: marcos
"""

import json
import os

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.pylab as pylab

from matplotlib import cm
import matplotlib.colors as colors

#from mpl_toolkits.mplot3d import Axes3D

import importlib
importlib.import_module('mpl_toolkits.mplot3d').Axes3D

from reportStats import RunBasicStats
from report import BasicHtml
from sphinx import BasicSphinx
from sphinx import ConfigSphinx

class VariantStats(object):
    """ Basic definition of variant statistics """
    
    def __init__(self):
        """ Initializates a dictionary of basic values """

        self.data = {}
        
    def addValues(self,dictionary):
        """ Adds to previous data new set of values
        
            dictionary - Dictionary of values to be added
        """
        for key,value in dictionary.iteritems():
            if key in self.data:
                self.data[key] = self.data[key] + value
            else:
                self.data[key] = value
                
                
    def getBasicStatsTable(self):
        """ Gets Basic Stats From Variants Reports Table 
        
            returns a vector of values
        """
        collection_values = []
        
        collection_values.append(["SNPs","%i"%(self.data["TotalSnps"]),"%i"%(self.data["q20Snps"]),"%.2f %%" %(self.getPercentage(self.data["TotalSnps"],self.data["q20Snps"]))])
        #collection_values.append(["InDels","%i"%(self.data["TotalInDels"]),"%i"%(self.data["q20InDels"]),"%.2f %%" %(self.getPercentage(self.data["TotalInDels"],self.data["q20InDels"]))])
        collection_values.append(["Multiallelic","%i"%(self.data["TotalMultiallelic"]),"%i"%(self.data["q20Multiallelic"]),
                                  "%.2f %%" %(self.getPercentage(self.data["TotalMultiallelic"],self.data["q20Multiallelic"]))])
                                  
        return collection_values
                  

        
    def getPercentage(self,total,subtotal):
        """Calculates the percentage for a given value"""
        if total != 0.0:
            return (float(subtotal)/float(total))*100
        else:
            return 0.0
            
    def getTiTvRatio(self):
        """ Calculates Transition Transversion Ratio 
        
            return list of Ts/Tv, Transitions,Transversions        
        """
        transition = 0
        transversion = 0
        for key,value in self.data.iteritems():
            if key == "A>G" or key == "T>C" or key == "C>T" or key == "G>A":
                transition = transition + value
            elif key == "A>C" or key == "T>G" or key == "A>T" or key == "T>A" or key == "C>A" or key == "G>T" or key == "C>G" or key == "G>C" :
                transversion = transversion + value
                
        return [(float(transition) / float(transversion)),transition,transversion]
        
    
    def coveragePlot(self,fileOutput=None,concept=None,title=None,xScale=100,nBins=10):
        """Builds a Density Plots and print it to a given file
        
           fileOutput - File Ouput to print the plot
           concept - X label Concept
           title - title of the plot
           xScale - Size of its bin for the X axis
           nBins - Number of bins
        """
        # Dictionary of contens sorted by concept
        dictContents = {int(k):v for k,v in self.data.items()}
        
        #Create vector of Y values, where its position corresponds to the total number of values for a given concept. Ex: Concept = 10 Coverage Values = 120 Snps
        yValues = [0 for x in range(nBins)]
        
        for key,value in sorted(dictContents.iteritems()):
            index_for_key = int ((nBins*key)/(nBins*xScale))                        
            if (index_for_key >= nBins):
                index_for_key = nBins -1
            
            yValues[index_for_key] = yValues[index_for_key] + value
            
        #Create Vector of X ticks
        xValues = [xScale * x for x in range(nBins)]
            
        width = xScale
        matplotlib.pyplot.ioff()
        _, ax = plt.subplots()
        
        ax.bar(xValues, yValues, width, color="blue")
        ax.set_ylabel("#SNPs")
        ax.set_xlabel(concept)
        ax.set_title(title)

        pylab.savefig(fileOutput)
        
            
    def densityplot(self,fileOutput=None,scale=20,concept=None,title=None,bandwidth=0.45):
        """Builds a Density Plots and print it to a given file
        
           fileOutput - File Ouput to print the plot
           scale - X scale
           concept - X label Concept
           title - title of the plot
           xLim - X Limit To tune the plot size
        """
        # Create a density estimate from our data
        dictContents = {int(k):v for k,v in self.data.items()}
        dataList = [0 for i in range(sorted(dictContents.keys())[-1]+1)]

        #Get Values In dataList
        for key,value in sorted(dictContents.iteritems()):
            dataList[key] = value           
             
        
        density_est = gaussian_kde(dataList)

        # Control the 'smoothness'of the estimate. Higher values give smoother estimates.
        density_est.covariance_factor = lambda : bandwidth
        density_est._compute_covariance()
        
        i = 0
        maxim = len(dataList)
        
        xdata = []
       
        while i < maxim:
            xdata.append(i)
            i = scale + i
            
        # Density Plot   
        matplotlib.pyplot.ioff()

        _, ax = plt.subplots()

        ax.plot(xdata, density_est(xdata), color = '#539caf', lw = 2)
        ax.set_ylabel("Density")
        ax.set_xlabel(concept)
        ax.set_title(title)
        
        pylab.savefig(fileOutput)
        


        
        
class HtmlVariantsReport(object):
    """ Builds Html Variants Report """
    
    def addHtmlReportHeader(self,vectorHtml=None):
        ''' Add HTML Report Header to a given vector'''
        vectorHtml.append("<HTML>\n")

        vectorHtml.append(" <HEAD>\n")
        vectorHtml.append(" <STYLE TYPE=\"text/css\">\n")
        vectorHtml.append("  <!--\n")
        vectorHtml.append("   @import url(\"style.css\"); \n")
        vectorHtml.append("  -->\n")
        vectorHtml.append(" </STYLE>\n")
        vectorHtml.append(" </HEAD>\n")
        vectorHtml.append(" <BODY>\n")
        
    def closeHtmlReport(self,vectorHtml=None):
        ''' Close HTML header to a given vector object '''  
        vectorHtml.append(" </BODY>\n")
        vectorHtml.append("</HTML>\n")
        
    def addSpaceSection(self,vectorHtml=None):
        ''' Add HTML Report Header to a given vector'''
        vectorHtml.append("<BR><BR><BR>\n")
        
    def addHtmlStartReport(self,vectorHtml=None,parentName=None,currentName=None,parentDocument=None):
        ''' Add HTML Report Header to a given vector'''
        self.addHtmlReportHeader(vectorHtml)

        vectorHtml.append("\n <P id='path'> /%s/%s </P> \n" %(parentName,currentName))
        vectorHtml.append("\n <a class=\"link\" href=\"%s\"><B>BACK</B></a> <br>\n" %(os.path.basename(parentDocument)))   
        vectorHtml.append("  <H1 id=\"title\"> <U> SAMPLE %s </U> </H1>\n" %(currentName))
        
        
    def addRow(self,vectorHtml=None,value=None,name_concept=None,isOdd=True):
        """ Add Row html tags
            
            vectorHtml - Vector of tags to add
            value - stats value
            name_concept - Name of the concept
        """
        if isOdd == True: 
            vectorHtml.append("   <TR class=\"odd\">\n")
        else:
            vectorHtml.append("   <TR>\n")
        
        vectorHtml.append("   <TD> %s </TD> <TD> %s </TD> \n" %(name_concept,value))
        vectorHtml.append("   </TR>\n")
        

    def addRowDynamic(self,vectorHtml=None,values=None,isOdd=True):
        """ Add Row html tags
            
            vectorHtml - Vector of tags to add
            values - stats value
            name_concept - Name of the concept
        """
        if isOdd == True: 
            vectorHtml.append("   <TR class=\"odd\">\n")
        else:
            vectorHtml.append("   <TR>\n")
        
        for field in values:        
            vectorHtml.append("   <TD> %s </TD> \n" %(field))
            
        vectorHtml.append("   </TR>\n")        
        
    def tableDefinition(self,vectorHtml=None,color=None):
        """ HTML Table tag definition
        
        vectorHtml - List of html to add table tags
        color - Table color, could be green or blue
        """
        if color == "green":
            vectorHtml.append("  <TABLE id=\"green\">\n")
        else:
            vectorHtml.append("  <TABLE id=\"hor-zebra\">\n")
        
    def createStatsTable(self,color=None,variantsValue=None,tableTitle=None,headers=None):
        """Create Statistics Table
        
           color - Table color, could be green or blue  
           variantsValue - Dictionary of values to be printed
           tableTitle - Table Title
           headers - list of values to be printed as headers
           returns vector of HTML tags
        """
        vTableHtml = []   
        
        #0.Table Title
        vTableHtml.append('<H3 id="section"> %s </H3>\n' %(tableTitle))
     
        #1.Table Definition
        self.tableDefinition(vectorHtml=vTableHtml,color=color)

        #2.Table header
        vTableHtml.append("   <TR>\n")
        for conceptHeader in headers:        
            vTableHtml.append("    <TH scope=\"col\">%s</TH>\n"%(conceptHeader))
            
        vTableHtml.append("   </TR>\n")
        
        #3.Fullfill table
        isOdd = True
        for key,value in variantsValue.iteritems():        
            self.addRow(vTableHtml,value,key,isOdd)
            isOdd = not isOdd
        
        vTableHtml.append(" </TABLE>\n")
        return vTableHtml
        
        
    def createGeneralTable(self,color=None,values=None,tableTitle=None,headers=None):
        """Create Statistics Table
        
           color - Table color, could be green or blue  
           values - Vector of values to be printed
           tableTitle - Table Title
           headers - list of values to be printed as headers
           returns vector of HTML tags
        """
        vTableHtml = []   
        
        #0.Table Title
        vTableHtml.append('<H3 id="section"> %s </H3>\n' %(tableTitle))
     
        #1.Table Definition
        self.tableDefinition(vectorHtml=vTableHtml,color=color)

        #2.Table header
        vTableHtml.append("   <TR>\n")
        for conceptHeader in headers:        
            vTableHtml.append("    <TH scope=\"col\">%s</TH>\n"%(conceptHeader))
            
        vTableHtml.append("   </TR>\n")
        
        #3.Fullfill table
        isOdd = True
        for fields in values:    
            self.addRowDynamic(vectorHtml=vTableHtml,values=fields,isOdd=isOdd)
            isOdd = not isOdd
        
        vTableHtml.append(" </TABLE>\n")
        return vTableHtml
   
        
    def buildPlot(self,color=None,tableTitle=None,pathImage=None):
        """Build Plot and locate on html report
        
          color - Color for the table were the plot is going to be located
          tableTitle - Table Title
          pathImage - Path to the image to be shown
          returns Vector of HTML tags
        """
        vTableHtml = []
        
        #1. HTML TAG BUILDING 
        vTableHtml = []
   
        self.tableDefinition(vectorHtml=vTableHtml,color=color)

        vTableHtml.append("   <TR class=\"odd\"> <TH scope=\"col\"> %s </TH> </TR> \n" %(tableTitle))
        vTableHtml.append('   <TR> <TD> <img src="%s" alt="%s"> </TD> </TR> \n' %(pathImage,pathImage))
        vTableHtml.append("  </TABLE>\n")
        
        return vTableHtml
        
        
    def buildTwoPlots(self,color=None,tableTitleFirst=None,tableTitleSecond=None,pathFirstImage=None,pathSecondImage=None):
        """Build Plot and locate on html report
        
          color - Color for the table were the plot is going to be located
          tableTitleFirst - Table First Title
          tableTitleSecond - Table Second Title
          pathFirstImage - Path to the first image to be shown
          pathSecondImage - Path to the second image to be shown
          returns Vector of HTML tags
        """
        vTableHtml = []
        
        #1. HTML TAG BUILDING 
        vTableHtml = []
   
        self.tableDefinition(vectorHtml=vTableHtml,color=color)

        vTableHtml.append("   <TR class=\"odd\"> <TH scope=\"col\"> %s </TH>  <TH scope=\"col\"> %s </TH> </TR> \n" %(tableTitleFirst,tableTitleSecond))
        vTableHtml.append('   <TR> <TD> <img src="%s" alt="%s"> </TD> <TD> <img src="%s" alt="%s"> </TD> </TR> \n' %(pathFirstImage,pathFirstImage,pathSecondImage,pathSecondImage))
        vTableHtml.append("  </TABLE>\n")
        
        return vTableHtml
        
        
        
    def buildThreePlots(self,color=None,tableTitleFirst=None,tableTitleSecond=None,tableTitleThird=None,
                       pathFirstImage=None,pathSecondImage=None,pathThirdImage=None):
        """Build Plot and locate on html report
        
          color - Color for the table were the plot is going to be located
          tableTitleFirst - Table First Title
          tableTitleSecond - Table Second Title
          tableTitleThird - Table Third Title
          pathFirstImage - Path to the first image to be shown
          pathSecondImage - Path to the second image to be shown
          pathThirdImage - Path to the third image to be shown
          returns Vector of HTML tags
        """
        vTableHtml = []
        
        #1. HTML TAG BUILDING 
        vTableHtml = []
   
        self.tableDefinition(vectorHtml=vTableHtml,color=color)

        vTableHtml.append("   <TR class=\"odd\"> <TH scope=\"col\"> %s </TH>  <TH scope=\"col\"> %s </TH> <TH scope=\"col\"> %s </TH> </TR> \n" %(tableTitleFirst,tableTitleSecond,tableTitleThird))
        vTableHtml.append('   <TR> <TD> <img src="%s" alt="%s"> </TD> <TD> <img src="%s" alt="%s"> </TD>  <TD> <img src="%s" alt="%s"> </TD> </TR> \n' 
                          %(pathFirstImage,pathFirstImage,pathSecondImage,pathSecondImage,pathThirdImage,pathThirdImage))
        vTableHtml.append("  </TABLE>\n")
        
        return vTableHtml
        
        
        
    def run(self,vectorHtml=None,listOfConcepts=None,output_dir=None,name=None,parentName=None,parentDocument=None):
        """ Run Html Variants Report 
        
            vectorHtml - Vector of html tags
            listOfConcepts - List Of VariantStats Concpets to be printed
            ouput_dir - Ouput Directory to store the report
            name - Basic name
            parentName - Parent Name Document
            parentDocument - HTML Parent Document
        """    
        
        #1. Add HTML Report Header
        self.addHtmlStartReport(vectorHtml=vectorHtml,parentName=parentName,currentName=name,parentDocument=parentDocument)
        
        #2. Add Basic Stats Table    
        self.addSpaceSection(vectorHtml=vectorHtml) 
        vectorHtml.extend(self.createGeneralTable(color='blue',values=listOfConcepts[0].getBasicStatsTable(),tableTitle="Variants Statistics",headers=["Type","Total","Q>20","%"]))       
              
        #3. Add Density plot Coverage   
        self.addSpaceSection(vectorHtml=vectorHtml) 
        coveragePlot = "%s/IMG/%s.coverage.png" %(output_dir,name)
        listOfConcepts[1].coveragePlot(fileOutput=coveragePlot,concept="Depth of Coverage",title="Variants Coverage",xScale=10,nBins=20)
        #listOfConcepts[1].densityplot(fileOutput=coveragePlot,scale=20,concept="Depth of Coverage",title="",bandwidth=0.45)        
        vectorHtml.extend(self.buildPlot(color="green",tableTitle="Variants Coverage",pathImage="IMG/%s"%(os.path.basename(coveragePlot))))
                
        #4. Genotype Quality plot Coverage   
        self.addSpaceSection(vectorHtml=vectorHtml) 
        genotypeQualityPlot = "%s/IMG/%s.genotype.qualities.png" %(output_dir,name)
        #listOfConcepts[2].densityplot(fileOutput=genotypeQualityPlot,scale=20,concept="Genotype Quality",title="",bandwidth=0.01)  
        listOfConcepts[2].coveragePlot(fileOutput=genotypeQualityPlot,concept="Genotype Quality",title="Variants Qualities",xScale=10,nBins=26)        
        vectorHtml.extend(self.buildPlot(color="blue",tableTitle="Variants Qualities",pathImage="IMG/%s"%(os.path.basename(genotypeQualityPlot))))
        
        #5. Mutation Profile    
        self.addSpaceSection(vectorHtml=vectorHtml) 
        vectorHtml.extend(self.createStatsTable(color='green',variantsValue=listOfConcepts[3].data,tableTitle="Mutation Profile",headers=["Type","Total"]))

        #6. Mutation Profile    
        self.addSpaceSection(vectorHtml=vectorHtml) 
        vectorHtml.extend(self.createStatsTable(color='green',variantsValue=listOfConcepts[4].data,tableTitle="Mutation Profile Q>20",headers=["Type","Total"]))
        
        #7. Transition / Transversion Ratio
        tiTvRatioList = listOfConcepts[3].getTiTvRatio() 
        tiTvRatioDict = {}
        tiTvRatioDict["Ts/Tv"] = tiTvRatioList[0] 
        tiTvRatioDict["Ts"] = tiTvRatioList[1] 
        tiTvRatioDict["Tv"] = tiTvRatioList[2] 
        self.addSpaceSection(vectorHtml=vectorHtml) 
        vectorHtml.extend(self.createStatsTable(color='blue',variantsValue=tiTvRatioDict,tableTitle="Transition Transversion",headers=["Type","Total"]))
        
        #8. Transition / Transversion Ratio
        tiTvRatioList = listOfConcepts[4].getTiTvRatio() 
        tiTvRatioDict = {}
        tiTvRatioDict["Ts/Tv"] = tiTvRatioList[0] 
        tiTvRatioDict["Ts"] = tiTvRatioList[1] 
        tiTvRatioDict["Tv"] = tiTvRatioList[2] 
        self.addSpaceSection(vectorHtml=vectorHtml) 
        vectorHtml.extend(self.createStatsTable(color='blue',variantsValue=tiTvRatioDict,tableTitle="Transition Transversion Q>20",headers=["Type","Total"]))
        
      
class IndexVariantsHtml(HtmlVariantsReport):
    """ Class which defines Index Sample Report """
    
    def __init__(self,output_dir=None,name_project=None,dictionary_samples=None):
        """  Class constructor
        
             output_dir -- Output directory to store HTML reports
             name_project -- Project name
             dictionary_samples -- Dictionary of samples and list of stats objects. [sample][ListStatsObject]
        """
        #Call Parent class constructor
        HtmlVariantsReport.__init__(self)        
        
        self.output_dir = output_dir
        self.name_project = name_project
        self.variants_html_document = "%s/%s.html" %(self.output_dir,self.name_project)
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
        
        vectorHtml.append('<H1 id="title"> <U> Variants Report from Methylation Pipeline Project %s </U> </H1>' %(self.name_project))
        
        vector_sample_links = []
        
        for sample, listStats in sorted(self.dictionary_samples.iteritems()):            
            sampleHtml = "%s/%s_variants.html" %(self.output_dir,sample)
            
            vSampleHtml = []            
            sample_report = HtmlVariantsReport()
            sample_report.run(vectorHtml=vSampleHtml,listOfConcepts=listStats,
                              output_dir=self.output_dir,name=sample,
                              parentName=self.name_project,parentDocument=self.variants_html_document)            
            
            self.closeHtmlReport(vectorHtml=vSampleHtml)
            RunBasicStats.saveDocument(file_name=sampleHtml,vectorContent=vSampleHtml)
            vector_sample_links.append([os.path.basename(sampleHtml),sample])
            
            
        #2.Links table
        self.createLinksTables(vector_html=vectorHtml,vector_links=vector_sample_links,title="Sample Reports",color="blue")
        
        #3.Close Html
        self.closeHtmlReport(vectorHtml=vectorHtml)
                  
      
class SphinxVariantsReport(BasicSphinx):
    """ Class responsable for printing Sphinx report per Lane """
    
    def __init__(self):
        """Sphinx Variants Report building"""
        
        #Call Parent class constructor
        BasicSphinx.__init__(self,mapping_stats=None,png_mapq_histogram="",png_insert_size_histogram="")        


    def addRow(self,ident=20,lenCell=30,vectorSphinx=None,statValue=None,name_concept=None):
        """ Add row sphinx
            
            ident - Characters to ident
            lenCell - Length of the Cell in each row
            vectorSphinx - Vector of tags to add
            statValue - Value to be built
            name_concept - Name of the concept
        """
        cells = 0
        cellContents = ""
        #Add ident chars
        for i in range(ident):
            cellContents += " "        
        
        cellContents += self.addCell(text="%s" %(name_concept),lenCell=lenCell,isLeft=True,isRight=False)
        cells += 1
        cellContents += self.addCell(text="%s" %(statValue),lenCell=lenCell,isLeft=False,isRight=True)
        cells += 1
        
        vectorSphinx.append(cellContents) 
            
        #Underline
        self.addSimpleLine(ident=ident,vectorSphinx=vectorSphinx,cells=cells,lenCell=lenCell)


    def addRowDynamic(self,ident=20,lenCell=30,vectorSphinx=None,rowValues=None): 
        """ Add row sphinx
            
            ident - Characters to ident
            lenCell - Length of the Cell in each row
            vectorSphinx - Vector of tags to add
            rowValues - Vector of values to be printed in a row
        """
        cells = 0
        cellContents = ""
        #Add ident chars
        for i in range(ident):
            cellContents += " " 
 
        lenFields = len(rowValues)
 
        for field in rowValues:
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


    def createStatsTable(self,ident=20,lenCell=30,dictValues=None,headers=None):
        """Create Statistics Table
        
           ident - Characters to ident 
           lenCell - Cell length
           dictValues - Dictionary of key values
           headers - List Of Headers
           returns vector of Sphinx tags
        """
        vTableSphinx = []   
        cellsHeader = 0
        
        #1.Table header
        cellsHeader = 2
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
        for key,value in dictValues.iteritems():        
            self.addRow(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,statValue=value,name_concept=key)    
               
        return vTableSphinx
        

    def createGeneralTable(self,ident=20,lenCell=30,values=[],headers=[]):
        """ Create General Variants Table

           ident - Characters to ident 
           lenCell - Cell length
           values - Array of values
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
        for rowValues in values:        
            self.addRowDynamic(ident=ident,lenCell=lenCell,vectorSphinx=vTableSphinx,rowValues=rowValues)    
               
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



    def run(self,vectorSphinx=None,ident=20,lenCell=30,listOfConcepts=None,output_dir=None,name=None):
        """ Run Lane SPHINX Documentation Building
        
            vectorSphinx - vector of Sphinx tags
            ident - Characters to ident
            lenCell - Length of table cells in characters
            listOfConcepts - List of variant report concepts
            ouput_dir - Ouput Directory to store the report
            name - Basic name
        """  
        #0. Add Chapter Line       
        self.addTopSection(ident=ident,vectorSphinx=vectorSphinx,title="Variants Statistics for %s" %(name))        
        vectorSphinx.append("\n")        
        
        #1. Add Basic Stas Table
        self.addSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Basic Stats")
        vectorSphinx.append("\n")
        vectorSphinx.extend(self.createGeneralTable(ident=ident,lenCell=lenCell,values=listOfConcepts[0].getBasicStatsTable(),headers=["Type","Total","Q>20","%"]))        
        vectorSphinx.append("\n")
        
        #4. Add Density plot Coverage   
        self.addSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Coverage density plot")
        vectorSphinx.append("\n")
        coveragePlot = "%s/IMG/%s_coverage.png" %(output_dir,name)
        #listOfConcepts[1].densityplot(fileOutput=coveragePlot,scale=20,concept="Depth of Coverage",title="",bandwidth=0.45)
        listOfConcepts[1].coveragePlot(fileOutput=coveragePlot,concept="Depth of Coverage",title="Variants Coverage",xScale=10,nBins=20)
        vectorSphinx.extend(self.buildImage(ident=20,pathImage="./IMG/%s"%(os.path.basename(coveragePlot))))
        vectorSphinx.append("\n")
        
        #5. Genotype Quality plot Coverage   
        self.addSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Genotype Quality density plot")
        vectorSphinx.append("\n")
        genotypeQualityPlot = "%s/IMG/%s_genotype_qualities.png" %(output_dir,name)
        #listOfConcepts[2].densityplot(fileOutput=genotypeQualityPlot,scale=20,concept="Genotype Quality",title="",bandwidth=0.01) 
        listOfConcepts[2].coveragePlot(fileOutput=genotypeQualityPlot,concept="Genotype Quality",title="Variants Qualities",xScale=10,nBins=26)
        vectorSphinx.extend(self.buildImage(ident=20,pathImage="./IMG/%s"%(os.path.basename(genotypeQualityPlot))))
        vectorSphinx.append("\n")
        
        #6. Mutation Profile    
        self.addSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Mutation Profile")
        
        self.addSubSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Total")
        vectorSphinx.append("\n")
        vectorSphinx.extend(self.createStatsTable(ident=ident,lenCell=lenCell,dictValues=listOfConcepts[3].data,headers=["Type","Total"]))
        vectorSphinx.append("\n")

        #7. Mutation Profile    
        self.addSubSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Quality >20")
        vectorSphinx.append("\n")
        vectorSphinx.extend(self.createStatsTable(ident=ident,lenCell=lenCell,dictValues=listOfConcepts[4].data,headers=["Type","Total"]))
        vectorSphinx.append("\n")

        #8. Transition / Transversion Ratio
        tiTvRatioList = listOfConcepts[3].getTiTvRatio() 
        tiTvRatioDict = {}
        tiTvRatioDict["Ts/Tv"] = "%.2f" %(tiTvRatioList[0]) 
        tiTvRatioDict["Ts"] = tiTvRatioList[1] 
        tiTvRatioDict["Tv"] = tiTvRatioList[2] 
        self.addSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Transition Transversion")   
        
        self.addSubSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Total")
        vectorSphinx.append("\n")
        vectorSphinx.extend(self.createStatsTable(ident=ident,lenCell=lenCell,dictValues=tiTvRatioDict,headers=["Type","Total"]))
        vectorSphinx.append("\n")  
        
        #9. Transition / Transversion Ratio
        tiTvRatioList = listOfConcepts[4].getTiTvRatio() 
        tiTvRatioDict = {}
        tiTvRatioDict["Ts/Tv"] = "%.2f" %(tiTvRatioList[0]) 
        tiTvRatioDict["Ts"] = tiTvRatioList[1] 
        tiTvRatioDict["Tv"] = tiTvRatioList[2] 
        self.addSubSubSection(ident=ident,vectorSphinx=vectorSphinx,title="Quality >20")

        vectorSphinx.append("\n")
        vectorSphinx.extend(self.createStatsTable(ident=ident,lenCell=lenCell,dictValues=tiTvRatioDict,headers=["Type","Total"]))
        vectorSphinx.append("\n")        
        

class SumupVariantsSphinx(BasicSphinx):
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
        self.variants_sphinx_document = "%s/%s.rst" %(self.output_dir,self.name_project)
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
        
        self.addPartLine(ident=0,vectorSphinx=vectorSphinx,title="Variants Report from Methylation Pipeline Project %s"%(scaped_name))

        vectorDocuments = []
        vectorDocuments.append("SUMUP")        
          
        #2. Create Sample and lanes Report 
        for sample, listStats in sorted(self.dictionary_samples.iteritems()):            
            sampleSphinx = "%s/%s.rst" %(self.output_dir,sample)
            vectorDocuments.append("%s" %(sample))
       
            #Create Build SPHINX
            vector_sphinx = []
            sphinxReport = SphinxVariantsReport()
            sphinxReport.run(vectorSphinx=vector_sphinx,ident=0,lenCell=30,
                             listOfConcepts=listStats,output_dir=self.output_dir,name=sample)

            #Store Sphinx file
            RunBasicStats.saveDocument(file_name=sampleSphinx,vectorContent=vector_sphinx)
             
        #3.Toc Table
        vectorSphinx.append("\n")
        self.addTocTree(vectorSphinx=vectorSphinx,vectorDocuments=vectorDocuments)
      
      
def buildVariantReport(inputs=None,output_dir=None,name=None):
    """ Build variant report.
    
        inputs -- Dictionary of samples and list of files. [sample] [barcode_chrN.json,barcode_chrN+1.json]
        output_dir -- Output directory to store html documents.
        name --  Name basic to build output results.
    """

    #Check output directory
    if not os.path.exists("%s/IMG/" %(output_dir)):
        os.makedirs("%s/IMG/" %(output_dir)) 

    #Proces list chromosome files
    dict_samples = {}
       
    for sample,chrom_json_files in inputs.iteritems():
        #Parsing json file
        basicStats = VariantStats()
        q20BasicStats = VariantStats()
        coverageStats = VariantStats()
        genotypeQualityStats = VariantStats()
        mutationChangesStats = VariantStats()
        mutationChangesQ20Stats = VariantStats()
        chromosomeStats = VariantStats()

        #Load al json files
        for json_file in chrom_json_files:    
            with open(json_file, 'r') as file_json:
                data = json.load(file_json)
                basicStats.addValues({"TotalSnps":data["TotalSnps"],"TotalInDels":data["TotalIndels"],"TotalMultiallelic":data["TotalMultiallelic"],
                                      "q20Snps":data["q20Snps"],"q20InDels":data["q20Indels"],"q20Multiallelic":data["q20Multiallelic"] })
                                      
                coverageStats.addValues(data["coverageVariants"])
                genotypeQualityStats.addValues(data["qualityVariants"])
                mutationChangesStats.addValues(data["mutations"])
                mutationChangesQ20Stats.addValues(data["mutationsQ20"])
                chromosomeStats.addValues(data["chromosomeVariants"])
                
        #Create Build HTML
        dict_samples[sample] = [basicStats,coverageStats,genotypeQualityStats,mutationChangesStats,mutationChangesQ20Stats,chromosomeStats]

    #IndexHtml object
    vector_index_html = []
    indexVariantsHtml = IndexVariantsHtml(output_dir=output_dir,name_project=name,dictionary_samples=dict_samples)
    indexVariantsHtml.run(vectorHtml=vector_index_html)
    RunBasicStats.saveDocument(file_name=indexVariantsHtml.variants_html_document,vectorContent=vector_index_html)
    
    #CSS object
    cssBuilder = BasicHtml()
    vector_css = []
    cssBuilder.buildStyleSheet(vector_css)
    RunBasicStats.saveDocument(file_name="%s/style.css" %(output_dir),vectorContent=vector_css) 
    
    
    
    
def buildSphinxVariantReport(inputs=None,output_dir=None,name=None):
    """ Build Sphinx Variant report
    
        inputs -- List of json files
        output_dir -- Output directory to store html documents.
        name --  Name basic to build output results.
    """

    #Check output directory
    if not os.path.exists("%s/IMG/" %(output_dir)):
        os.makedirs("%s/IMG/" %(output_dir))  
        

    #Proces list chromosome files
    dict_samples = {}
    for sample,chrom_json_files in inputs.iteritems():
        #Parsing json file
        basicStats = VariantStats()
        q20BasicStats = VariantStats()
        coverageStats = VariantStats()
        genotypeQualityStats = VariantStats()
        mutationChangesStats = VariantStats()
        mutationChangesQ20Stats = VariantStats()
        chromosomeStats = VariantStats()

        #Load al json files
        for json_file in chrom_json_files:    
            with open(json_file, 'r') as file_json:
                data = json.load(file_json)
                basicStats.addValues({"TotalSnps":data["TotalSnps"],"TotalInDels":data["TotalIndels"],"TotalMultiallelic":data["TotalMultiallelic"],
                                      "q20Snps":data["q20Snps"],"q20InDels":data["q20Indels"],"q20Multiallelic":data["q20Multiallelic"] })           
                                
                coverageStats.addValues(data["coverageVariants"])
                genotypeQualityStats.addValues(data["qualityVariants"])
                mutationChangesStats.addValues(data["mutations"])
                mutationChangesQ20Stats.addValues(data["mutationsQ20"])
                chromosomeStats.addValues(data["chromosomeVariants"])
                
        #Create Build HTML
        dict_samples[sample] = [basicStats,coverageStats,genotypeQualityStats,mutationChangesStats,mutationChangesQ20Stats,chromosomeStats]




    #SumupSphinx object
    vector_sumup_sphinx = []
    sumupVariantsSphinx = SumupVariantsSphinx(output_dir=output_dir,name_project=name,dictionary_samples=dict_samples)
    sumupVariantsSphinx.run(vectorSphinx=vector_sumup_sphinx)
    RunBasicStats.saveDocument(file_name=sumupVariantsSphinx.variants_sphinx_document,vectorContent=vector_sumup_sphinx)
    
    #Config python file
    cfgFile = ConfigSphinx(path_config_file="%s/conf.py" %(output_dir),path_makefile_file="%s/Makefile" %(output_dir),master_file=name,project_name=name,main_title='VARIANTS REPORT')
    cfgFile.run()


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
