/*
 * main.c
 *
 *  Created on: 20 Mai, 2016
 *      Author: marcos
 */
#include <stdio.h>
#include <unistd.h>
#include "parseArgs.h"
#include "parseInput.h"
#include "methBed.h"
#include "counts.h"
#include "intersection.h"
#include <string.h>


/**
 * \brief Run Bed annotation
 * \param cpgInputFile CpG Input file
 * \param bedFile BedFile to be methylation annotated
 * \param annotatedFile File output annotation
 * \param jsonFile JSON format input file
 * \param isGzip 1 if input is zipped otherwise 0
 * \returns return value
 */
int bedAnnotation(char * cpgInputFile,char * bedFile,char * annotatedFile,char * jsonFile, int isGzip)
{
	/*1. Configure Input Data*/
	if (configInput(cpgInputFile,isGzip) < 1)
	{
	    return 0;
	}

	/*2. Configure Bed Input File*/
	if(!checkInputBed(bedFile))
	{
		return 0;
	}

	/*3. Initializations of counts */
	initCounts();

	/*4. Parse Each Window BED */
	struct Bed bedWindow;
	unsigned int homozygous;
	Dinucleotide * headDino, * currentDino;
	headDino = NULL;
	currentDino = NULL;

	if (checkFileOutput(annotatedFile) != 1)
	{
		return 0;
	}

	while (getNewWindow(&bedWindow) != 0)
	{
        struct Record record;
        while (getNewRecord(&record,isGzip) != 0)
        {
        	/*4.1.0 Add Record Stats*/
        	addRecordStats(&record);
        	/*4.1.1 If record is in windows range then add it to the memory heap*/
        	if( (strcmp(record.contig,bedWindow.contig)==0) && (record.position >= bedWindow.start) && (record.position <= bedWindow.end))
        	{
        		homozygous = 1;
        		if(strcmp(record.callContext,record.referenceContext) != 0)
        		{
        			homozygous = 0;
        		}
        		newNode(&headDino,&currentDino,record.contig,record.position,record.methValue,homozygous);
        	}
        	/*4.1.2 If record is over windows range then add it to the memory heap and stop looping over dinucleotides record*/
        	else if( (strcmp(record.contig,bedWindow.contig)!=0) || ( (strcmp(record.contig,bedWindow.contig)==0) && (record.position >= bedWindow.end)) )
        	{
           		homozygous = 1;
				if(strcmp(record.callContext,record.referenceContext) != 0)
				{
					homozygous = 0;
				}
        		newNode(&headDino,&currentDino,record.contig,record.position,record.methValue,homozygous);
				break;
        	}
        }
        /*4.2 Get Window Methylation */
        getWindowMethylation(headDino, &bedWindow);
        /*4.3 Remove Dinucleotides in memory that are not useful*/
        toRemove(&headDino,&bedWindow);
        /*4.4 Add window bed File*/
        addBedWindow(&bedWindow);
	}

	/*5. Close file descriptors*/
	closeInput(isGzip);
	closeFileOutput();

	/*6. PRINT RESULTS */
	printCounts();

	/*7. JSON FILE OUTPUT */
	if(jsonFile != NULL)
	{
		saveCounts(jsonFile);
	}

	return 1;
}

/**
 * \brief Get Statistics from CpG Input file
 * \param cpgInputFile CpG Input file
 * \param jsonFile JSON format input file
 * \param methJsonFile Methylation Values JSON output file
 * \param isGzip 1 if input is zipped otherwise 0
 * \return 1 if everything goes well otherwise 0
 */
int getStats(char * cpgInputFile, char * jsonFile,char * methJsonFile, int isGzip)
{
	/*1. READ INPUT DATA LINE TO LINE TO GET RECORDS*/
	struct Record record;

	/*1. Configure Input Data*/
	if (configInput(cpgInputFile,isGzip) < 1)
	{
		return 0;
	}

	/*1.1 Initializations of counts*/
	initCounts();

	while (getNewRecord(&record,isGzip)!=0)
	{
		addRecordStats(&record);
	}

	/*2. CLOSE FILE*/
	closeInput(isGzip);

	/*3. PRINT RESULTS */
	printCounts();

	/*3.1 JSON FILE OUTPUT */
	if(jsonFile != NULL)
	{
		saveCounts(jsonFile);
	}

	/*4. METHYLATION JASON FILE VALUES */
	if(jsonFile != NULL)
	{
		saveJsonMethylationCounts(methJsonFile);
	}

	return 1;
}

int main (int argc, char *argv[])
{
    /*0. PARSE ARGUMENTS*/
	struct Args arguments;

	if (getArgs(&arguments,argc,argv) < 1)
	{
		return 1;
	}

	/*1. GET STATS*/
	if(arguments.bedFile == NULL && arguments.cpgInputFile != NULL)
	{
		if (getStats(arguments.cpgInputFile,arguments.jsonFile,arguments.methJsonFile,arguments.isZipped) < 1)
		{
			return 1;
		}
	}

	/*5. BED FILE METHYLATION ANNOTATION*/
	if (arguments.bedFile != NULL && arguments.annotatedFile != NULL && arguments.cpgInputFile != NULL)
	{
		return bedAnnotation(arguments.cpgInputFile,arguments.bedFile,arguments.annotatedFile,arguments.jsonFile,arguments.isZipped);
	}
	else if( (arguments.bedFile != NULL && arguments.annotatedFile == NULL) || (arguments.bedFile == NULL && arguments.annotatedFile != NULL) )
	{
		printf("Sorry! Not mandatory arguments for bed annotation has been specified. Please use -b file.sorted.bed and -a file.annotated.bed \n");
	}

	/*6. INTERSECTION OF DINUCLEOTIDE FILES*/
	if (arguments.firstCpGIsecFile != NULL && arguments.secondCpGIsecFile != NULL)
	{
		return runIsec(arguments.firstCpGIsecFile, arguments.secondCpGIsecFile,arguments.areIsecZipped);
	}

	return 0;
}
