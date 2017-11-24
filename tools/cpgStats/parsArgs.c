/*
 * parsArgs.c
 *
 *  Created on: 20 Mai, 2016
 *      Author: marcos
 */

#include "parseArgs.h"
#include <stdio.h>
#include <unistd.h>

/**
 * \brief Print Program help
 */
void printHelp()
{
	printf("cpgStats Parses a CpG Dinucleotide file and outputs its statistical results. \n");
	printf("It can produce a json output and annotate a bed file. \n");
	printf("Two CpG Dinucleotides files could be compared using Intersection parameters. \n");
	printf("cpgStats -i cpgFile [-z] [-o results.json] [-s meth.values.json] [-r inforeads.values.json] [-b file.sorted.bed -a file.output.bed [-g] ] \n");
	printf("STATS:\n");
	printf("\t-i \t CpG Input file. \n");
	printf("\t-z \t CpG File is gzipped. \n");
	printf("\t-o \t JSON Output file. \n");
	printf("\t-s \t JSON Output Methylation Values File. \n");
	printf("\t-r \t JSON Output Information Reads Values File. \n");
	printf("ANNOTATION:\n");
	printf("\t-b \t SORTED BED Input file to annotate methylation per each window. \n");
	printf("\t-a \t Output annotated methylation file.\n");
	printf("INTERSECTION:\n");
    printf("\t-x \t CpG First File.\n");
	printf("\t-y \t CpG Second File.\n");
	printf("\t-g \t CpG Intersection Files are gzipped. \n");
	printf("VERSION:\n");
	printf("\t-v \t Program Version. \n");
}

/**
 * \brief Print program version
 */
void printVersion()
{
	printf("cpgStats version %s \n",VERSION);
}



/**
 * \brief Initialization of arguments
 */
void initArgs(struct Args * arguments)
{
	arguments->cpgInputFile = NULL;
	arguments->jsonFile = NULL;
	arguments->methJsonFile = NULL;
	arguments->bedFile = NULL;
	arguments->annotatedFile = NULL;
	arguments->isZipped = 0;
	arguments->firstCpGIsecFile = NULL;
	arguments->secondCpGIsecFile = NULL;
	arguments->areIsecZipped = 0;
}


/**
 * \brief Check Arguments
 * \params Args - struct of arguments
 * \return 1 if everything is ok  0 if not ok
 */
int checkArguments(struct Args * arguments)
{
	if (arguments->cpgInputFile == NULL)
	{
		printf("Sorry!! No CpG file input. Please specify it with -i parameter. \n");
		return 0;
	}
	return 1;
}

/**
 * \brief Get arguments
 * \params arguments - struct of arguments
 * \params argc - number of arguments
 * \params argv - arguments passed
 * \return 1 if everything is ok  0 if not ok
 */
int getArgs (struct Args * arguments, int argc, char *argv[])
{
	int opt = 0;

    initArgs(arguments);

	/*If not arguments Print help*/
	if (argc == 1)
	{
		printHelp();
		return 0;
	}

	while ((opt = getopt(argc, argv, "i:o:s:r:b:a:zhvx:y:g")) != -1)
	{
	    switch(opt)
	    {
	        case 'i':
	        	arguments->cpgInputFile = optarg;
	            break;
	        case 'o':
	            arguments->jsonFile = optarg;
	            break;
	        case 's':
	        	arguments->methJsonFile = optarg;
	        	break;
	        case 'r':
	        	arguments->infoReadsJsonFile = optarg;
	        	break;
	        case 'h':
	        	printHelp();
	        	return 0;
	        	break;
	        case 'b':
	        	arguments->bedFile = optarg;
	            break;
	        case 'a':
	            arguments->annotatedFile = optarg;
	            break;
	        case 'z':
	        	arguments->isZipped = 1;
	        	break;
	        case 'x':
	        	arguments->firstCpGIsecFile = optarg;
	        	break;
	        case 'y':
	        	arguments->secondCpGIsecFile = optarg;
	        	break;
	        case 'g':
	        	arguments->areIsecZipped = 1;
	        	break;
	        case 'v':
	        	printVersion();
	        	return 0;
	        	break;
	        case '?':
	            if (optopt == 'i')
	            {
	                printf("Missing mandatory input option \n");
	            }
	            else if (optopt == 'o')
	            {
	                printf("Missing mandatory output option \n");
	            }
	            else if (optopt == 'b')
	            {
	                printf("Missing bed file option \n");
	            }
	            else if (optopt == 'a')
	            {
	                printf("Missing annotation file option \n");
	            }
	            else if (optopt == 'x')
				{
					printf("Missing CpG Dinucleotides File One \n");
				}
				else if (optopt == 'y')
				{
					printf("Missing CpG Dinucleotides File Two \n");
				}
	            else
	            {
	                printf("Invalid option received \n");
	            }
	            break;
	    }
	}


	/*return checkArguments(arguments);*/
	return 1;
}


