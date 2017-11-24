/*
 * parseInput.c
 *
 *  Created on: 20 Mai, 2016
 *      Author: marcos
 */

#include "parseInput.h"
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>

/**
 * \brief Record struct to represent the read parsed
 * \param Record structure to be initializated
 */
void initRecord(struct Record * record)
{
	record->contig = NULL;
	record->position = 0;
	record->referenceContext = NULL;
	record->callContext = NULL;
	record->phredScore = 0;
	record->methValue = 0.0;
	record->methDev = 0.0;
	record->noValue = 0;
	record->infoReads = 0;
}

/**
 * \brief From Read line to record read
 * \params line - line to transform
 * \params record - record to create from line
 */
void fromLineToRecord (char *line,struct Record * record)
{
	char * pch;
	pch = strtok (line," \t");

	unsigned int field = 0;

	while (pch != NULL)
	{
	    char * content = pch;

	    switch(field)
	    {
	        case 0:
	        	record->contig = content;
	        	break;
	        case 1:
	        	record->position = atoi(content);
	        	break;
	        case 2:
	        	record->referenceContext = content;
	        	break;
	        case 3:
	        	record->callContext = content;
	        	break;
	        case 4:
	        	record->phredScore = atoi(content);
	        	break;
	        case 5:
	        	record->methValue = atof(content);

	        	if(strcmp(content, "-") == 0)
	        	{
	        	    record->noValue = 1;
	        	}
	        	break;
	        case 6:
	        	record->methDev = atof(content);

	        	if(strcmp(content, "-") == 0)
	        	{
	        	    record->noValue = 1;
	        	}
	        	break;
	        case 7:
	        case 8:
	        	record->infoReads = record->infoReads + atoi(content);
	        	break;
	    }

	    field = field + 1;
	    pch = strtok (NULL, " \t");
	}
}

/********************************************************************************************************************/
/***************************************      GENERAL INPUT CALLS    ************************************************/
/********************************************************************************************************************/

/**
 * \brief Config Input file
 * \param inputName Input file name
 * \param is input zipped 1 yes 0 not
 * \returns 1 in success otherwise 0
 */
int configInput(char * inputName,int isZip)
{
	if(isZip)
	{
        return gzConfigInput(inputName);
	}
	else
	{
		return textConfigInput(inputName);
	}
	return 0;
}

/**
 * \brief Close File
 * \param is input zipped 1 yes 0 not
 */
void closeInput(int isZip)
{
	if(isZip)
	{
		return gzCloseInput();
	}
	else
	{
		return textCloseInput();
	}
}

/**
 * \brief Read a new line from file or input stream and returns a record
 * \param record - Record to be returned
 * \param is input zipped 1 yes 0 not
 * \returns 1 if there is more reads to read otherwise 0
 */
int getNewRecord(struct Record * record,int isZip)
{
	if(isZip)
	{
		return gzGetNewRecord(record);
	}
	else
	{
		return textGetNewRecord(record);
	}
}

/********************************************************************************************************************/
/***************************************     NOT COMPRESSED FILES    ************************************************/
/********************************************************************************************************************/

/**
 * \brief Config Input file
 * \param inputName Input file name
 * \returns 1 in success otherwise 0
 */
int textConfigInput(char * inputName)
{
	inputData = fopen(inputName, "r");

	if (inputData == NULL)
	{
		printf("Sorry!! Not possible to read file: %s \n",inputName);
		printf("Errno:%i \n",errno);
		exit(EXIT_FAILURE);
	}
    return 1;
}

/**
 * \brief Close File
 * \param inputName -  File input name
 */
void textCloseInput()
{
	fclose(inputData);
}

/**
 * \brief Read a new line from file or input stream and returns a record
 * \param record - Record to be returned
 * \returns 1 if there is more reads to read otherwise 0
 */
int textGetNewRecord(struct Record * record)
{
	char * line = NULL;
	size_t len = 0;
	ssize_t read;

	/*0. INIT record STRUCTURE */
	initRecord(record);

	/*1. GET read LINE */
	read = getline(&line, &len, inputData);

	if (read != -1)
	{
		fromLineToRecord(line,record);
		return 1;
	}

    return 0;
}

/********************************************************************************************************************/
/***************************************    GZIP COMPRESSED FILES    ************************************************/
/********************************************************************************************************************/

int gzConfigInput(char * inputName)
{
	gzInputFile = gzopen (inputName, "r");
	if (!inputName)
	{
		fprintf (stderr, "gzopen of '%s' failed: %s.\n", inputName,strerror (errno));
		exit (EXIT_FAILURE);
	}
	return 1;
}

void gzCloseInput()
{
	gzclose(gzInputFile);
}

/**
 * \brief Read a new line from Gzipped file and returns a record
 * \param record - Record to be returned
 * \returns 1 if there is more reads to read otherwise 0
 */
int gzGetNewRecord(struct Record * record)
{
	char * readGz;
	char buffer[LENGTH];

	/*0. INIT record STRUCTURE */
	initRecord(record);

	/*1. GET read LINE */
	readGz = gzgets (gzInputFile, buffer, LENGTH - 1);

	if (readGz != 0)
	{
		fromLineToRecord(buffer,record);
		return 1;
	}

	return 0;
}

/********************************************************************************************************************/
/****************************************        BED FILE INPUT    **************************************************/
/********************************************************************************************************************/

/**
 * \brief Removes carry return character
 * \param string to remove last carry return character
 */
void chomp(const char *s)
{
    char *p;
    while (NULL != s && NULL != (p = strrchr(s, '\n')))
    {
        *p = '\0';
    }
}


/**
 * \brief Record struct to represent the window bed parsed
 * \param Record structure to be initializated
 */
void initBed(struct Bed * window)
{
	window->contig = NULL;
	window->start = 0;
	window->end = 0;
	window->extra = NULL;
	window->meanMeth = 0.0;
	window->medianMeth = 0.0;
	window->stDevMeth = 0.0;
	window->cpgDinucleotides = 0;
	window->snps = 0;
}

/**
 * \brief Check Input Bed File
 * \param inputName Input file name
 * \returns 1 in success otherwise 0
 */
int checkInputBed(char * inputName)
{
	inputBed = fopen(inputName, "r");

	if (inputBed == NULL)
	{
		printf("Sorry!! Not possible to read file: %s \n",inputName);
		exit(EXIT_FAILURE);
	}

	return 1;
}

/**
 * \brief Close File
 */
void closeInputBed()
{
	fclose(inputBed);
}

/**
 * \brief From Read line to window bed
 * \params line - line to transform
 * \params Bed - window bed to create from line
 */
void fromLineToBed(char *line,struct Bed * window)
{
	char * pch;
	pch = strtok (line," \t");

	unsigned int field = 0;

	while (pch != NULL)
	{
	    char * content = pch;

	    switch(field)
	    {
	        case 0:
	        	window->contig = content;
	        	break;
	        case 1:
	        	window->start = atoi(content);
	        	break;
	        case 2:
	        	window->end = atoi(content);
	        	break;
	        default:
	        	chomp(content);

	        	if (field == 3)
	        	{
	        		window->extra = malloc (strlen (content) + 1);   // Space for length plus nul
	        		strcpy (window->extra,content);                        // Copy the characters
	        	}
	        	else
	        	{
	        		window->extra = realloc(window->extra,strlen(window->extra) + strlen("\t") + 1);
	        		strcat(window->extra, "\t");
	        		window->extra = realloc(window->extra,strlen(window->extra) + strlen(content) + 1);
	        		strcat(window->extra, content);
	        	}
	        	break;
	    }

	    field = field + 1;
	    pch = strtok (NULL, " \t");
	}
}

/**
 * \brief Read a new line from file and returns a record
 * \param window - window to be returned
 * \returns 1 if there is more reads to read otherwise 0
 */
int getNewWindow(struct Bed * window)
{
	char * line = NULL;
	size_t len = 0;
	ssize_t read;

	/*0. INIT window STRUCTURE */
	initBed(window);

	/*1. GET read LINE */
	read = getline(&line, &len, inputBed);

	if (read != -1)
	{
		fromLineToBed(line,window);
		return 1;
	}

    return 0;
}



/********************************************************************************************************************/
/**************************************        GENERAL INPUT FILES    ***********************************************/
/********************************************************************************************************************/


/**
 * \brief Config Input file
 * \param inputName Input file name
 * \param void * fileDescriptor File Descriptor to be open
 * \param is input zipped 1 yes 0 not
 * \returns 1 in success otherwise 0
 */
int setupInput(char * inputName,void ** fileDescriptor,int isZip)
{
	if(isZip == 0)
	{
		(*fileDescriptor) = fopen(inputName, "r");
	}
	else
	{
		(*fileDescriptor) = gzopen (inputName, "r");
	}

	if ((*fileDescriptor) == NULL)
	{
		printf("Sorry!! Not possible to read file: %s \n",inputName);
		printf("Errno:%i \n",errno);
		exit(EXIT_FAILURE);
	}
	return 1;
}

/**
 * \brief Close File
 * \param void * fileDescriptor File Descriptor to be open
 * \param is input zipped 1 yes 0 not
 */
void fileCloseInput(void * fileDescriptor,int isZip)
{
	if(isZip == 1)
	{
		gzclose(fileDescriptor);
	}
	else
	{
		fclose(fileDescriptor);
	}
}

/**
 * \brief Read a new line from file or input stream and returns a record
 * \param void * fileDescriptor File Descriptor to be open
 * \param record - Record to be returned
 * \param is input zipped 1 yes 0 not
 * \returns 1 if there is more reads to read otherwise 0
 */
int readNewRecord(void * fileDescriptor,struct Record * record,int isZip)
{
	/*0. INIT record STRUCTURE */
	initRecord(record);

	/*1. GET read LINE */
	if(isZip == 1)
	{
		char * readGz;
		char buffer[LENGTH];

		readGz = gzgets (fileDescriptor, buffer, LENGTH - 1);

		if (readGz != 0)
		{
			fromLineToRecord(buffer,record);
			return 1;
		}
	}
	else
	{
		char * line = NULL;
		size_t len = 0;
		ssize_t read;


		read = getline(&line, &len, fileDescriptor);

		if (read != -1)
		{
			fromLineToRecord(line,record);
			return 1;
		}
	}
    return 0;
}
