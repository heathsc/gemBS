/*
 * Intersection.c
 *
 *  Created on: 6 Xuñ, 2016
 *      Author: marcos
 */
#include "intersection.h"
#include "parseInput.h"
#include "counts.h"
#include <stdlib.h>
#include <string.h>


/**
 * \brief Initialization of intersection results
 */
void initIsec()
{
	resultsIsec.n_total_high_quality_one = 0;
	resultsIsec.n_total_low_quality_one = 0;
	resultsIsec.n_total_high_quality_two = 0;
	resultsIsec.n_total_low_quality_two = 0;

	resultsIsec.n_shared_high_quality_one = 0;
	resultsIsec.n_shared_low_quality_one = 0;
	resultsIsec.n_shared_high_quality_two = 0;
	resultsIsec.n_shared_low_quality_two = 0;

	resultsIsec.n_private_high_quality_one = 0;
	resultsIsec.n_private_low_quality_one = 0;

	resultsIsec.n_private_high_quality_two = 0;
	resultsIsec.n_private_low_quality_two = 0;
}

/**
 * \brief Initialization of intersection Stack
 */
void initStack()
{
	unsigned int i = 0;
	for (i=0; i < STACK_SIZE; i++)
	{
		stackInfo[i].hasValue = 0;
		if (stackInfo[i].record != NULL)
		{
			free(stackInfo[i].record);
			stackInfo[i].record = NULL;
		}
		/*initRecord(&stackInfo[i].record);*/
	}
}

float testMe(unsigned int concept, unsigned int totalValue)
{
	return ((float)concept/totalValue)*100;
}

/**
 * \brief Print Intersection Results
 * \param nameFileOne File Name first file
 * \param nameFileTwo File Name second file
 */
void printIsecResults(char * nameFileOne,char * nameFileTwo)
{
    unsigned int nTotalDinucleotides = resultsIsec.n_total_high_quality_one + resultsIsec.n_total_low_quality_one + resultsIsec.n_total_high_quality_two +	resultsIsec.n_total_low_quality_two;
    printf("Totals \n");
    printf("\tTotal Dinucleotides Processed:\t%i\t(100 %%)\n",nTotalDinucleotides);

    unsigned int nTotalOne = resultsIsec.n_total_high_quality_one + resultsIsec.n_total_low_quality_one;
    printf("\tTotal Dinucleotides file %s:\t%i\t(%.2f %%) \n",nameFileOne,nTotalOne,getPercentage(nTotalOne,nTotalDinucleotides));
    printf("\t\tTotal Dinucleotides High Quality file %s:\t%i\t(%.2f %%)\n",nameFileOne,resultsIsec.n_total_high_quality_one,getPercentage(resultsIsec.n_total_high_quality_one,nTotalOne));
    printf("\t\tTotal Dinucleotides Low Quality file %s:\t%i\t(%.2f %%)\n",nameFileOne,resultsIsec.n_total_low_quality_one,getPercentage(resultsIsec.n_total_low_quality_one,nTotalOne));

    unsigned int nTotalTwo = resultsIsec.n_total_high_quality_two + resultsIsec.n_total_low_quality_two;
    printf("\tTotal Dinucleotides file %s:\t%i\t(%.2f %%) \n",nameFileTwo,nTotalTwo,getPercentage(nTotalTwo,nTotalDinucleotides));
    printf("\t\tTotal Dinucleotides High Quality file %s:\t%i\t(%.2f %%)\n",nameFileTwo,resultsIsec.n_total_high_quality_two,getPercentage(resultsIsec.n_total_high_quality_two,nTotalTwo));
    printf("\t\tTotal Dinucleotides Low Quality file %s:\t%i\t(%.2f %%)\n",nameFileTwo,resultsIsec.n_total_low_quality_two,getPercentage(resultsIsec.n_total_low_quality_two,nTotalTwo));

    unsigned int nTotalShared = resultsIsec.n_shared_high_quality_one + resultsIsec.n_shared_low_quality_one + resultsIsec.n_shared_high_quality_two + resultsIsec.n_shared_low_quality_two;
    printf("Shared between %s and %s:\t%i\t(%.2f %%) \n ",nameFileOne,nameFileTwo,nTotalShared,getPercentage(nTotalShared,nTotalDinucleotides));
    unsigned int nSharedOne = resultsIsec.n_shared_high_quality_one + resultsIsec.n_shared_low_quality_one;
    printf("\tTotal Shared file %s:\t%i\t(%.2f %%) \n",nameFileOne,nSharedOne,getPercentage(nSharedOne,nTotalOne));
    printf("\t\tTotal Shared High Quality file %s:\t%i\t(%.2f %%)\n",nameFileOne,resultsIsec.n_shared_high_quality_one,getPercentage(resultsIsec.n_shared_high_quality_one,nSharedOne));
    printf("\t\tTotal Shared Low Quality file %s:\t%i\t(%.2f %%)\n",nameFileOne,resultsIsec.n_shared_low_quality_one,getPercentage(resultsIsec.n_shared_low_quality_one,nSharedOne));

    unsigned int nSharedTwo = resultsIsec.n_shared_high_quality_two + resultsIsec.n_shared_low_quality_two;
    printf("\tTotal Shared file %s:\t%i\t(%.2f %%) \n",nameFileTwo,nSharedTwo,getPercentage(nSharedTwo,nTotalTwo));
    printf("\t\tTotal Shared High Quality file %s:\t%i\t(%.2f %%)\n",nameFileTwo,resultsIsec.n_shared_high_quality_two,getPercentage(resultsIsec.n_shared_high_quality_two,nSharedTwo));
    printf("\t\tTotal Shared Low Quality file %s:\t%i\t(%.2f %%)\n",nameFileTwo,resultsIsec.n_shared_low_quality_two,getPercentage(resultsIsec.n_shared_low_quality_two,nSharedTwo));

    unsigned int nPrivateOne = resultsIsec.n_private_high_quality_one + resultsIsec.n_private_low_quality_one;
    printf("Private to %s \n",nameFileOne);
    printf("\tTotal Private %s:\t%i\t(%.2f %%) \n",nameFileOne,nPrivateOne,getPercentage(nPrivateOne,nTotalOne));
    printf("\t\tTotal Private High Quality %s:\t%i\t(%.2f %%) \n",nameFileOne,resultsIsec.n_private_high_quality_one,getPercentage(resultsIsec.n_private_high_quality_one,nPrivateOne));
    printf("\t\tTotal Private Low Quality %s:\t%i\t(%.2f %%) \n",nameFileOne,resultsIsec.n_private_low_quality_one,getPercentage(resultsIsec.n_private_low_quality_one,nPrivateOne));

    unsigned int nPrivateTwo = resultsIsec.n_private_high_quality_two + resultsIsec.n_private_low_quality_two;
    printf("Private to %s \n",nameFileTwo);
    printf("\tTotal Private %s:\t%i\t(%.2f %%) \n",nameFileTwo,nPrivateTwo,getPercentage(nPrivateTwo,nTotalTwo));
    printf("\t\tTotal Private High Quality %s:\t%i\t(%.2f %%) \n",nameFileTwo,resultsIsec.n_private_high_quality_two,getPercentage(resultsIsec.n_private_high_quality_two,nPrivateTwo));
    printf("\t\tTotal Private Low Quality %s:\t%i\t(%.2f %%) \n",nameFileTwo,resultsIsec.n_private_low_quality_two,getPercentage(resultsIsec.n_private_low_quality_two,nPrivateTwo));
}

/**
 * \brief Load file STACK_SIZE elemenets from one of the reads
 * \param void fileDescriptor File Descriptor to process
 * \param elements Number of elements in the vector
 * \param isGzip 1 if file is zipped otherwise 0
 * \return no more registers to read
 */
int loadStack(void * fileDescriptor,int isGzip,int * elements)
{
	unsigned int nRecord;

	initStack();
	(*elements) = 0;

	for (nRecord = 0; nRecord < STACK_SIZE;nRecord = nRecord+1)
	{
		struct Record record;
		if ( readNewRecord(fileDescriptor,&record,isGzip) != 0)
		{
			stackInfo[nRecord].hasValue = 1;

			stackInfo[nRecord].record = malloc(sizeof(record));

			memcpy(stackInfo[nRecord].record, &record, sizeof(record));

			(*elements) =  (*elements) +1;
		}
		else
		{
			return 0;
		}
	}

	return 1;
}

/**
 * \brief Add counts according to Phreed Score value
 * \param onePhredScore Phred Score Read File one
 * \param twoPhredScore Phred Score Read File two
 * \param isEqual One if One and two are equal 0 Different
 */
void addBothCounts(unsigned int onePhredScore,unsigned int twoPhredScore, unsigned int isEqual)
{
	if(isEqual)
	{
		if(onePhredScore < 30)
		{
			/*1.1.1 RECORD ONE LOW QUALITY*/
			resultsIsec.n_shared_low_quality_one = resultsIsec.n_shared_low_quality_one + 1;
			resultsIsec.n_total_low_quality_one = resultsIsec.n_total_low_quality_one + 1;
		}
		else
		{
			/*1.1.2 RECORD ONE HIGH QUALITY*/
			resultsIsec.n_shared_high_quality_one = resultsIsec.n_shared_high_quality_one + 1;
			resultsIsec.n_total_high_quality_one = resultsIsec.n_total_high_quality_one + 1;
		}

		if(twoPhredScore <30)
		{
			/*1.1.3 RECORD TWO LOW QUALITY*/
			resultsIsec.n_shared_low_quality_two = resultsIsec.n_shared_low_quality_two + 1;
			resultsIsec.n_total_low_quality_two = resultsIsec.n_total_low_quality_two +1;
		}
		else
		{
			/*1.1.4 RECORD TWO HIGH QUALITY*/
			resultsIsec.n_shared_high_quality_two = resultsIsec.n_shared_high_quality_two + 1;
			resultsIsec.n_total_high_quality_two = resultsIsec.n_total_high_quality_two + 1;
		}
	}
	else
	{
		if(onePhredScore < 30)
		{
			/*1.1.1 RECORD ONE LOW QUALITY*/
			resultsIsec.n_private_low_quality_one = resultsIsec.n_private_low_quality_one + 1;
			resultsIsec.n_total_low_quality_one = resultsIsec.n_total_low_quality_one + 1;
		}
		else
		{
			/*1.1.2 RECORD ONE HIGH QUALITY*/
			resultsIsec.n_private_high_quality_one = resultsIsec.n_private_high_quality_one + 1;
			resultsIsec.n_total_high_quality_one = resultsIsec.n_total_high_quality_one + 1;
		}

		if(twoPhredScore <30)
		{
			/*1.1.3 RECORD TWO LOW QUALITY*/
			resultsIsec.n_private_low_quality_two = resultsIsec.n_private_low_quality_two + 1;
			resultsIsec.n_total_low_quality_two = resultsIsec.n_total_low_quality_two +1;
		}
		else
		{
			/*1.1.4 RECORD TWO HIGH QUALITY*/
			resultsIsec.n_private_high_quality_two = resultsIsec.n_private_high_quality_two + 1;
			resultsIsec.n_total_high_quality_two = resultsIsec.n_total_high_quality_two + 1;
		}
	}
}


/**
 * \brief Add counts to one side according to Phreed Score value
 * \param phredScore Phred Score Of the site to modify counts
 * \param site to modify 1 or 2
 */
void addOneSiteCounts(unsigned int phredScore, unsigned int site)
{
	if(site == 1)
	{
		if(phredScore < 30)
		{
			/*1.1.1 RECORD ONE LOW QUALITY*/
			resultsIsec.n_private_low_quality_one = resultsIsec.n_private_low_quality_one + 1;
			resultsIsec.n_total_low_quality_one = resultsIsec.n_total_low_quality_one + 1;
		}
		else
		{
			/*1.1.2 RECORD ONE HIGH QUALITY*/
			resultsIsec.n_private_high_quality_one = resultsIsec.n_private_high_quality_one + 1;
			resultsIsec.n_total_high_quality_one = resultsIsec.n_total_high_quality_one + 1;
		}


	}
	else if (site == 2)
	{
		if(phredScore <30)
		{
			/*1.1.3 RECORD TWO LOW QUALITY*/
			resultsIsec.n_private_low_quality_two = resultsIsec.n_private_low_quality_two + 1;
			resultsIsec.n_total_low_quality_two = resultsIsec.n_total_low_quality_two +1;
		}
		else
		{
			/*1.1.4 RECORD TWO HIGH QUALITY*/
			resultsIsec.n_private_high_quality_two = resultsIsec.n_private_high_quality_two + 1;
			resultsIsec.n_total_high_quality_two = resultsIsec.n_total_high_quality_two + 1;
		}
	}
}


/**
 * \brief Compare records from file 1 and file 2
 * \param recordTwo Record from file to to compare with the stack of reads from file one
 * \param start_position_stack First position to be checked at the stack
 * \param processed_read_two Processed second read 0 no processed 1 processed
 */
void compareRecords(struct Record * recordTwo,int * start_positon_stack,int * processed_read_two)
{
	unsigned int nRecord;

	for(nRecord = (*start_positon_stack); nRecord < STACK_SIZE; nRecord = nRecord + 1)
	{
		struct StackPosition currentOne = stackInfo[nRecord];

		if (currentOne.hasValue)
		{
			/*1. POSITONS ARE EQUAL*/
			if((strcmp(recordTwo->contig,currentOne.record->contig)==0) && (recordTwo->position == currentOne.record->position))
			{
				/*1.1 CALLS ARE EQUAL */
				if( (strcmp(recordTwo->referenceContext,currentOne.record->referenceContext) == 0) && (strcmp(recordTwo->callContext,currentOne.record->callContext) == 0) )
				{
					addBothCounts(currentOne.record->phredScore,recordTwo->phredScore, 1);
					(*processed_read_two) = 1;
				}
				else
				{
                /*1.2 CALLS ARE DIFFERENT*/
					addBothCounts(currentOne.record->phredScore,recordTwo->phredScore, 0);
					(*processed_read_two) = 1;
				}
		     	(*start_positon_stack) = (*start_positon_stack) + 1;
		     	break;
			}
			else if( ( (strcmp(recordTwo->contig,currentOne.record->contig)==0) && (recordTwo->position > currentOne.record->position) ) ||
					 (strcmp(recordTwo->contig,currentOne.record->contig) > 0) )
			{
			/*2. READ TWO GREATER THAN READ ONE*/
				addOneSiteCounts(currentOne.record->phredScore, 1);
				(*start_positon_stack) = (*start_positon_stack) + 1;
			}
			else if( ( (strcmp(recordTwo->contig,currentOne.record->contig)==0) && (recordTwo->position < currentOne.record->position)) ||
					 (strcmp(recordTwo->contig,currentOne.record->contig) < 0) )
			{
			/*3. READ TWO SMALLER THAN READ TWO*/
				addOneSiteCounts(recordTwo->phredScore, 2);
				(*processed_read_two) = 1;
				break;
			}
		}
	}
}

/**
 * \brief Process Remaining Stack Reads from File 1
 * \param start_position_stack first postion to check remaing reads in stack
 */
void processRemainingStack(unsigned int * start_positon_stack)
{
	unsigned int i;

	for (i=(*start_positon_stack);i < STACK_SIZE;i=i+1)
	{
		struct StackPosition currentOne = stackInfo[i];

		if (currentOne.hasValue)
		{
			addOneSiteCounts(currentOne.record->phredScore, 1);
			(*start_positon_stack) = (*start_positon_stack) + 1;
		}
	}
}

/**
 * \brief Process remaining reads from file 1
 * \param void * fileDescriptor File Descriptor from Read 1 file
 * \param isGzip 1 if file is zipped otherwise 0
 */
void processRemainingReadsFileOne(void * fileDescriptor,int isGzip)
{
	struct Record recordOne;

	while (readNewRecord(fileDescriptor,&recordOne,isGzip) != 0)
	{
		addOneSiteCounts(recordOne.phredScore, 1);
	}
}

/**
 * \brief Process remaining reads from file 2
 * \param int processed_read_two 1 if current read two was already processed
 * \param fileDescriptor to get remaining reads
 * \param isGzip 1 if file is zipped otherwise 0
 * \param recordTwo last record read from file two
 */
void processRemainingReadsFileTwo(int processed_read_two,void * fileDescriptor,int isGzip,struct Record * recordTwo)
{
	if(processed_read_two == 0)
	{
		addOneSiteCounts(recordTwo->phredScore, 2);
	}

	struct Record recordRemainTwo;

	while (readNewRecord(fileDescriptor,&recordRemainTwo,isGzip) != 0)
	{
		addOneSiteCounts(recordRemainTwo.phredScore, 2);
	}
}

/**
 * \brief Run Intersection between two dinucleotide files
 * \param nameFileOne Name first Dinucleotide file
 * \param nameFileTwo Name second Dinucleotide file
 * \param isGzip 1 if input is zipped otherwise 0
 * \return 1 if everything goes well otherwise 0
 */
int runIsec(char * nameFileOne, char * nameFileTwo,int isGzip)
{
	/*1. INIT INPUT FILES*/
	void * fileDescriptor_one;
	void * fileDescriptor_two;
	unsigned int stack_elements;
	unsigned int hasToReadOne;
	unsigned int start_positon_stack;
	unsigned int processed_read_two;
	unsigned int current_elements;

	if(isGzip)
	{
		gzFile * zDescriptor_one;
		fileDescriptor_one = zDescriptor_one;
		gzFile * zDescriptor_two;
		fileDescriptor_two = zDescriptor_two;
	}
	else
	{
		FILE * textDescriptor_one;
		fileDescriptor_one = textDescriptor_one;
		FILE * textDescriptor_two;
		fileDescriptor_two = textDescriptor_two;
	}

	if (setupInput(nameFileOne,&fileDescriptor_one,isGzip) < 1)
	{
    	return 0;
	}

	if (setupInput(nameFileTwo,&fileDescriptor_two,isGzip) < 1)
	{
		return 0;
	}

	/*2. PARSE BOTH FILES*/
	start_positon_stack = 0;

	hasToReadOne = 0;
	hasToReadOne = loadStack(fileDescriptor_one,isGzip,&stack_elements);

	processed_read_two = -1;
	struct Record recordTwo;
	int checkReadFileTwo = 0;

	/*3. COMPARE FILE TO FILE */
	while (stack_elements)
	{
		if (processed_read_two != 0)
		{
			/*If previous record from second file  was already computed or is the first to be computed then get read*/
			checkReadFileTwo = readNewRecord(fileDescriptor_two,&recordTwo,isGzip);
		}
		else
		{
			/*If previous record was not computed (processed_read_two == 0) then keep on check the read on the Stack*/
			checkReadFileTwo = 1;
		}

		/*3.1. CHECK READ TWO IN THE STACK OF READS COMMING FROM READ ONE*/
		if (checkReadFileTwo != 0)
		{
			processed_read_two = 0;
			compareRecords(&recordTwo,&start_positon_stack,&processed_read_two);
			current_elements = stack_elements - start_positon_stack;

			if (current_elements == 0)
			{
				stack_elements = 0;
				/*3.2. STACK COMPLETED PROCESSED*/
				if (hasToReadOne)
				{
					/* 3.2.1 Get Stack of Reads */
					start_positon_stack = 0;
					hasToReadOne = loadStack(fileDescriptor_one,isGzip,&stack_elements);
					if(stack_elements == 0)
					{
						processRemainingReadsFileTwo(processed_read_two,fileDescriptor_two,isGzip,&recordTwo);
					}
				}
				else
				{
				    /*3.3 NO MORE READS FROM FILE ONE TO BE PROCESSED*/
					processRemainingReadsFileTwo(processed_read_two,fileDescriptor_two,isGzip,&recordTwo);
				}
			}
		}
		else
		{
			/*3.2.1 Process remaining reads in Stack*/
			processRemainingStack(&start_positon_stack);
            /*3.2.1 Process remaining reads from file 1*/
			processRemainingReadsFileOne(fileDescriptor_one,isGzip);
			stack_elements = 0;
		}
	}

	/*4. CLOSE FILES*/
	fileCloseInput(fileDescriptor_one,isGzip);
	fileCloseInput(fileDescriptor_two,isGzip);

	/*5. PRINT STATS RESULTS*/
	printIsecResults(nameFileOne,nameFileTwo);

	return 1;
}
