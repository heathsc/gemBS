/*
 * counts.cpp
 *
 *  Created on: 20 Mai, 2016
 *      Author: marcos
 */

#include "counts.h"
#include <stdlib.h>
#include <string.h>

/**
 * \brief Counts Initialization
 */
void initCounts()
{
	unsigned int i;

	for (i=0; i < 2; i=i+1)
	{
	    counts.referenceCGsMethylated[i] = 0;
	    counts.referenceCGsInterMethylated[i] = 0;
	    counts.referenceCGsUnMethylated[i] = 0;

	    counts.nonReferenceCGsMethylated[i] = 0;
	    counts.nonReferenceCGsInterMethylated[i] = 0;
	    counts.nonReferenceCGsUnMethylated[i] = 0;
	}

	//Initialize vector methylation value
	vector_init(&vReferenceCGsMethValuesQ20);
	vector_init(&vNonReferenceCGsMethValuesQ20);

	//Initialize Vector of Information Reads
	info_reads_init(&vInfoReadsReferenceCGs);
	info_reads_init(&vInfoReadsNonReferenceCGs);
	info_reads_init(&vInfoReadsReferenceCGsQ20);
	info_reads_init(&vInfoReadsNonReferenceCGsQ20);
}

/**
 * \brief Add a new record on a given truct
 * \param Record stats
 * \param Counter vector to Update
 */
void addCounts(struct Record * record,unsigned int * cnt)
{
	if(record->phredScore < 20)
	{
		/*Quality Under 20*/
		cnt[0] = cnt[0] + 1;
	}
	else if(record->phredScore >= 20)
	{
		/*Quality over 20*/
		cnt[1] = cnt[1] + 1;
	}
}

/**
 * \brief Updates a given info reads register with a given number of information reads
 * \param VectorInfoReads Vector which stores Information Reads and number of CpGs
 * \param infoReads Number of Information reads for a given CpG
 */
void updateInfoReadsRegister(VectorInfoReads *vector,unsigned int infoReads)
{
	unsigned int infoReadsRecorded = 0;
    unsigned int i = 0;

    InfoReads infReads;

    for(i=0; i < vector->size;i++)
    {
    	infReads = info_reads_get(vector,i);

    	if(infReads.informationReads == infoReads)
		{
    		infoReadsRecorded = 1;
		    infReads.cpgs = infReads.cpgs + 1;
		    info_reads_set(vector,i,infReads);
		    return;
		}
    }

    if(infoReadsRecorded == 0)
    {
    	infReads.informationReads = infoReads;
		infReads.cpgs = 1;
		info_reads_append(vector,infReads);
	}
}

/**
 * \brief Add record to the global stats
 * \params record -  Record Stats to add
 */
void addRecordStats(struct Record * record)
{
	//TODO METH VALUES VECTOR https://www.happybearsoftware.com/implementing-a-dynamic-array


    /* 1. Check contig existance */
	if(record->contig == NULL)
	{
		return;
	}

	/* 2. Check Methylation Value existance */
	if(record->noValue == 1)
	{
		return;
	}


	if(strcmp(record->referenceContext, record->callContext) == 0)
	{
        /*3.1 Reference CGs Dinucleotide */
		if (record->methValue <= 0.3)
		{
			/*3.1.1 UnMethylated*/
			addCounts(record,counts.referenceCGsUnMethylated);
		}
		else if(record->methValue > 0.3 && record->methValue <= 0.7)
		{
			/*3.1.2 Intermediate Methylation*/
			addCounts(record,counts.referenceCGsInterMethylated);
		}
		else
		{
			/*3.1.3 Methylated*/
			addCounts(record,counts.referenceCGsMethylated);
		}
	}
	else
	{
		/*3.3 Non Reference CpGs Dinucleotide */
		if (record->methValue <= 0.3)
		{
			/*3.3.1 Un Methylated*/
			addCounts(record,counts.nonReferenceCGsUnMethylated);
		}
		else if(record->methValue > 0.3 && record->methValue <= 0.7)
		{
			/*3.3.2 Intermediate Methylation*/
			addCounts(record,counts.nonReferenceCGsInterMethylated);
		}
		else
		{
			/*3.3.3 Methylated*/
			addCounts(record,counts.nonReferenceCGsMethylated);
		}
	}


	/*4. Append Mehylation Levels*/
	if(record->phredScore >= 20)
	{
        if(strcmp(record->referenceContext,"CG") == 0 && strcmp(record->referenceContext,record->callContext) == 0)
	    {
	        vector_append(&vReferenceCGsMethValuesQ20, record->methValue);
		}
		else
		{
		   	vector_append(&vNonReferenceCGsMethValuesQ20, record->methValue);
		}
	}

    /*5. Append InfoReads Levels*/
	if(strcmp(record->referenceContext,"CG") == 0 && strcmp(record->referenceContext,record->callContext) == 0)
	{
		updateInfoReadsRegister(&vInfoReadsReferenceCGs,record->infoReads);
		if (record->phredScore >= 20)
		{
			updateInfoReadsRegister(&vInfoReadsReferenceCGsQ20,record->infoReads);
		}
	}
	else
	{
		updateInfoReadsRegister(&vInfoReadsNonReferenceCGs,record->infoReads);
		if (record->phredScore >= 20)
		{
			updateInfoReadsRegister(&vInfoReadsNonReferenceCGsQ20,record->infoReads);
		}
	}
}


/**
 * \brief Get Total from a given vector
 * \param Get counter vector
 * \return Return total value
 */
unsigned int getTotalCount(unsigned int * cnt)
{
	unsigned int i = 0;
	unsigned int total = 0;

	for (i=0; i< 2; i=i+1)
	{
		total = total + cnt[i];
	}
	return total;
}



/* Get Total Reference CGs Counts */
unsigned int getTotalReferenceCGs(){return (getTotalCount(counts.referenceCGsMethylated) + getTotalCount(counts.referenceCGsInterMethylated) + getTotalCount(counts.referenceCGsUnMethylated));}

/* Get Total Reference CGs Counts Quality 20*/
unsigned int getTotalReferenceCGsQuality20()
{
	return counts.referenceCGsMethylated[1] + counts.referenceCGsInterMethylated[1] + counts.referenceCGsUnMethylated[1];
}

/* Get Total Non Reference CGs Counts */
unsigned int getTotalNonReferenceCGs(){return (getTotalCount(counts.nonReferenceCGsMethylated) + getTotalCount(counts.nonReferenceCGsInterMethylated) + getTotalCount(counts.nonReferenceCGsUnMethylated));}

/* Get Total Non Reference CGs Quality 20*/
unsigned int getTotalNonReferenceCGsQuality20()
{
	return counts.nonReferenceCGsMethylated[1] + counts.nonReferenceCGsInterMethylated[1] + counts.nonReferenceCGsUnMethylated[1];
}

/* Get Total Methylated */
unsigned int getTotalMethylated(){ return (getTotalCount(counts.referenceCGsMethylated) + getTotalCount(counts.nonReferenceCGsMethylated));}

/* Get Total Methylated Quality 20*/
unsigned int getTotalMethylatedQuality20(){ return counts.referenceCGsMethylated[1] + counts.nonReferenceCGsMethylated[1];}


/* Get Total Intermediate Methylated */
unsigned int getTotalIntermediateMethylated(){ return (getTotalCount(counts.referenceCGsInterMethylated) + getTotalCount(counts.nonReferenceCGsInterMethylated));}

/* Get Total Intermediate Methylated Quality 20*/
unsigned int getTotalIntermediateMethylatedQuality20(){ return counts.referenceCGsInterMethylated[1] + counts.nonReferenceCGsInterMethylated[1];}

/* Get Total UnMethylated */
unsigned int getTotalUnMethylated(){ return (getTotalCount(counts.referenceCGsUnMethylated) + getTotalCount(counts.nonReferenceCGsUnMethylated));}

/* Get Total UnMethylated Quality 20*/
unsigned int getTotalUnMethylatedQuality20(){ return counts.referenceCGsUnMethylated[1] + counts.nonReferenceCGsUnMethylated[1];}

/* Get Total Methylated Reference CGs */
unsigned int getTotalMethylatedReferenceCGs(){ return getTotalCount(counts.referenceCGsMethylated);}

/* Get Total Methylated Homozigous Quality 20 */
unsigned int getTotalMethylatedReferenceCGsQuality20(){ return counts.referenceCGsMethylated[1];}

/* Get Total Methylated NonReference CGs */
unsigned int getTotalMethylatedNonReferenceCGs(){ return getTotalCount(counts.nonReferenceCGsMethylated);}

/* Get Total Methylated NonReference CGs Quality 20 */
unsigned int getTotalMethylatedNonReferenceCGsQuality20(){ return counts.nonReferenceCGsMethylated[1];}


/* Get Total Intermediate Methylated Reference CGs */
unsigned int getTotalIntermediateMethylatedReferenceCGs(){ return getTotalCount(counts.referenceCGsInterMethylated);}

/* Get Total Intermediate Methylated ReferenceCGs Quality 20 */
unsigned int getTotalIntermediateMethylatedReferenceCGsQuality20(){ return counts.referenceCGsInterMethylated[1];}

/* Get Total Intermediate Methylated NonReference CGs */
unsigned int getTotalIntermediateMethylatedNonReferenceCGs(){ return getTotalCount(counts.nonReferenceCGsInterMethylated);}

/* Get Total Intermediate Methylated NonReference CGs Quality 20 */
unsigned int getTotalIntermediateMethylatedNonReferenceCGsQuality20(){ return counts.nonReferenceCGsInterMethylated[1];}


/* Get Total UnMethylated ReferenceCGs */
unsigned int getTotalUnMethylatedReferenceCGs(){ return getTotalCount(counts.referenceCGsUnMethylated);}

/* Get Total UnMethylated ReferenceCGs Quality 20 */
unsigned int getTotalUnMethylatedReferenceCGsQuality20(){ return counts.referenceCGsUnMethylated[1];}

/* Get Total UnMethylated NonReference CGs */
unsigned int getTotalUnMethylatedNonReferenceCGs(){ return getTotalCount(counts.nonReferenceCGsUnMethylated);}

/* Get Total UnMethylated NonReference CGs Quality 20 */
unsigned int getTotalUnMethylatedNonReferenceCGsQuality20(){ return counts.nonReferenceCGsUnMethylated[1];}

/* Get Total number of Dinucleotides*/
unsigned int getTotalDinucleotides(){ return getTotalReferenceCGs() + getTotalNonReferenceCGs();}

/* Get Total Quality 20 */
unsigned int getTotalQuality20()
{
	return counts.referenceCGsMethylated[1] + counts.referenceCGsInterMethylated[1] + counts.referenceCGsUnMethylated[1] +
		   counts.nonReferenceCGsMethylated[1] + counts.nonReferenceCGsInterMethylated[1] + counts.nonReferenceCGsUnMethylated[1];
}


/**
 * \brief get percentage value
 * \param concept to calculate the percentage
 * \param total Value from which estimate the percentage
 */
float getPercentage(unsigned int concept, unsigned int totalValue)
{
	return ((float)concept/totalValue)*100;
}

/**
 * \brief Print Counters
 */
void printCounts()
{
	printf("CpG Stats \n");
	printf("Total Dinucleotides: \t %i \t (%.2f %%) \n", getTotalDinucleotides(),getPercentage(getTotalDinucleotides(), getTotalDinucleotides()));
	printf("\t Total Quality > 20: \t %i \t (%.2f %%) \n", getTotalQuality20(),getPercentage(getTotalQuality20(), getTotalDinucleotides()));

	printf("\n");
	printf("Total Reference CGs: \t %i \t (%.2f %%) \n",getTotalReferenceCGs(),getPercentage(getTotalReferenceCGs(), getTotalDinucleotides()));
	printf("\t Total Reference CGs Quality > 20: \t %i \t (%.2f %%) \n",getTotalReferenceCGsQuality20(),getPercentage(getTotalReferenceCGsQuality20(), getTotalDinucleotides()));
	printf("Total Non Reference CGs: \t %i \t (%.2f %%) \n",getTotalNonReferenceCGs(),getPercentage(getTotalNonReferenceCGs(), getTotalDinucleotides()));
	printf("\t Total Non Reference CGs Quality > 20: \t %i \t (%.2f %%) \n",getTotalNonReferenceCGsQuality20(),getPercentage(getTotalNonReferenceCGsQuality20(), getTotalDinucleotides()));

	printf("\n");
	printf("Total Methylated: \t %i \t (%.2f %%) \n", getTotalMethylated(),getPercentage(getTotalMethylated(), getTotalDinucleotides()));
	printf("\t Total Methylated Quality > 20: \t %i \t (%.2f %%) \n", getTotalMethylatedQuality20(),getPercentage(getTotalMethylatedQuality20(), getTotalDinucleotides()));
    printf("Total Intermediate Methylated: \t %i \t (%.2f %%) \n", getTotalIntermediateMethylated(),getPercentage(getTotalIntermediateMethylated(), getTotalDinucleotides()));
	printf("\t Total Intermediate Methylated Quality > 20: \t %i \t (%.2f %%) \n", getTotalIntermediateMethylatedQuality20(),getPercentage(getTotalIntermediateMethylatedQuality20(), getTotalDinucleotides()));
	printf("Total UnMethylated: \t %i \t (%.2f %%) \n", getTotalUnMethylated(),getPercentage(getTotalUnMethylated(), getTotalDinucleotides()));
	printf("\t Total Unmethylated Quality > 20: \t %i \t (%.2f %%) \n", getTotalUnMethylatedQuality20(),getPercentage(getTotalUnMethylatedQuality20(), getTotalDinucleotides()));

	printf("\n");
	printf("Total Reference CGs Methylated: \t %i \t (%.2f %%)  \n", getTotalMethylatedReferenceCGs(),getPercentage(getTotalMethylatedReferenceCGs(),getTotalMethylated()));
	printf("\t Total Reference CGs Methylated Quality > 20: \t %i \t (%.2f %%)  \n", getTotalMethylatedReferenceCGsQuality20(),getPercentage(getTotalMethylatedReferenceCGsQuality20(),getTotalMethylated()));
	printf("Total Non Reference CGs Methylated: \t %i \t (%.2f %%)  \n", getTotalMethylatedNonReferenceCGs(),getPercentage(getTotalMethylatedNonReferenceCGs(),getTotalMethylated()));
	printf("\t Total Non Reference CGs Methylated Quality > 20: \t %i \t (%.2f %%)  \n", getTotalMethylatedNonReferenceCGsQuality20(),getPercentage(getTotalMethylatedNonReferenceCGsQuality20(),getTotalMethylated()));

	printf("\n");
	printf("Total Reference Intermediate Methylated: \t %i \t (%.2f %%)  \n", getTotalIntermediateMethylatedReferenceCGs(),getPercentage(getTotalIntermediateMethylatedReferenceCGs(),getTotalIntermediateMethylated()));
	printf("\t Total Reference Intermediate Methylated Quality > 20: \t %i \t (%.2f %%)  \n", getTotalIntermediateMethylatedReferenceCGsQuality20(),getPercentage(getTotalIntermediateMethylatedReferenceCGsQuality20(),getTotalIntermediateMethylated()));
	printf("Total Non Reference Intermediate Methylated: \t %i \t (%.2f %%)  \n", getTotalIntermediateMethylatedNonReferenceCGs(),getPercentage(getTotalIntermediateMethylatedNonReferenceCGs(),getTotalIntermediateMethylated()));
	printf("\t Total Non Reference Intermediate Methylated Quality > 20: \t %i \t (%.2f %%)  \n", getTotalIntermediateMethylatedNonReferenceCGsQuality20(),getPercentage(getTotalIntermediateMethylatedNonReferenceCGsQuality20(),getTotalIntermediateMethylated()));

	printf("\n");
	printf("Total Reference CGs UnMethylated: \t %i \t (%.2f %%) \n", getTotalUnMethylatedReferenceCGs(),getPercentage(getTotalUnMethylatedReferenceCGs(),getTotalUnMethylated()));
	printf("\t Total Reference CGs UnMethylated Quality > 20: \t %i \t (%.2f %%) \n", getTotalUnMethylatedReferenceCGsQuality20(),getPercentage(getTotalUnMethylatedReferenceCGsQuality20(),getTotalUnMethylated()));
	printf("Total Non Reference CGs UnMethylated: \t %i \t (%.2f %%) \n", getTotalUnMethylatedNonReferenceCGs(),getPercentage(getTotalUnMethylatedNonReferenceCGs(),getTotalUnMethylated()));
	printf("\t Total Non Reference CGs UnMethylated Quality > 20: \t %i \t (%.2f %%) \n", getTotalUnMethylatedNonReferenceCGsQuality20(),getPercentage(getTotalUnMethylatedNonReferenceCGsQuality20(),getTotalUnMethylated()));
}

/**
 * \brief Print to JSON file
 * \param json file to store the data
 */
void saveCounts(char * fileName)
{
	FILE *fp;

	fp = fopen(fileName, "w");

	fprintf(fp,"{\n");

	fprintf(fp,"  \"TotalDinucleotides\":%i,\n", getTotalDinucleotides());
	fprintf(fp,"  \"TotalQuality20\":%i,\n", getTotalQuality20());

	fprintf(fp,"  \"TotalReferenceCGs\":%i,\n",getTotalReferenceCGs());
	fprintf(fp,"  \"TotalReferenceCGsQuality20\":%i,\n",getTotalReferenceCGsQuality20());

	fprintf(fp,"  \"TotalNonReferenceCGs\":%i,\n",getTotalNonReferenceCGs());
	fprintf(fp,"  \"TotalNonReferenceCGsQuality20\":%i,\n",getTotalNonReferenceCGsQuality20());

	fprintf(fp,"  \"TotalMethylated\":%i,\n", getTotalMethylated());
	fprintf(fp,"  \"TotalMethylatedQuality20\":%i,\n", getTotalMethylatedQuality20());
	fprintf(fp,"  \"TotalIntermediateMethylated\":%i,\n", getTotalIntermediateMethylated());
	fprintf(fp,"  \"TotalIntermediateMethylatedQuality20\":%i,\n", getTotalIntermediateMethylatedQuality20());
	fprintf(fp,"  \"TotalUnMethylated\":%i,\n", getTotalUnMethylated());
	fprintf(fp,"  \"TotalUnmethylatedQuality20\":%i,\n", getTotalUnMethylatedQuality20());

	fprintf(fp,"  \"TotalReferenceCGsMethylated\":%i,\n", getTotalMethylatedReferenceCGs());
	fprintf(fp,"  \"TotalReferenceCGsMethylatedQuality20\":%i,\n", getTotalMethylatedReferenceCGsQuality20());
	fprintf(fp,"  \"TotalNonReferenceCGsMethylated\":%i,\n", getTotalMethylatedNonReferenceCGs());
	fprintf(fp,"  \"TotalNonReferenceCGsMethylatedQuality20\":%i,\n", getTotalMethylatedNonReferenceCGsQuality20());

	fprintf(fp,"  \"TotalReferenceCGsIntermediateMethylated\":%i,\n", getTotalIntermediateMethylatedReferenceCGs());
	fprintf(fp,"  \"TotalReferenceCGsIntermediateMethylatedQuality20\":%i,\n", getTotalIntermediateMethylatedReferenceCGsQuality20());
	fprintf(fp,"  \"TotalNonReferenceCGsIntermediateMethylated\":%i,\n", getTotalIntermediateMethylatedNonReferenceCGs());
	fprintf(fp,"  \"TotalNonReferenceCGsIntermediateMethylatedQuality20\":%i,\n", getTotalIntermediateMethylatedNonReferenceCGsQuality20());

	fprintf(fp,"  \"TotalReferenceCGsUnMethylated\":%i,\n", getTotalUnMethylatedReferenceCGs());
	fprintf(fp,"  \"TotalReferenceCGsUnMethylatedQuality20\":%i,\n", getTotalUnMethylatedReferenceCGsQuality20());
	fprintf(fp,"  \"TotalNonReferenceCGsUnMethylated\":%i,\n", getTotalUnMethylatedNonReferenceCGs());
	fprintf(fp,"  \"TotalNonReferenceCGsUnMethylatedQuality20\":%i\n", getTotalUnMethylatedNonReferenceCGsQuality20());

	fprintf(fp,"}\n");

	fclose(fp);


}

/**
 * \brief Print Methylation values to JSON file
 * \param json file to store the methylation data
 */
void saveJsonMethylationCounts(char * fileName)
{
	FILE *fp;

	fp = fopen(fileName, "w");

	fprintf(fp,"{\n");

	unsigned int i = 0;

	/*Reference CGs Methylation Values Quality Over 20 */
	fprintf(fp,"  \"referenceCGsMethValuesQ20\":[");

	for(i=0; i < vReferenceCGsMethValuesQ20.size;i++)
	{
		if(i == (vReferenceCGsMethValuesQ20.size -1))
		{
			fprintf(fp," %.3f",vector_get(&vReferenceCGsMethValuesQ20,i));
		}
		else
		{
			fprintf(fp," %.3f,",vector_get(&vReferenceCGsMethValuesQ20,i));
		}
	}

	fprintf(fp,"],\n");

	/*Non Reference CGs Methylation Values Quality Over 20 */
	fprintf(fp,"  \"nonReferenceCGsMethValuesQ20\":[");

	for (i=0; i < vNonReferenceCGsMethValuesQ20.size;i++)
	{
		if(i == (vNonReferenceCGsMethValuesQ20.size -1))
		{
			fprintf(fp," %.3f",vector_get(&vNonReferenceCGsMethValuesQ20,i));
		}
		else
		{
			fprintf(fp," %.3f,",vector_get(&vNonReferenceCGsMethValuesQ20,i));
		}
	}

	fprintf(fp,"]\n");

	fprintf(fp,"}\n");

	fclose(fp);

	/*Close Vector*/
    vector_free(&vReferenceCGsMethValuesQ20);
    vector_free(&vNonReferenceCGsMethValuesQ20);
}

/**
 * Print Dictionary of Information Reads Vector
 *
 * fp -- File pointer descriptor
 * vector -- Information Reads Vector
 * text -- Text to print
 */
void printVectorInformationReads(FILE * fp,VectorInfoReads * vector,char * text)
{
	unsigned int i = 0;

	fprintf(fp,"  \"%s\":{",text);

	for (i=0; i < vector->size;i++)
	{
		InfoReads ir = info_reads_get(vector,i);

		if(i == (vector->size -1))
		{
			fprintf(fp,"   \"%i\":%i \n",ir.informationReads,ir.cpgs);
		}
		else
		{
			fprintf(fp,"   \"%i\":%i ,\n",ir.informationReads,ir.cpgs);
		}
	}
}


/**
 * \brief Print Information Reads values to JSON file
 * \param json file to store the methylation data
 */
void saveJsonInformationReads(char * fileName)
{
	FILE *fp;
    fp = fopen(fileName, "w");
    fprintf(fp,"{\n");

    printVectorInformationReads(fp,&vInfoReadsReferenceCGs,"informationReadsReferenceCGs");
    fprintf(fp," },\n");
    printVectorInformationReads(fp,&vInfoReadsNonReferenceCGs,"informationReadsNonReferenceCGs");
    fprintf(fp," },\n");
    printVectorInformationReads(fp,&vInfoReadsReferenceCGsQ20,"informationReadsReferenceCGsQ20");
    fprintf(fp," },\n");
    printVectorInformationReads(fp,&vInfoReadsNonReferenceCGsQ20,"informationReadsNonReferenceCGsQ20");
    fprintf(fp," }\n");

    fprintf(fp,"}\n");
	fclose(fp);

	/*Close Vector*/
	info_reads_free(&vInfoReadsReferenceCGs);
	info_reads_free(&vInfoReadsNonReferenceCGs);
	info_reads_free(&vInfoReadsReferenceCGsQ20);
	info_reads_free(&vInfoReadsNonReferenceCGsQ20);
}









