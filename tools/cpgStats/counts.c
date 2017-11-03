/*
 * counts.cpp
 *
 *  Created on: 20 Mai, 2016
 *      Author: marcos
 */

#include "counts.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/**
 * \brief Counts Initialization
 */
void initCounts()
{
	unsigned int i;

	for (i=0; i < 2; i=i+1)
	{
	    counts.homozygousMethylated[i] = 0;
	    counts.homozygousInterMethylated[i] = 0;
	    counts.homozygousUnMethylated[i] = 0;

	    counts.alternativeCXMethylated[i] = 0;
	    counts.alternativeCXInterMethylated[i] = 0;
	    counts.alternativeCXUnMethylated[i] = 0;

	    counts.nonReferenceCpgMethylated[i] = 0;
	    counts.nonReferenceCpgInterMethylated[i] = 0;
	    counts.nonReferenceCpgUnMethylated[i] = 0;

	    counts.snpsReferenceCpGsMethylated[i] = 0;
	    counts.snpsReferenceCpGsInterMethylated[i] = 0;
	    counts.snpsReferenceCpGsUnMethylated[i] = 0;
	}

	//Initialize vector methylation value
	vector_init(&vectorMethValues);
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
        /*3.1 Homozygous Dinucleotide */
		if (record->methValue <= 0.3)
		{
			/*3.1.1 UnMethylated*/
			addCounts(record,counts.homozygousUnMethylated);
		}
		else if(record->methValue > 0.3 && record->methValue <= 0.7)
		{
			/*3.1.2 Intermediate Methylation*/
			addCounts(record,counts.homozygousInterMethylated);
		}
		else
		{
			/*3.1.3 Methylated*/
			addCounts(record,counts.homozygousMethylated);
		}

		/*3.2 Append Mehylation Values for those CGs Homozygous and High Quality*/
		if(strcmp(record->referenceContext,"CG") == 0)
		{
			if(record->phredScore > 20)
			{
				vector_append(&vectorMethValues, record->methValue);
			}
		}
	}
	else
	{
		/*3.3 AlternativeCX Dinucleotide */
		if (record->methValue <= 0.3)
		{
			/*3.3.1 Un Methylated*/
			addCounts(record,counts.alternativeCXUnMethylated);
		}
		else if(record->methValue > 0.3 && record->methValue <= 0.7)
		{
			/*3.3.2 Intermediate Methylation*/
			addCounts(record,counts.alternativeCXInterMethylated);
		}
		else
		{
			/*3.3.3 Methylated*/
			addCounts(record,counts.alternativeCXMethylated);
		}

		/*3.4 NonReference CpGs Detectected*/
		if(strcmp(record->referenceContext, "CG") != 0 && strcmp(record->callContext, "CG") == 0 )
		{
			if (record->methValue <= 0.3)
			{
				/*3.4.1 Methylated*/
				addCounts(record,counts.nonReferenceCpgUnMethylated);
			}
			else if(record->methValue > 0.3 && record->methValue <= 0.7)
			{
				/*3.4.2 Intermediate Methylation*/
				addCounts(record,counts.nonReferenceCpgInterMethylated);
			}
			else
			{
				/*3.4.3 UnMethylated*/
				addCounts(record,counts.nonReferenceCpgMethylated);
			}
		}

		/*3.5 CpG With SNPs*/
		if(strcmp(record->referenceContext, "CG") == 0 && strcmp(record->callContext, "CG") != 0 )
		{
			if (record->methValue <= 0.3)
			{
				/*3.5.1 Methylated*/
				addCounts(record,counts.snpsReferenceCpGsUnMethylated);
			}
			else if(record->methValue > 0.3 && record->methValue <= 0.7)
			{
				/*3.5.2 Intermediate Methylation*/
				addCounts(record,counts.snpsReferenceCpGsInterMethylated);
			}
			else
			{
				/*3.5.3 UnMethylated*/
				addCounts(record,counts.snpsReferenceCpGsMethylated);
			}
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



/* Get Total Homozygous Counts */
unsigned int getTotalHomozygous(){return (getTotalCount(counts.homozygousMethylated) + getTotalCount(counts.homozygousInterMethylated) + getTotalCount(counts.homozygousUnMethylated));}

/* Get Total Homozygous Counts Quality 20*/
unsigned int getTotalHomozygousQuality20()
{
	return counts.homozygousMethylated[1] + counts.homozygousInterMethylated[1] + counts.homozygousUnMethylated[1];
}

/* Get Total AlternativeCX Counts */
unsigned int getTotalAlternativeCX(){return (getTotalCount(counts.alternativeCXMethylated) + getTotalCount(counts.alternativeCXInterMethylated) + getTotalCount(counts.alternativeCXUnMethylated));}

/* Get Total AlternativeCX Counts Quality 20*/
unsigned int getTotalAlternativeCXQuality20()
{
	return counts.alternativeCXMethylated[1] + counts.alternativeCXInterMethylated[1] + counts.alternativeCXUnMethylated[1];
}

/* Get Total Methylated */
unsigned int getTotalMethylated(){ return (getTotalCount(counts.alternativeCXMethylated) + getTotalCount(counts.homozygousMethylated));}

/* Get Total Methylated Quality 20*/
unsigned int getTotalMethylatedQuality20(){ return counts.alternativeCXMethylated[1] + counts.homozygousMethylated[1];}


/* Get Total Intermediate Methylated */
unsigned int getTotalIntermediateMethylated(){ return (getTotalCount(counts.alternativeCXInterMethylated) + getTotalCount(counts.homozygousInterMethylated));}

/* Get Total Intermediate Methylated Quality 20*/
unsigned int getTotalIntermediateMethylatedQuality20(){ return counts.alternativeCXInterMethylated[1] + counts.homozygousInterMethylated[1];}

/* Get Total UnMethylated */
unsigned int getTotalUnMethylated(){ return (getTotalCount(counts.homozygousUnMethylated) + getTotalCount(counts.alternativeCXUnMethylated));}

/* Get Total UnMethylated Quality 20*/
unsigned int getTotalUnMethylatedQuality20(){ return counts.homozygousUnMethylated[1] + counts.alternativeCXUnMethylated[1];}

/* Get Total Methylated Homozygous */
unsigned int getTotalMethylatedHomozygous(){ return getTotalCount(counts.homozygousMethylated);}

/* Get Total Methylated Homozigous Quality 20 */
unsigned int getTotalMethylatedHomozygousQuality20(){ return counts.homozygousMethylated[1];}

/* Get Total Methylated AlternativeCX */
unsigned int getTotalMethylatedAlternativeCX(){ return getTotalCount(counts.alternativeCXMethylated);}

/* Get Total Methylated AlternativeCX Quality 20 */
unsigned int getTotalMethylatedAlternativeCXQuality20(){ return counts.alternativeCXMethylated[1];}


/* Get Total Intermediate Methylated Homozygous */
unsigned int getTotalIntermediateMethylatedHomozygous(){ return getTotalCount(counts.homozygousInterMethylated);}

/* Get Total Intermediate Methylated Homozigous Quality 20 */
unsigned int getTotalIntermediateMethylatedHomozygousQuality20(){ return counts.homozygousInterMethylated[1];}

/* Get Total Intermediate Methylated AlternativeCX */
unsigned int getTotalIntermediateMethylatedAlternativeCX(){ return getTotalCount(counts.alternativeCXInterMethylated);}

/* Get Total Intermediate Methylated AlternativeCX Quality 20 */
unsigned int getTotalIntermediateMethylatedAlternativeCXQuality20(){ return counts.alternativeCXInterMethylated[1];}


/* Get Total UnMethylated Homozygous */
unsigned int getTotalUnMethylatedHomozygous(){ return getTotalCount(counts.homozygousUnMethylated);}

/* Get Total UnMethylated Homozigous Quality 20 */
unsigned int getTotalUnMethylatedHomozygousQuality20(){ return counts.homozygousUnMethylated[1];}

/* Get Total UnMethylated AlternativeCX */
unsigned int getTotalUnMethylatedAlternativeCX(){ return getTotalCount(counts.alternativeCXUnMethylated);}

/* Get Total UnMethylated AlternativeCX Quality 20 */
unsigned int getTotalUnMethylatedAlternativeCXQuality20(){ return counts.alternativeCXUnMethylated[1];}

/* Get Total number of Dinucleotides*/
unsigned int getTotalDinucleotides(){ return getTotalHomozygous() + getTotalAlternativeCX();}

/* Get Total Quality 20 */
unsigned int getTotalQuality20()
{
	return counts.homozygousMethylated[1] + counts.homozygousInterMethylated[1] + counts.homozygousUnMethylated[1] +
		   counts.alternativeCXMethylated[1] + counts.alternativeCXInterMethylated[1] + counts.alternativeCXUnMethylated[1];
}


/**
 * \brief Get Total Non Reference CpGs
 */
unsigned int getTotalNonReferenceCpgs()
{
	return getTotalCount(counts.nonReferenceCpgMethylated) + getTotalCount(counts.nonReferenceCpgInterMethylated) + getTotalCount(counts.nonReferenceCpgUnMethylated);
}

/**
 * \brief Get Total Non Reference CpG quality over 20
 */
unsigned int getTotalNonReferenceCpgsQuality20()
{
	return counts.nonReferenceCpgMethylated[1] + counts.nonReferenceCpgInterMethylated[1] + counts.nonReferenceCpgUnMethylated[1];
}



/**
 * \brief Total Non Reference CpG Methylated
 */
unsigned int getTotalNonReferenceCpgsMethylated(){ return getTotalCount(counts.nonReferenceCpgMethylated); }


/**
 * \brief Total Non Reference Quality 20 Methylated
 */
unsigned int getTotalNonReferenceCpgsMethylatedQuality20(){	return counts.nonReferenceCpgMethylated[1];}


/**
 * \brief Total Non Reference CpG Intermediate Methylated
 */
unsigned int getTotalNonReferenceCpgsIntermediateMethylated(){ return getTotalCount(counts.nonReferenceCpgInterMethylated);}

/**
 * \brief Total Non Reference Quality 20 Intermediate Methylated
 */
unsigned int getTotalNonReferenceCpgsIntermediateMethylatedQuality20(){ return counts.nonReferenceCpgInterMethylated[1];}


/**
 * \brief Total Non Reference CpGs UnMethylated
 */
unsigned int getTotalNonReferenceCpgsUnMethylated(){ return getTotalCount(counts.nonReferenceCpgUnMethylated);}


/**
 * \brief Total Non Reference Quality 20 UnMethylated
 */
unsigned int getTotalNonReferenceCpgsUnMethylatedQuality20(){ return counts.nonReferenceCpgUnMethylated[1]; }


/**
 * \brief Total CpG with SNP Called
 */
unsigned int getTotalSnpsReferenceCpGs()
{
	return getTotalCount(counts.snpsReferenceCpGsMethylated) + getTotalCount(counts.snpsReferenceCpGsInterMethylated) + getTotalCount(counts.snpsReferenceCpGsUnMethylated);
}

/**
 * \brief Total CpG with SNP Called Quality 20
 */
unsigned int getTotalSnpsReferenceCpGsQuality20()
{
	return  counts.nonReferenceCpgMethylated[1] + counts.nonReferenceCpgInterMethylated[1] + counts.nonReferenceCpgUnMethylated[1];
}



/**
 * \brief Total CpGs With SNPs Methylated
 */
unsigned int getTotalSnpsReferenceCpGsMethylated(){ return getTotalCount(counts.snpsReferenceCpGsMethylated);}


/**
 * \brief Total CpG with SNP Called Methylated Quality 20
 */
unsigned int getTotalSnpsReferenceCpGsMethylatedQuality20(){ return counts.snpsReferenceCpGsMethylated[1]; }


/**
 * \brief Total CpGs With SNPs Intermediate Methylated
 */
unsigned int getTotalSnpsReferenceCpGsIntermediateMethylated(){ return getTotalCount(counts.snpsReferenceCpGsInterMethylated);}


/**
 * \brief Total CpGs With SNPs Intermediate Methylated Quality 20
 */
unsigned int getTotalSnpsReferenceCpGsIntermediateMethylatedQuality20(){ return counts.snpsReferenceCpGsInterMethylated[1]; }



/**
 * \brief Get Total Cpg Snp UnMethylated
 */
unsigned int getTotalSnpsReferenceCpGsUnMethylated(){	return getTotalCount(counts.snpsReferenceCpGsUnMethylated);}

/**
 * \brief Get Total Cpg Snp UnMethylated Quality 20
 */
unsigned int getTotalSnpsReferenceCpGsUnMethylatedQuality20(){ return counts.snpsReferenceCpGsUnMethylated[1];}


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
	printf("Total Homozygous: \t %i \t (%.2f %%) \n",getTotalHomozygous(),getPercentage(getTotalHomozygous(), getTotalDinucleotides()));
	printf("\t Total Homozygous Quality > 20: \t %i \t (%.2f %%) \n",getTotalHomozygousQuality20(),getPercentage(getTotalHomozygousQuality20(), getTotalDinucleotides()));
	printf("Total AlternativeCX: \t %i \t (%.2f %%) \n",getTotalAlternativeCX(),getPercentage(getTotalAlternativeCX(), getTotalDinucleotides()));
	printf("\t Total AlternativeCX Quality > 20: \t %i \t (%.2f %%) \n",getTotalAlternativeCXQuality20(),getPercentage(getTotalAlternativeCXQuality20(), getTotalDinucleotides()));

	printf("\n");
	printf("Total Methylated: \t %i \t (%.2f %%) \n", getTotalMethylated(),getPercentage(getTotalMethylated(), getTotalDinucleotides()));
	printf("\t Total Methylated Quality > 20: \t %i \t (%.2f %%) \n", getTotalMethylatedQuality20(),getPercentage(getTotalMethylatedQuality20(), getTotalDinucleotides()));
    printf("Total Intermediate Methylated: \t %i \t (%.2f %%) \n", getTotalIntermediateMethylated(),getPercentage(getTotalIntermediateMethylated(), getTotalDinucleotides()));
	printf("\t Total Intermediate Methylated Quality > 20: \t %i \t (%.2f %%) \n", getTotalIntermediateMethylatedQuality20(),getPercentage(getTotalIntermediateMethylatedQuality20(), getTotalDinucleotides()));
	printf("Total UnMethylated: \t %i \t (%.2f %%) \n", getTotalUnMethylated(),getPercentage(getTotalUnMethylated(), getTotalDinucleotides()));
	printf("\t Total Unmethylated Quality > 20: \t %i \t (%.2f %%) \n", getTotalUnMethylatedQuality20(),getPercentage(getTotalUnMethylatedQuality20(), getTotalDinucleotides()));

	printf("\n");
	printf("Total Homozygous Methylated: \t %i \t (%.2f %%)  \n", getTotalMethylatedHomozygous(),getPercentage(getTotalMethylatedHomozygous(),getTotalMethylated()));
	printf("\t Total Homozygous Methylated Quality > 20: \t %i \t (%.2f %%)  \n", getTotalMethylatedHomozygousQuality20(),getPercentage(getTotalMethylatedHomozygousQuality20(),getTotalMethylated()));
	printf("Total AlternativeCX Methylated: \t %i \t (%.2f %%)  \n", getTotalMethylatedAlternativeCX(),getPercentage(getTotalMethylatedAlternativeCX(),getTotalMethylated()));
	printf("\t Total AlternativeCX Methylated Quality > 20: \t %i \t (%.2f %%)  \n", getTotalMethylatedAlternativeCXQuality20(),getPercentage(getTotalMethylatedAlternativeCXQuality20(),getTotalMethylated()));

	printf("\n");
	printf("Total Homozygous Intermediate Methylated: \t %i \t (%.2f %%)  \n", getTotalIntermediateMethylatedHomozygous(),getPercentage(getTotalIntermediateMethylatedHomozygous(),getTotalIntermediateMethylated()));
	printf("\t Total Homozygous Intermediate Methylated Quality > 20: \t %i \t (%.2f %%)  \n", getTotalIntermediateMethylatedHomozygousQuality20(),getPercentage(getTotalIntermediateMethylatedHomozygousQuality20(),getTotalIntermediateMethylated()));
	printf("Total AlternativeCX Intermediate Methylated: \t %i \t (%.2f %%)  \n", getTotalIntermediateMethylatedAlternativeCX(),getPercentage(getTotalIntermediateMethylatedAlternativeCX(),getTotalIntermediateMethylated()));
	printf("\t Total AlternativeCX Intermediate Methylated Quality > 20: \t %i \t (%.2f %%)  \n", getTotalIntermediateMethylatedAlternativeCXQuality20(),getPercentage(getTotalIntermediateMethylatedAlternativeCXQuality20(),getTotalIntermediateMethylated()));

	printf("\n");
	printf("Total Homozygous UnMethylated: \t %i \t (%.2f %%) \n", getTotalUnMethylatedHomozygous(),getPercentage(getTotalUnMethylatedHomozygous(),getTotalUnMethylated()));
	printf("\t Total Homozygous UnMethylated Quality > 20: \t %i \t (%.2f %%) \n", getTotalUnMethylatedHomozygousQuality20(),getPercentage(getTotalUnMethylatedHomozygousQuality20(),getTotalUnMethylated()));
	printf("Total AlternativeCX UnMethylated: \t %i \t (%.2f %%) \n", getTotalUnMethylatedAlternativeCX(),getPercentage(getTotalUnMethylatedAlternativeCX(),getTotalUnMethylated()));
	printf("\t Total AlternativeCX UnMethylated Quality > 20: \t %i \t (%.2f %%) \n", getTotalUnMethylatedAlternativeCXQuality20(),getPercentage(getTotalUnMethylatedAlternativeCXQuality20(),getTotalUnMethylated()));

	printf("\n");
	printf("Total NonReference CpGs: \t %i \t (%.2f %%) \n", getTotalNonReferenceCpgs(),getPercentage(getTotalNonReferenceCpgs(),getTotalDinucleotides()));
	printf("\t Total NonReference CpGs Quality > 20: \t %i \t (%.2f %%) \n", getTotalNonReferenceCpgsQuality20(),getPercentage(getTotalNonReferenceCpgsQuality20(),getTotalDinucleotides()));

	printf("\t Total NonReference CpGs Methylated: \t %i \t (%.2f %%) \n",getTotalNonReferenceCpgsMethylated(),getPercentage(getTotalNonReferenceCpgsMethylated(),getTotalNonReferenceCpgs()));
	printf("\t\t Total NonReference CpGs Methylated Quality > 20: \t %i \t (%.2f %%) \n",getTotalNonReferenceCpgsMethylatedQuality20(),getPercentage(getTotalNonReferenceCpgsMethylatedQuality20(),getTotalNonReferenceCpgs()));

	printf("\t Total NonReference CpGs Intermediate Methylated: \t %i \t (%.2f %%) \n",getTotalNonReferenceCpgsIntermediateMethylated(),getPercentage(getTotalNonReferenceCpgsIntermediateMethylated(),getTotalNonReferenceCpgs()));
	printf("\t\t Total NonReference CpGs Intermediate Methylated Quality > 20: \t %i \t (%.2f %%) \n",getTotalNonReferenceCpgsIntermediateMethylatedQuality20(),getPercentage(getTotalNonReferenceCpgsIntermediateMethylatedQuality20(),getTotalNonReferenceCpgs()));

	printf("\t Total NonReference CpGs UnMethylated: \t %i \t (%.2f %%) \n",getTotalNonReferenceCpgsUnMethylated(),getPercentage(getTotalNonReferenceCpgsUnMethylated(),getTotalNonReferenceCpgs()));
	printf("\t\t Total NonReference CpGs UnMethylated Quality > 20: \t %i \t (%.2f %%) \n",getTotalNonReferenceCpgsUnMethylatedQuality20(),getPercentage(getTotalNonReferenceCpgsUnMethylatedQuality20(),getTotalNonReferenceCpgs()));

	printf("\n");
	printf("Total CpG with SNP Called: \t %i \t (%.2f %%) \n", getTotalSnpsReferenceCpGs(),getPercentage(getTotalSnpsReferenceCpGs(),getTotalDinucleotides()));
	printf("\t Total CpG with SNP Called Quality > 20: \t %i \t (%.2f %%) \n", getTotalSnpsReferenceCpGsQuality20(),getPercentage(getTotalSnpsReferenceCpGsQuality20(),getTotalDinucleotides()));

	printf("\t Total CpG with SNP Called Methylated: \t %i \t (%.2f %%) \n",getTotalSnpsReferenceCpGsMethylated(),getPercentage(getTotalSnpsReferenceCpGsMethylated(),getTotalSnpsReferenceCpGs()));
	printf("\t\t Total CpG with SNP Called Methylated Quality > 20: \t %i \t (%.2f %%) \n",getTotalSnpsReferenceCpGsMethylatedQuality20(),getPercentage(getTotalSnpsReferenceCpGsMethylatedQuality20(),getTotalSnpsReferenceCpGs()));

	printf("\t Total CpG with SNP Called Intermediate Methylated: \t %i \t (%.2f %%) \n",getTotalSnpsReferenceCpGsIntermediateMethylated(),getPercentage(getTotalSnpsReferenceCpGsIntermediateMethylated(),getTotalSnpsReferenceCpGs()));
	printf("\t\t Total CpG with SNP Called Intermediate Methylated Quality > 20: \t %i \t (%.2f %%) \n",getTotalSnpsReferenceCpGsIntermediateMethylatedQuality20(),getPercentage(getTotalSnpsReferenceCpGsIntermediateMethylatedQuality20(),getTotalSnpsReferenceCpGs()));

	printf("\t Total CpG with SNP Called UnMethylated: \t %i \t (%.2f %%) \n",getTotalSnpsReferenceCpGsUnMethylated(),getPercentage(getTotalSnpsReferenceCpGsUnMethylated(),getTotalSnpsReferenceCpGs()));
	printf("\t\t Total CpG with SNP Called UnMethylated Quality > 20: \t %i \t (%.2f %%) \n",getTotalSnpsReferenceCpGsUnMethylatedQuality20(),getPercentage(getTotalSnpsReferenceCpGsUnMethylatedQuality20(),getTotalSnpsReferenceCpGs()));

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

	fprintf(fp,"  \"TotalHomozygous\":%i,\n",getTotalHomozygous());
	fprintf(fp,"  \"TotalHomozygousQuality20\":%i,\n",getTotalHomozygousQuality20());

	fprintf(fp,"  \"TotalAlternativeCX\":%i,\n",getTotalAlternativeCX());
	fprintf(fp,"  \"TotalAlternativeCXQuality20\":%i,\n",getTotalAlternativeCXQuality20());

	fprintf(fp,"  \"TotalMethylated\":%i,\n", getTotalMethylated());
	fprintf(fp,"  \"TotalMethylatedQuality20\":%i,\n", getTotalMethylatedQuality20());
	fprintf(fp,"  \"TotalIntermediateMethylated\":%i,\n", getTotalIntermediateMethylated());
	fprintf(fp,"  \"TotalIntermediateMethylatedQuality20\":%i,\n", getTotalIntermediateMethylatedQuality20());
	fprintf(fp,"  \"TotalUnMethylated\":%i,\n", getTotalUnMethylated());
	fprintf(fp,"  \"TotalUnmethylatedQuality20\":%i,\n", getTotalUnMethylatedQuality20());

	fprintf(fp,"  \"TotalHomozygousMethylated\":%i,\n", getTotalMethylatedHomozygous());
	fprintf(fp,"  \"TotalHomozygousMethylatedQuality20\":%i,\n", getTotalMethylatedHomozygousQuality20());
	fprintf(fp,"  \"TotalAlternativeCXMethylated\":%i,\n", getTotalMethylatedAlternativeCX());
	fprintf(fp,"  \"TotalAlternativeCXMethylatedQuality20\":%i,\n", getTotalMethylatedAlternativeCXQuality20());

	fprintf(fp,"  \"TotalHomozygousIntermediateMethylated\":%i,\n", getTotalIntermediateMethylatedHomozygous());
	fprintf(fp,"  \"TotalHomozygousIntermediateMethylatedQuality20\":%i,\n", getTotalIntermediateMethylatedHomozygousQuality20());
	fprintf(fp,"  \"TotalAlternativeCXIntermediateMethylated\":%i,\n", getTotalIntermediateMethylatedAlternativeCX());
	fprintf(fp,"  \"TotalAlternativeCXIntermediateMethylatedQuality20\":%i,\n", getTotalIntermediateMethylatedAlternativeCXQuality20());

	fprintf(fp,"  \"TotalHomozygousUnMethylated\":%i,\n", getTotalUnMethylatedHomozygous());
	fprintf(fp,"  \"TotalHomozygousUnMethylatedQuality20\":%i,\n", getTotalUnMethylatedHomozygousQuality20());
	fprintf(fp,"  \"TotalAlternativeCXUnMethylated\":%i,\n", getTotalUnMethylatedAlternativeCX());
	fprintf(fp,"  \"TotalAlternativeCXUnMethylatedQuality20\":%i,\n", getTotalUnMethylatedAlternativeCXQuality20());

	fprintf(fp,"  \"TotalNonReferenceCpGs\":%i,\n", getTotalNonReferenceCpgs());
	fprintf(fp,"  \"TotalNonReferenceCpGsQuality20\":%i,\n", getTotalNonReferenceCpgsQuality20());

	fprintf(fp,"  \"TotalNonReferenceCpGsMethylated\":%i,\n",getTotalNonReferenceCpgsMethylated());
	fprintf(fp,"  \"TotalNonReferenceCpGsMethylatedQuality20\":%i,\n",getTotalNonReferenceCpgsMethylatedQuality20());

	fprintf(fp,"  \"TotalNonReferenceCpGsIntermediateMethylated\":%i,\n",getTotalNonReferenceCpgsIntermediateMethylated());
	fprintf(fp,"  \"TotalNonReferenceCpGsIntermediateMethylatedQuality20\":%i,\n",getTotalNonReferenceCpgsIntermediateMethylatedQuality20());

	fprintf(fp,"  \"TotalNonReferenceCpGsUnMethylated\":%i,\n",getTotalNonReferenceCpgsUnMethylated());
	fprintf(fp,"  \"TotalNonReferenceCpGsUnMethylatedQuality20\":%i,\n",getTotalNonReferenceCpgsUnMethylatedQuality20());

	fprintf(fp,"  \"TotalSnpsReferenceCpGsCalled\":%i,\n", getTotalSnpsReferenceCpGs());
	fprintf(fp,"  \"TotalSnpsReferenceCpGsCalledQuality20\":%i,\n", getTotalSnpsReferenceCpGsQuality20());

	fprintf(fp,"  \"TotalSnpsReferenceCpGsCalledMethylated\":%i,\n",getTotalSnpsReferenceCpGsMethylated());
	fprintf(fp,"  \"TotalSnpsReferenceCpGsCalledMethylatedQuality20\":%i,\n",getTotalSnpsReferenceCpGsMethylatedQuality20());

	fprintf(fp,"  \"TotalSnpsReferenceCpGsCalledIntermediateMethylated\":%i,\n",getTotalSnpsReferenceCpGsIntermediateMethylated());
	fprintf(fp,"  \"TotalSnpsReferenceCpGsCalledIntermediateMethylatedQuality20\":%i,\n",getTotalSnpsReferenceCpGsIntermediateMethylatedQuality20());

	fprintf(fp,"  \"TotalSnpsReferenceCpGsCalledUnMethylated\":%i,\n",getTotalSnpsReferenceCpGsUnMethylated());
	fprintf(fp,"  \"TotalSnpsReferenceCpGsCalledUnMethylatedQuality20\":%i\n",getTotalSnpsReferenceCpGsUnMethylatedQuality20());

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

	/*Vector of Methylation Value */
	unsigned int i = 0;

	fprintf(fp,"  \"MethValues\":[");

	for(i=0; i < vectorMethValues.size;i++)
	{
		if(i == (vectorMethValues.size -1))
		{
			fprintf(fp," %.3f",vector_get(&vectorMethValues,i));
		}
		else
		{
			fprintf(fp," %.3f,",vector_get(&vectorMethValues,i));
		}
	}

	fprintf(fp,"]\n");

	fprintf(fp,"}\n");


	fclose(fp);

	/*Close Vector*/
    vector_free(&vectorMethValues);
}















