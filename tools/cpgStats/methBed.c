/*
 * methBed.c
 *
 *  Created on: 24 Mai, 2016
 *      Author: marcos
 */

#include "methBed.h"
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>

/**
 * \brief Get Mean Value
 * \param Get set of elements from which calculate a mean value
 * \param Array of methylation elements
 * \returns Mean Value
 */
float mean(unsigned int m, float * methValues)
{
    float sum;
    unsigned i;

    sum = 0.0;

    for(i=0; i<m; i++)
    {
        sum+=methValues[i];
    }

    return(sum/(float)m);
}

/**
 * \brief Median get Median value from a collection of methylation values
 * \param n number of methylation values
 * \param array of methylation values
 */
float median(unsigned int n, float * methValues)
{
    float temp;
    unsigned int i, j;
    // the following two loops sort the array x in ascending order
    for(i=0; i<n-1; i++) {
        for(j=i+1; j<n; j++) {
            if(methValues[j] < methValues[i]) {
                // swap elements
                temp = methValues[i];
                methValues[i] = methValues[j];
                methValues[j] = temp;
            }
        }
    }

    if(n%2==0)
    {
        // if there is an even number of elements, return mean of the two elements in the middle
        return((methValues[n/2] + methValues[n/2 - 1]) / 2.0);
    } else {
        // else return the element in the middle
        return methValues[n/2];
    }
}

/**
 * \brief Standard Deviation from a set of methylation values
 * \param n number of methylation values
 * \param array of methylation values
 */
float standardDeviation(int n,float * methValues)
{
    float mean, sumDeviation;
    int i;

    mean=0.0;
    sumDeviation=0.0;

    for(i=0; i<n;++i)
    {
        mean+=methValues[i];
    }

    mean=mean/n;
    for(i=0; i<n;++i)
    {
    	sumDeviation = sumDeviation + (methValues[i]-mean)*(methValues[i]-mean);
    }

    return sqrt(sumDeviation/(float)n);
}


/**
 * \brief Get Methylation Statistics
 * \param Windows in which region a set of stats should be calculated
 * \param nLen number of dinucleotides to calculate methylation values
 */
void getMethylationStats(struct Bed * window,unsigned int nLen,float vector[])
{
	window->meanMeth = mean(nLen,vector);
	window->medianMeth = median(nLen,vector);
	window->stDevMeth = standardDeviation(nLen,vector);
}

/**
 * \brief Get Windows Methylation Value
 * \param Dinucleotide first element in memory
 * \param contig Contig Name
 * \param Position start
 * \param Positione end
 * \param struct Bed * window
 */
void getWindowMethylation(Dinucleotide *head, struct Bed * window)
{
	Dinucleotide *current;

	unsigned int range = window->end - window->start;
	float methValues[range];
    unsigned int nDinucleotides;
    unsigned int nSNPs;

    nDinucleotides = 0;
    nSNPs = 0;

	if (head != NULL)
	{
        current = head;

        while (current!=NULL)
        {
            if ( (strcmp(current->contig, window->contig) == 0) && (current->position >= window->start) && (current->position <= window->end) )
            {
            	methValues[nDinucleotides] = current->methValue;
            	nDinucleotides = nDinucleotides + 1;

            	if (current->homozygous == 0)
    	        {
    	        	nSNPs = nSNPs + 1;
    	        }
            }
            else if((strcmp(current->contig, window->contig) != 0) || ((strcmp(current->contig, window->contig) == 0) && (current->position > window->end)) )
            {
        	    break;
            }

            current = current->next;
        }
	}

	window->cpgDinucleotides = nDinucleotides;
	window->snps = nSNPs;

    if (nDinucleotides == 0)
    {
    	window->meanMeth = -1;
    	window->medianMeth = -1;
    	window->stDevMeth = -1;
    }
    else
    {
    	getMethylationStats(window,nDinucleotides,methValues);
    }
}



/**
 * \brief Get Windows Methylation Value
 * \param Dinucleotide first element in memory
 */
void printNodes(Dinucleotide *head)
{
	Dinucleotide *previous;
	Dinucleotide *current;


	if (head != NULL)
	{
        previous = head;
        current = head->next;

        while (current!=NULL)
        {
            printf("[%s\t",current->contig);
            printf("%i\t",current->position);
            printf("%i\t",current->homozygous);
            printf("%.2f] --> \n",current->methValue);


            previous = current;
            current = current->next;
        }
	}

	printf("################################\n");
}


/**
 * \brief Find Dinucleotides in memory to be removed
 * \param Dinucleotide first element in memory
 * \param window reference to remove dinucleotides from main memory
 */
void toRemove (Dinucleotide ** head,struct Bed * window)
{
	Dinucleotide *current,*next;

	if((*head) == NULL)
	{
	    return;
	}

	current = *head;
	next = (*head)->next;

    /* Find rest nodes that should be removed */
    while (next!=NULL)
    {
    	if( ((strcmp(current->contig, window->contig) == 0) && (current->position < window->start)) || (strcmp(current->contig, window->contig) != 0) )
    	{
    		free(current->contig);
    		current->contig = NULL;
    		free(current);
    		current = next;
    		next = current->next;
    		*head = current;
    	}
    	else
    	{
        	return;
    	}
    }

    /*Remove last element*/
    if(next == NULL)
    {
    	if( ((strcmp(current->contig, window->contig) == 0) && (current->position < window->start)) || (strcmp(current->contig, window->contig) != 0) )
		{
			if (current != NULL)
			{
    		    free(current);
			    current = NULL;
			    *head = current;
			}
		}
    }
}

/**
 * \brief Create a new Node in the linked listed memory representation
 * \param Dinucleotide first element head
 * \param Dinucleotide current element
 * \param contig Contig Name
 * \param Position start
 * \param methValue Methylation Value
 * \param homozygous 1 If dinucleotide is homozygous otherwise 0
 */
void newNode(Dinucleotide** head,Dinucleotide** current,char * contig,unsigned int position,float methValue,int homozygous)
{
	Dinucleotide * newNode = (Dinucleotide*) malloc(sizeof(Dinucleotide));

	newNode->contig = strdup(contig);
    newNode->position = position;
	newNode->methValue = methValue;
	newNode->homozygous = homozygous;
	newNode->next = NULL;

	if (*head == NULL)
	{
		*head = newNode;
		*current = *head;
	}
	else
	{
		(*current)->next = newNode;
		*current = newNode;
	}
}

/****************************************************************************************************/
/*****************************                FILE OUTPUT                ****************************/
/****************************************************************************************************/

/**
 * \brief Check File Output
 * \param char * fileName File Name to write the set of data
 * \returns 1 if all was OK otherwise the program wil be quited
 */
int checkFileOutput(char * fileName)
{
	fileOutput = fopen(fileName, "w");

	if(fileOutput == NULL)
	{
	    printf("Sorry!! Something went wrong with file %s which outputs error: %s \n",fileName,strerror(errno));
	    exit(EXIT_FAILURE);
	}

	return 1;
}

/**
 * \brief Add Bed Window
 * \param contig Contig Name
 * \param Position start
 * \param Position End
 * \param extra bed fields
 * \param methValue Methylation Value
 * \param snps Number of snps in the window
 */
/*void addBedWindow(char * contig, unsigned int start, unsigned int end,char * extra,float methValue, unsigned int snps)*/
void addBedWindow(struct Bed * window)
{
	if(window->extra != NULL)
	{
		fprintf(fileOutput,"%s\t%i\t%i\t%s\t%.2f\t%.2f\t%.2f\t%i\t%i\n",
				window->contig,window->start,window->end,window->extra,window->meanMeth,window->medianMeth,window->stDevMeth,window->cpgDinucleotides,window->snps);
	}
	else
	{
		fprintf(fileOutput,"%s\t%i\t%i\t%.2f\t%.2f\t%.2f\t%i\t%i\n",
				window->contig,window->start,window->end,window->meanMeth,window->medianMeth,window->stDevMeth,window->cpgDinucleotides,window->snps);
	}

}

/**
 * \brief Close File Output
 */
void closeFileOutput()
{
	fclose(fileOutput);
}
