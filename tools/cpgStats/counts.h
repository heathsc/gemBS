/*
 * counts.h
 *
 *  Created on: 20 Mai, 2016
 *      Author: marcos
 */

#ifndef COUNTS_H_
#define COUNTS_H_

#include "common.h"

struct Counts {
   unsigned int homozygousMethylated[2];
   unsigned int homozygousInterMethylated[2];
   unsigned int homozygousUnMethylated[2];

   unsigned int alternativeCXMethylated[2];
   unsigned int alternativeCXInterMethylated[2];
   unsigned int alternativeCXUnMethylated[2];

   unsigned int nonReferenceCpgMethylated[2];
   unsigned int nonReferenceCpgInterMethylated[2];
   unsigned int nonReferenceCpgUnMethylated[2];

   unsigned int snpsReferenceCpGsMethylated[2];
   unsigned int snpsReferenceCpGsInterMethylated[2];
   unsigned int snpsReferenceCpGsUnMethylated[2];
} counts;

//Declare a new vector
Vector vectorMethValues;


void initCounts();
void addRecordStats(struct Record * record);

unsigned int getTotalHomozygous();
unsigned int getTotalHomozygousQuality20();
unsigned int getTotalAlternativeCX();
unsigned int getTotalAlternativeCXQuality20();

unsigned int getTotalMethylated();
unsigned int getTotalMethylatedQuality20();

unsigned int getTotalIntermediateMethylated();
unsigned int getTotalIntermediateMethylatedQuality20();

unsigned int getTotalUnMethylated();
unsigned int getTotalUnMethylatedQuality20();

unsigned int getTotalMethylatedHomozygous();
unsigned int getTotalMethylatedHomozygousQuality20();

unsigned int getTotalIntermediateMethylatedHomozygous();
unsigned int getTotalIntermediateMethylatedHomozygousQuality20();

unsigned int getTotalUnMethylatedHomozygous();
unsigned int getTotalUnMethylatedHomozygousQuality20();

unsigned int getTotalMethylatedAlternativeCX();
unsigned int getTotalMethylatedAlternativeCXQuality20();

unsigned int getTotalIntermediateMethylatedAlternativeCX();
unsigned int getTotalIntermediateMethylatedAlternativeCXQuality20();

unsigned int getTotalUnMethylatedAlternativeCX();
unsigned int getTotalUnMethylatedAlternativeCXQuality20();

unsigned int getTotalQuality20();

/*Non References Cpgs Stats*/
unsigned int getTotalNonReferenceCpgs();
unsigned int getTotalNonReferenceCpgsQuality20();

unsigned int getTotalNonReferenceCpgsMethylated();
unsigned int getTotalNonReferenceCpgsMethylatedQuality20();

unsigned int getTotalNonReferenceCpgsIntermediateMethylated();
unsigned int getTotalNonReferenceCpgsIntermediateMethylatedQuality20();

unsigned int getTotalNonReferenceCpgsUnMethylated();
unsigned int getTotalNonReferenceCpgsUnMethylatedQuality20();

/*SNPs (CX) at Reference CpGs*/
unsigned int getTotalSnpsReferenceCpGs();
unsigned int getTotalSnpsReferenceCpGsQuality20();

unsigned int getTotalSnpsReferenceCpGsMethylated();
unsigned int getTotalSnpsReferenceCpGsMethylatedQuality20();

unsigned int getTotalSnpsReferenceCpGsIntermediateMethylated();
unsigned int getTotalSnpsReferenceCpGsIntermediateMethylatedQuality20();

unsigned int getTotalSnpsReferenceCpGsUnMethylated();
unsigned int getTotalSnpsReferenceCpGsUnMethylatedQuality20();

float getPercentage(unsigned int concept, unsigned int totalValue);

void printCounts();
void saveCounts(char * fileName);
void saveJsonMethylationCounts(char * fileName);



#endif /* COUNTS_H_ */
