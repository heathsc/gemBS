/*
 * counts.h
 *
 *  Created on: 20 Mai, 2016
 *      Author: marcos
 */

#ifndef COUNTS_H_
#define COUNTS_H_

#include "common.h"
#include <stdio.h>

struct Counts {
   unsigned int referenceCGsMethylated[2];
   unsigned int referenceCGsInterMethylated[2];
   unsigned int referenceCGsUnMethylated[2];

   unsigned int nonReferenceCGsMethylated[2];
   unsigned int nonReferenceCGsInterMethylated[2];
   unsigned int nonReferenceCGsUnMethylated[2];
} counts;

//Declare a new vector
Vector vReferenceCGsMethValuesQ20;
Vector vNonReferenceCGsMethValuesQ20;

VectorInfoReads vInfoReadsReferenceCGs;
VectorInfoReads vInfoReadsNonReferenceCGs;
VectorInfoReads vInfoReadsReferenceCGsQ20;
VectorInfoReads vInfoReadsNonReferenceCGsQ20;

void initCounts();
void addRecordStats(struct Record * record);

void updateInfoReadsRegister(VectorInfoReads *vector,unsigned int infoReads);

unsigned int getTotalReferenceCGs();
unsigned int getTotalReferenceCGsQuality20();
unsigned int getTotalNonReferenceCGs();
unsigned int getTotalNonReferenceCGsQuality20();

unsigned int getTotalMethylated();
unsigned int getTotalMethylatedQuality20();

unsigned int getTotalIntermediateMethylated();
unsigned int getTotalIntermediateMethylatedQuality20();

unsigned int getTotalUnMethylated();
unsigned int getTotalUnMethylatedQuality20();

unsigned int getTotalMethylatedReferenceCGs();
unsigned int getTotalMethylatedReferenceCGsQuality20();

unsigned int getTotalIntermediateMethylatedReferenceCGs();
unsigned int getTotalIntermediateMethylatedReferenceCGsQuality20();

unsigned int getTotalUnMethylatedReferenceCGs();
unsigned int getTotalUnMethylatedReferenceCGsQuality20();

unsigned int getTotalMethylatedNonReferenceCGs();
unsigned int getTotalMethylatedNonReferenceCGsQuality20();

unsigned int getTotalIntermediateMethylatedNonReferenceCGs();
unsigned int getTotalIntermediateMethylatedNonReferenceCGsQuality20();

unsigned int getTotalUnMethylatedNonReferenceCGs();
unsigned int getTotalUnMethylatedNonReferenceCGsQuality20();

unsigned int getTotalQuality20();

float getPercentage(unsigned int concept, unsigned int totalValue);

void printCounts();
void saveCounts(char * fileName);
void saveJsonMethylationCounts(char * fileName);

void printVectorInformationReads(FILE * fp,VectorInfoReads * vector,char * text);
void saveJsonInformationReads(char * fileName);


#endif /* COUNTS_H_ */
