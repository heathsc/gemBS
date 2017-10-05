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
   unsigned int homozygousMethylated[4];
   unsigned int homozygousInterMethylated[4];
   unsigned int homozygousUnMethylated[4];

   unsigned int heterozygousMethylated[4];
   unsigned int heterozygousInterMethylated[4];
   unsigned int heterozygousUnMethylated[4];

   unsigned int deNovoCpgMethylated[4];
   unsigned int deNovoCpgInterMethylated[4];
   unsigned int deNovoCpgUnMethylated[4];

   unsigned int cpgWithSnpMethylated[4];
   unsigned int cpgWithSnpInterMethylated[4];
   unsigned int cpgWithSnpUnMethylated[4];
} counts;

//Declare a new vector
Vector vectorMethValues;


void initCounts();
void addRecordStats(struct Record * record);

unsigned int getTotalHomozygous();
unsigned int getTotalHomozygousHighQuality();
unsigned int getTotalHeterozygous();
unsigned int getTotalHeterozygousHighQuality();

unsigned int getTotalMethylated();
unsigned int getTotalMethylatedHighQuality();

unsigned int getTotalIntermediateMethylated();
unsigned int getTotalIntermediateMethylatedHighQuality();

unsigned int getTotalUnMethylated();
unsigned int getTotalUnMethylatedHighQuality();

unsigned int getTotalMethylatedHomozygous();
unsigned int getTotalMethylatedHomozygousHighQuality();

unsigned int getTotalIntermediateMethylatedHomozygous();
unsigned int getTotalIntermediateMethylatedHomozygousHighQuality();

unsigned int getTotalUnMethylatedHomozygous();
unsigned int getTotalUnMethylatedHomozygousHighQuality();

unsigned int getTotalMethylatedHeterozygous();
unsigned int getTotalMethylatedHeterozygousHighQuality();

unsigned int getTotalIntermediateMethylatedHeterozygous();
unsigned int getTotalIntermediateMethylatedHeterozygousHighQuality();

unsigned int getTotalUnMethylatedHeterozygous();
unsigned int getTotalUnMethylatedHeterozygousHighQuality();

unsigned int getTotalQualityUnder10();
unsigned int getTotalQualityBetween10_20();
unsigned int getTotalQualityBetween20_30();
unsigned int getTotalHighQuality();

/*De Novo Cpgs Stats*/
unsigned int getTotalDeNovoCpgs();
unsigned int getTotalDeNovoCpgsHighQuality();
unsigned int getTotalDeNovoCpgsQuality_over20();

unsigned int getTotalDeNovoCpgsMethylated();
unsigned int getTotalDeNovoCpgsMethylatedHighQuality();
unsigned int getTotalDeNovoCpgsMethylatedQuality_over20();

unsigned int getTotalDeNovoCpgsIntermediateMethylated();
unsigned int getTotalDeNovoCpgsIntermediateMethylatedHighQuality();
unsigned int getTotalDeNovoCpgsIntermediateMethylatedQuality_over20();

unsigned int getTotalDeNovoCpgsUnMethylated();
unsigned int getTotalDeNovoCpgsUnMethylatedHighQuality();
unsigned int getTotalDeNovoCpgsUnMethylatedQuality_over20();


/*Reference CpG With Snps Detected*/
unsigned int getTotalCpgSnp();
unsigned int getTotalCpgSnpHighQuality();
unsigned int getTotalCpgSnpQuality_over20();

unsigned int getTotalCpgSnpMethylated();
unsigned int getTotalCpgSnpMethylatedHighQuality();
unsigned int getTotalCpgSnpMethylatedQuality_over20();

unsigned int getTotalCpgSnpIntermediateMethylated();
unsigned int getTotalCpgSnpIntermediateMethylatedHighQuality();
unsigned int getTotalCpgSnpIntermediateMethylatedQuality_over20();

unsigned int getTotalCpgSnpUnMethylated();
unsigned int getTotalCpgSnpUnMethylatedHighQuality();
unsigned int getTotalCpgSnpUnMethylatedQuality_over20();


float getPercentage(unsigned int concept, unsigned int totalValue);

void printCounts();
void saveCounts(char * fileName);
void saveJsonMethylationCounts(char * fileName);



#endif /* COUNTS_H_ */
