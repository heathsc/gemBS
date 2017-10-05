/*
 * methBed.h
 *
 *  Created on: 24 Mai, 2016
 *      Author: marcos
 */

#ifndef METHBED_H_
#define METHBED_H_

#include <stdio.h>
#include "common.h"

typedef struct dinucleotideInfo Dinucleotide;

struct dinucleotideInfo
{
	char * contig;             /*Contig name*/
	unsigned int position;     /*Position at contig*/
	float methValue;           /*Methylation Value*/
	int homozygous;            /*1 If dinucleotide is homozygous otherwise 0*/
	Dinucleotide * next;       /*ptr to next car in list */
};

void getWindowMethylation(Dinucleotide *head, struct Bed * window);
void toRemove (Dinucleotide ** head,struct Bed * window);
void newNode(Dinucleotide** head,Dinucleotide** current,char * contig,unsigned int position,float methValue,int homozygous);

FILE *fileOutput;

int checkFileOutput(char * fileName);
void addBedWindow(struct Bed * window);
void closeFileOutput();

void printNodes(Dinucleotide *head);

void getMethylationStats(struct Bed * window,unsigned int nLen,float vector[] );

#endif /* METHBED_H_ */
