/*
 * Intersection.h
 *
 *  Created on: 6 Xuñ, 2016
 *      Author: marcos
 */

#ifndef INTERSECTION_H_
#define INTERSECTION_H_

#include "common.h"
#define STACK_SIZE 100

struct ResultsIsec
{
	/* TOTAL PROCESSED*/
	unsigned int n_total_high_quality_one;
	unsigned int n_total_low_quality_one;
	unsigned int n_total_high_quality_two;
	unsigned int n_total_low_quality_two;

	/* SHARED */
	unsigned int n_shared_high_quality_one;
	unsigned int n_shared_low_quality_one;
	unsigned int n_shared_high_quality_two;
	unsigned int n_shared_low_quality_two;

	/* PRIVATE ONE*/
	unsigned int n_private_high_quality_one;
	unsigned int n_private_low_quality_one;

	/* PRIVATE SECOND*/
	unsigned int n_private_high_quality_two;
	unsigned int n_private_low_quality_two;
} resultsIsec;

struct StackPosition
{
	int hasValue;
	struct Record * record;
};

struct StackPosition stackInfo [STACK_SIZE];

void initIsec();
void initStack();

void printIsecResults(char * nameFileOne,char * nameFileTwo);

int runIsec(char * nameFileOne, char * nameFileTwo,int isGzip);


#endif /* INTERSECTION_H_ */
