/*
 * common.h
 *
 *  Created on: 20 Mai, 2016
 *      Author: marcos
 */

#ifndef COMMON_H_
#define COMMON_H_

struct Record
{
	char * contig;             /*Contig name*/
	unsigned int position;     /*Position at contig*/
	char * referenceContext;   /*Dinucleotide base call according to reference*/
	char * callContext;        /*Dinucleotide base call according to genotype call*/
	unsigned int phredScore;   /*Probability call*/
	float methValue;           /*Methylation Value*/
	float methDev;             /*Methilation Deviation*/
	int noValue;               /*1 If there is no methylation Information*/
};

struct Bed
{
	char * contig;                   /*Contig name*/
	unsigned int start;              /*Star Position*/
	unsigned int end;                /*End Position*/
	char * extra;                    /*Extra BED Fields*/
	float meanMeth;                  /*Mean Meth*/
	float medianMeth;                /*Median Meth*/
	float stDevMeth;                 /*St Dev Methilation*/
	unsigned int cpgDinucleotides;   /*Number of CpG Dinucleotides In a given region*/
	unsigned int snps;               /*Number of snps in the bed Region*/
};


/**********************************
 *********VECTOR DEFINITIONS*******
 **********************************/

#define VECTOR_INITIAL_CAPACITY 100

// Define a vector type
typedef struct {
  int size;      // slots used so far
  int capacity;  // total available slots
  float *data;     // array of floats we're storing
} Vector;

/**
 *\brief vector_init is a function that initializes a vector struct.
 *\brief It sets size to 0, capacity to VECTOR_INITIAL_CAPACITY and allocates an appropriate amount of memory (vector->capacity * sizeof(int))
 *\brief for the underlying data array.
*/
void vector_init(Vector *vector);

/**
 *\brief vector_append appends the given value to the vector. If the underlying data array is full, then calling this function should cause
 *\brief vector->data to expand to accept this value. Increments vector->size.
*/
void vector_append(Vector *vector, float value);

/**
 *\brief vector_get returns a value out of a vector at the given index. If the index is below 0 or greater than vector->size - 1, this function
 *\brief should complain about the index being out of bounds.
*/
float vector_get(Vector *vector, int index);

/**
 *\brief vector_set sets the value at the given index to the given value. If the index is greater than the vector->size, this function should
 *\brief expand the vector until it is big enough to contain the index and set the value at that index. It should zero-fill all values in between.
 *\brief vector->size should be incremented accordingly.
*/
void vector_set(Vector *vector, int index, float value);

/**
 *\brief vector_double_capacity_if_full doubles the underlying data array capacity if vector->size >= vector->capacity.
 *\brief We'll find out later that changing the size of the array is expensive, so in order to minimize the number of times we need to resize,
 *\brief we double the capacity each time.
*/
void vector_double_capacity_if_full(Vector *vector);

/**
 *\brief vector_free frees the memory allocated for the data array.
 *\brief We leave freeing of the Vector struct itself to client code
 *\brief (so they can use any sort of pointer they like, be it stack or heap, and then clean up after themselves).
*/
void vector_free(Vector *vector);



#endif /* COMMON_H_ */
