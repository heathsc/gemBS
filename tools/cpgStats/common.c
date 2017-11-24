/*
 * common.c
 *
 *  Created on: 9 Mai, 2017
 *      Author: marcos
 */


#include <stdio.h>
#include <stdlib.h>
#include "common.h"

void vector_init(Vector *vector) {
  // initialize size and capacity
  vector->size = 0;
  vector->capacity = VECTOR_INITIAL_CAPACITY;

  // allocate memory for vector->data
  vector->data = malloc(sizeof(float) * vector->capacity);
}

void vector_append(Vector *vector, float value) {
  // make sure there's room to expand into
  vector_double_capacity_if_full(vector);

  // append the value and increment vector->size
  vector->data[vector->size++] = value;
}

float vector_get(Vector *vector, int index) {
  if (index >= vector->size || index < 0) {
    printf("Index %d out of bounds for vector of size %d\n", index, vector->size);
    exit(1);
  }
  return vector->data[index];
}

void vector_set(Vector *vector, int index, float value) {
  // zero fill the vector up to the desired index
  while (index >= vector->size) {
    vector_append(vector, 0);
  }

  // set the value at the desired index
  vector->data[index] = value;
}

void vector_double_capacity_if_full(Vector *vector) {
  if (vector->size >= vector->capacity) {
    // double vector->capacity and resize the allocated memory accordingly
    vector->capacity *= 2;
    vector->data = realloc(vector->data, sizeof(float) * vector->capacity);
  }
}

void vector_free(Vector *vector) {
  free(vector->data);
}


/**********************************
 *****VECTOR INFORMATION READS*****
 **********************************/

void info_reads_init(VectorInfoReads *vector) {
  // initialize size and capacity
  vector->size = 0;
  vector->capacity = VECTOR_INITIAL_CAPACITY;

  // allocate memory for vector->data
  vector->data = malloc(vector->capacity * sizeof(InfoReads));
}

void info_reads_append(VectorInfoReads *vector, InfoReads value) {
  // make sure there's room to expand into
  info_reads_double_capacity_if_full(vector);

  // append the value and increment vector->size
  vector->data[vector->size++] = value;
}

InfoReads info_reads_get(VectorInfoReads *vector, int index) {
  if (index >= vector->size || index < 0) {
    printf("Index %d out of bounds for info reads vector of size %d\n", index, vector->size);
    exit(1);
  }

  return vector->data[index];
}

void info_reads_set(VectorInfoReads *vector, int index, InfoReads value) {
  // zero fill the vector up to the desired index
  while (index >= vector->size) {
	  InfoReads current;
	  current.cpgs = 0;
	  current.informationReads = 0,
	  info_reads_append(vector, current);
  }

  // set the value at the desired index
  vector->data[index] = value;
}

void info_reads_double_capacity_if_full(VectorInfoReads *vector) {
  if (vector->size >= vector->capacity) {
    // double vector->capacity and resize the allocated memory accordingly
    vector->capacity *= 2;
    vector->data = realloc(vector->data, vector->capacity * sizeof(InfoReads));
  }
}

void info_reads_free(VectorInfoReads *vector) {
	free(vector->data);
}





