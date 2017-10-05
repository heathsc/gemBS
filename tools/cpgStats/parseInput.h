/*
 * parseInput.h
 *
 *  Created on: 20 Mai, 2016
 *      Author: marcos
 */

#ifndef PARSEINPUT_H_
#define PARSEINPUT_H_

#include <stdio.h>
#include <zlib.h>
#include "common.h"

#define LENGTH 0x1000

FILE *inputData;

void initRecord(struct Record * record);

int configInput(char * inputName,int isZip);
void closeInput(int isZip);
int getNewRecord(struct Record * record,int isZip);

int textConfigInput(char * inputName);
void textCloseInput();
int textGetNewRecord(struct Record * record);

gzFile * gzInputFile;
int gzConfigInput(char * inputName);
void gzCloseInput();
int gzGetNewRecord(struct Record * record);


FILE *inputBed;

void chomp(const char *s);
void initBed(struct Bed * window);
int checkInputBed(char * inputName);
void closeInputBed();
int getNewWindow(struct Bed * window);


int setupInput(char * inputName,void ** fileDescriptor,int isZip);
void fileCloseInput(void * fileDescriptor,int isZip);
int readNewRecord(void * fileDescriptor,struct Record * record,int isZip);

#endif /* PARSEINPUT_H_ */
