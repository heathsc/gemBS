/*
 * parseArgs.h
 *
 *  Created on: 20 Mai, 2016
 *      Author: marcos
 */

#ifndef PARSEARGS_H_
#define PARSEARGS_H_

#define VERSION "0.1"

struct Args
{
	char * cpgInputFile;
	char * jsonFile;
	char * methJsonFile;
	char * bedFile;
	char * annotatedFile;
	int isZipped;
	char * firstCpGIsecFile;
	char * secondCpGIsecFile;
	int areIsecZipped;
};

int checkArguments(struct Args * arguments);

int getArgs (struct Args * arguments,int argc, char *argv[]);

void initArgs(struct Args * arguments);

#endif /* PARSEARGS_H_ */
