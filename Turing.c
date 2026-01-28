// -------------------------------------------------------------------
// Program last modified January 28, 2026. 
// -------------------------------------------------------------------

/*
MIT License

Copyright (c) 2025-2026 Terrence P. Murphy

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <time.h>
#include <stdio.h>
#include <math.h>			
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
#include <pthread.h>
#include <unistd.h>			// for command line argument processing
#include <mpfr.h>

#include "hgt.h"
#include "turing.h"


const char sUsage[] = "Command Line Parameters\n" \
 "-t [positive number]	Location of t along the critical line - this parameter is required. (Digits and '.' only).\n" \
 "-g [positive integer]	Count of the number of Gram intervals to check - between 1 and 24, defaults to 8.\n" \
 "-c [positive integer]	Count of the number of Z(t) values to check in a Gram interval - between 8 and 128, defaults to 8.\n" \
 "-p [positive integer]	Decimal point digits of 't values to show in report - between 2 and 60, defaults to 6.\n" \
 "-b [positive integer]	Floating point bits: 128 <= b <= 1024 - defaults to 256.\n" \
 "-d [positive integer]	Used for debugging only.  Please disregard.\n" \
 "-k [positive integer]	Number of threads to use - defaults to 1, maximum of 8.\n" \
 "-h			Show command line parameters.  All other parameters will be ignored.\n" \
 "-s			Report the total seconds taken to compute the Hardy Z values.\n"\
 "-v			Verbose report (provides additional useful information -- highly recommended)."; 

const char sCopyright[] = "Copyright 2025-2026 by Terrence P. Murphy." \
" Licensed under MIT License.\n\n"; 


int main( int argc, char *argv[] )  
{
int				c;
struct TURING	tur;
int				tDecimalDigits = -1;

tur.Verbose 	= false;
tur.ShowSeconds	= false;
tur.CountZ 		= 8;
tur.CountGram 	= 8;
tur.DebugFlags	= 2311; // (2 * 3 * 5 * 7 * 11) + 1
tur.OutputDP 	= 6;
tur.DefaultBits	= HGT_PRECISION_DEFAULT;
tur.Threads		= 1;

// strcpy(tur.incrBuf, "1");

fprintf(stderr, "%s", sCopyright);
if(argc == 1) {
	printf("%s\n", sUsage);
	exit(EXIT_FAILURE);
	}
opterr = 0; // To prevent _getopt from printing an error message on standard error

while ((c = getopt (argc, argv, "t:g:c:k:p:b:d:hvs")) != -1)
	switch (c)
		{
		case 'h':
			printf("%s\n", sUsage);		
//			intmax_t	iMax;
//			int64_t		i64;
//			printf("Sizeof intmax_t: %zu, Sizeof int64_t: %zu", sizeof(iMax), sizeof(i64));
			exit(EXIT_SUCCESS);
		case 'v':
			tur.Verbose = true;
			break;
		case 's':
			tur.ShowSeconds = true;
			break;			
		case 't':
			if (ValidateHardyT (optarg) < 1) {
				printf("Invalid argument to -t \n");
				return(EXIT_FAILURE);
				}
			strcpy(tur.tBuf, optarg);
			tDecimalDigits =  GetDecimalDigits(tur.tBuf);
			break;
		case 'g':
			tur.CountGram = ValidateTuringGramPoints(optarg);
			if(tur.CountGram < 1){
				printf("Invalid argument to -c \n");
				return(EXIT_FAILURE);
				}
			break;
		case 'c':
			tur.CountZ = ValidateTuringSubIntervals(optarg);
			if(tur.CountZ < 1){
				printf("Invalid argument to -c \n");
				return(EXIT_FAILURE);
				}
			break;
		case 'k':
			tur.Threads = ValidateThreads(optarg);	
			if(tur.Threads < 1){
				printf("Invalid argument to -k \n");
				return(EXIT_FAILURE);
				}
			break;	
		case 'd':
			tur.DebugFlags = ValidateDebugFlags(optarg);	
			if( tur.DebugFlags < 1){
				printf("Invalid argument to -d \n");
				return(EXIT_FAILURE);
				}
			break;
		case 'p':
			tur.OutputDP = ValidateReportDecimalPlaces(optarg);
			if( tur.OutputDP < 1){
				printf("Invalid argument to -z \n");
				return(EXIT_FAILURE);
				}
			break;
		case 'b':
			tur.DefaultBits = ValidatePrecisionMPFR(optarg);	
			if(tur.DefaultBits < 1) {
				printf("Invalid argument to -b \n");
				return(EXIT_FAILURE);
				}
			break;
		case '?':
		case ':':
		default:
			printf("Option -%c is either unknown or missing its argument\n", optopt);
			return (EXIT_FAILURE);
		}
if(tDecimalDigits == -1) {
	printf("The t parameter is required.\n");
	return(EXIT_FAILURE);
	}
else if(tDecimalDigits > tur.OutputDP) {
	tur.OutputDP	= tDecimalDigits;
	}

// -------------------------------------------------------------------
// We have finished validating the command line parameters.  Now compute 
// and printf our HardyZ results.  Report the time it takes to do the
// computations if the user enters the -s command line parameter.
// -------------------------------------------------------------------
struct timespec start, end;
double 			time_taken;

clock_gettime(CLOCK_MONOTONIC, &start);
ComputeTuring(tur);
clock_gettime(CLOCK_MONOTONIC, &end);
time_taken =  (end.tv_sec - start.tv_sec);
time_taken += (end.tv_nsec - start.tv_nsec) / 1000000000.0;

if(tur.ShowSeconds){
	printf("Compute took %f seconds to execute \n", time_taken); 
	}
return(EXIT_SUCCESS);
}
