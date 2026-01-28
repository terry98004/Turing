// -------------------------------------------------------------------
// Program last modified January 28, 2026. 
// Copyright (c) 2025-2026 Terrence P. Murphy
// MIT License -- see turing.c for details.
// -------------------------------------------------------------------

#include <stdio.h>
#include <math.h>			
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdbool.h>
#include <ctype.h>
#include <pthread.h>
#include <mpfr.h>

#include "hgt.h"
#include "turing.h"

extern struct	HGT_INIT	hgt_init;

int    CountZ;   	// number of Z(t) per Gram interval (excluding right endpoint)
struct GRAMLIST		gList[HGT_TUR_GRAM_PTS_MAX + 1];
struct HARDYINFO	hInfo[(HGT_TUR_GRAM_PTS_MAX * HGT_TUR_SUBINTVL_MAX) + 2];

// *******************************************************************
// We 
// *******************************************************************
int ComputeTuring(struct TURING tur)
{
mpfr_t				t, nOfGram, Accuracy, HardyZ, Temp1, Temp2;
int					i, j, hzNum;
int					MinusOneToN;
uint64_t			ui64N;

// -------------------------------------------------------------------
// Initialize the MPFR system and the variables that will hold the 
// Hardy Z remainder coefficients.
// -------------------------------------------------------------------
InitMPFR(tur.DefaultBits, tur.Threads, tur.DebugFlags, true);
InitCoeffMPFR(hgt_init.DefaultBits);	


mpfr_inits2 (hgt_init.DefaultBits, t, nOfGram, Accuracy, HardyZ, Temp1, Temp2, (mpfr_ptr) 0);
CountZ = tur.CountZ;	// at least for now, save globally for callbacks

// -------------------------------------------------------------------
// Do an MPFR initialization of the MPFR elemenrts in the GRAMLIST 
// structure.  We go "one past" our Gram count for reasons stated below.
// -------------------------------------------------------------------
for(i=0; i <= tur.CountGram; i++) {
	mpfr_inits2 (hgt_init.DefaultBits, gList[i].Gram, gList[i].n, 
		gList[i].lenInterval, gList[i].lenSubInterval, (mpfr_ptr) 0);
	}

// -------------------------------------------------------------------
// Convert 't' from a string to MPFR.  Then locate the largest 'n' 
// such that G(n) <= 't'.  As part of this, we set the accuracy
// desired (for all code below) in locating Gram points.
// -------------------------------------------------------------------
mpfr_set_str (t, tur.tBuf, 10, MPFR_RNDN);

mpfr_set_str (Accuracy, "0.1", 10, MPFR_RNDN);
mpfr_pow_ui(Accuracy, Accuracy, 16, MPFR_RNDN);
mpfr_div_ui(Accuracy, Accuracy, 2, MPFR_RNDN);

GramNearT(&nOfGram, t);
ui64N = mpfr_get_uj (nOfGram, MPFR_RNDN);
MinusOneToN = (ui64N % 2 == 0) ? 1 : -1;

// -------------------------------------------------------------------
// In our GRAMLIST structure, we save the Count (plus 1) Gram 
// points beginning with the given 'n' and also save each 'n'.  The 
// "plus 1" is because we also do this for the first 'n' AFTER our 
// Count of Gram points -- this allows us to calculate the interval 
// length of our last Gram point, and the Hardy Z value at that
// next Gram point.
// -------------------------------------------------------------------
for(i=0; i <= tur.CountGram; i++) {
	GramAtN(&gList[i].Gram, nOfGram, Accuracy);
	mpfr_set (gList[i].n, nOfGram, MPFR_RNDN);
	gList[i].MinusOneToN = MinusOneToN;
	MinusOneToN *= -1;
	mpfr_add_ui (nOfGram, nOfGram, 1, MPFR_RNDN);
	}

// -------------------------------------------------------------------
// We next update our GRAMLIST structure for the first Count 
// Gram points by calculating and saving the Gram interval lengths
// and sub-interval lengths.
// -------------------------------------------------------------------
for(i=0; i < tur.CountGram; i++) {
	mpfr_sub (gList[i].lenInterval, gList[i+1].Gram, gList[i].Gram, MPFR_RNDN);
	mpfr_div_si (gList[i].lenSubInterval, gList[i].lenInterval, tur.CountZ, MPFR_RNDN);
	}

// -------------------------------------------------------------------
// We next compute and  update our GRAMLIST structure with the Hardy Z 
// values at the interval points for that Gram point.  
//
// By way of example for a given Gram point, assume we have 8 interval
// point per Gram point.  The Hardy Z value at index 0 is the value at
// the Gram point.  Indexes 1 through 7 are the Hardy Z values at the
// remaining 7 interval points.  The 7th interval poiunt is 
// "lenSubInterval" before the next Gram point.
// -------------------------------------------------------------------
fprintf(stderr, "Processing Gram interval: ");
for(i=0; i < tur.CountGram; i++) {
	HardyZWithCount(gList[i].Gram, gList[i].lenSubInterval, tur.CountZ, i, HardyZCallbackA);
	fprintf(stderr, "%d..", i);
	}
fprintf(stderr, "\n\n");

// -------------------------------------------------------------------
// At this point, we need two more Hardy Z values: (1) the value at the
// point that is "lenSubInterval" before the first Gram point, and (2)
// the value at the first Gram point after our "CountGram" Gram points.
// -------------------------------------------------------------------
mpfr_sub (Temp1, gList[0].Gram, gList[0].lenSubInterval, MPFR_RNDN);
mpfr_set_ui(Temp2, 1, MPFR_RNDN);
HardyZWithCount(Temp1, Temp2, 1, 0, HardyZCallbackB);
HardyZWithCount(gList[tur.CountGram].Gram, Temp2, 1, tur.CountGram, HardyZCallbackA);

// -------------------------------------------------------------------
// We now have the information to determine whether our Gram points are 
// "good" or "bad". 
// -------------------------------------------------------------------
double	GoodTest;
for(i=0; i <= tur.CountGram; i++) {
	GoodTest = gList[i].MinusOneToN * hInfo[(i * CountZ) + 1].hzValue;
	gList[i].Good = GoodTest > 0 ? true : false;
	}

// -------------------------------------------------------------------
// With the "good" or "bad" information in hand, we now have the 
// information needed to determine whether each of our Gram intervals 
// has an even or odd number of zeros.
// -------------------------------------------------------------------
for(i=0; i < tur.CountGram; i++) {
	gList[i].OddZeros = gList[i].Good == gList[i+1].Good ? true : false;
	}

// -------------------------------------------------------------------
// Compute the rise (fall) between computed Hardy Z values. With
// that, compute whether the interval is moving towards or away
// from zero.
// -------------------------------------------------------------------
hInfo[0].hzRise = 0;
hzNum = (tur.CountGram * CountZ) + 1;
for (i = 1; i <= hzNum; i++) {
	hInfo[i].hzRise     = hInfo[i].hzValue - hInfo[i-1].hzValue;
	hInfo[i].TowardZero = hInfo[i].hzRise * hInfo[i].hzValue 
		> 0 ? false : true;
	}

// -------------------------------------------------------------------
// Locate all zero crossings between intervals.
// -------------------------------------------------------------------
int		idx;
double	TestCross;

hInfo[0].ZeroCross = hInfo[1].ZeroCross = false;
for(i=0; i < tur.CountGram; i++) {
	gList[i].ZerosFound = 0;
	idx = (i * CountZ) + 2;
	for(j = 0; j < tur.CountZ; j++) {
		TestCross = hInfo[idx + j].hzValue * hInfo[idx + j - 1].hzValue;
		if(TestCross < 0) {
			hInfo[idx + j].ZeroCross = true;
			gList[i].ZerosFound += 1;
			}
		else {
			hInfo[idx + j].ZeroCross = false;
			}
		}
	}

// -------------------------------------------------------------------
// Look for a change in slope that suggests we may be missing a zero
// crossing. (Or, we have a Lehmer failure).  NOTE: the below
// setting of hInfo[i].Lehmer works because Lehmer is a bool.
// -------------------------------------------------------------------
hInfo[0].Lehmer = false;
hzNum = (tur.CountGram * CountZ) + 1;
for (i = 1; i <= hzNum; i++) {
	hInfo[i].Lehmer =
		(hInfo[i].TowardZero == false 
		&& hInfo[i].ZeroCross == false  
		&& hInfo[i-1].TowardZero == true);
	}

TuringReport(tur);

mpfr_clears (t, nOfGram, Accuracy, HardyZ, Temp1, Temp2, (mpfr_ptr) 0);

for(i=0; i <= tur.CountGram; i++) {
	mpfr_clears (gList[i].Gram, gList[i].lenInterval, gList[i].n, 
	gList[i].lenSubInterval, (mpfr_ptr) 0);
	}

// -------------------------------------------------------------------
// Free the MPFR variables that held the coefficients.
// -------------------------------------------------------------------
CloseCoeffMPFR();
CloseMPFR();
return(1);	
}


// *******************************************************************
// This function provides our report of the Gram points and Hardy Z
// values needed to apply Turing's Method.
// *******************************************************************
int TuringReport(struct TURING tur)
{
int		i, j, idx;
mpfr_t	Temp1;

mpfr_inits2 (hgt_init.DefaultBits, Temp1, (mpfr_ptr) 0);

// -------------------------------------------------------------------
// In a verbose report, we show: (1) the 'n' and Gram point closst to 
// but less than or equal to the command line 't' value, (2) the 'n' and
// Gram point and Hardy Z values (and whether the Gram point is "good")
// for the "CountGram" (plus one) Gram points  that begin at the Gram 
// point associated with 'n', and (3) the interval length and sub-interval
// length for the "CountGram" Gram points.
// -------------------------------------------------------------------
if(tur.Verbose == true) {
	mpfr_printf("This report begins with n = %.0Rf and t = g(n) = %.*Rf\n\n", gList[0].n, 
		tur.OutputDP, gList[0].Gram);

	printf("We must find K = %u consecutive Gram blocks that satisfy Rosser's Rule\n\n", 
			ComputeTuringK(gList[0].Gram));

	for(i=0; i <= tur.CountGram; i++) {
		mpfr_printf("For n = %.0Rf, Gram = %.*Rf, -1^{n} = %2d, Hardy Z = %*.*f, Gram good = %s \n", 
			gList[i].n, tur.OutputDP, gList[i].Gram, gList[i].MinusOneToN, TUR_HARDY_WIDTH, 
			TUR_HARDY_DECIMALS, 
			hInfo[(i * CountZ) + 1].hzValue, 
			gList[i].Good == true ? "true" : "false" );						
		}
	
	printf("\nThe Gram interval lengths and sub-interval lengths are as follows:\n\n");
	for(i=0; i < tur.CountGram; i++) {
		mpfr_printf("For Gram = %.16Rf, Interval length = %.16Rf, SubInterval len = %.16Rf \n", 
			gList[i].Gram, gList[i].lenInterval, gList[i].lenSubInterval);
		}

// -------------------------------------------------------------------
// Continuing our verbose report, we show the number of zeros located
// in each Gram interval, and whether that count of zeros per interval 
// is "as expected" (an even or odd number).
// -------------------------------------------------------------------
	char	ExpectOdd[]  = "Missing at least one zero (expected an odd number)";
	char	ExpectEven[] = "Missing at least one zero (expected an even number)";
	char	AsExpected[] = "Even/Odd As Expected";
	char *  Message;

	printf("\n");
	for(i=0; i < tur.CountGram; i++) {
		Message = AsExpected;
		if(gList[i].Good == gList[i+1].Good) {  // so, expecting odd number of zeros
			if((gList[i].ZerosFound % 2) == 0)	// actual found = even
			Message = ExpectOdd;
			}
		else if((gList[i].ZerosFound % 2) != 0)	{ // expecting even, but actual found = odd
			Message = ExpectEven;
			}
		printf("Gram idx = %d, Zeros Found = %d, %s \n", i, gList[i].ZerosFound, Message);
		}
	printf("\nThe requested Turing Method data is as follows:\n\n");
	}

// -------------------------------------------------------------------
// We show the 't' and Hardy Z values of the sub-interval before
// the first Gram point.
// -------------------------------------------------------------------
mpfr_sub (Temp1, gList[0].Gram, gList[0].lenSubInterval, MPFR_RNDN);
mpfr_printf("G( 0) -1, %.*Rf, %*.*f \n", tur.OutputDP, Temp1, TUR_HARDY_WIDTH, 
			TUR_HARDY_DECIMALS, hInfo[0].hzValue);

// -------------------------------------------------------------------
// For the "CountGram" Gram points, we show the 't' and Hardy Z
// values of the "CountZ" sub-intervals. We also show: (1) whether 
// there was a "zero crossing" in the interval, and (2) whether
// there is a "Lehmer problem".
// -------------------------------------------------------------------
for(i=0; i < tur.CountGram; i++) {
	for(j = 0; j < tur.CountZ; j++) {
		mpfr_mul_si (Temp1, gList[i].lenSubInterval, j, MPFR_RNDN);
		mpfr_add (Temp1, gList[i].Gram, Temp1, MPFR_RNDN);
		mpfr_printf("G(%2d) %2d, %.*Rf, %*.*f,  %10.6f, %9s, %7s \n", i, j, 
			tur.OutputDP, Temp1, TUR_HARDY_WIDTH, TUR_HARDY_DECIMALS, 
			hInfo[(i * CountZ) + 1 +j ].hzValue,
			hInfo[(i * CountZ) + 1 +j ].hzRise,
			hInfo[(i * CountZ) + 1 +j ].ZeroCross == true ? "Crossing" : " ",
			hInfo[(i * CountZ) + 1 +j ].Lehmer == true ? "Lehmer" : " ");
		}
	}
// -------------------------------------------------------------------
// We show the 't' and Hardy Z values (and Crossing / Lehmer info) of 
//the Gram point that immediately follows our "CountGram" Gram points.
// -------------------------------------------------------------------	
idx = (tur.CountGram * CountZ) + 1;
mpfr_printf("G(%2d)  0, %.*Rf, %*.*f,  %10.6f, %9s, %7s \n", tur.CountGram, tur.OutputDP, 
			gList[tur.CountGram].Gram, TUR_HARDY_WIDTH, TUR_HARDY_DECIMALS, 
			hInfo[idx].hzValue,
			hInfo[idx].hzRise,
			hInfo[idx].ZeroCross == true ? "Crossing" : " ",
			hInfo[idx].Lehmer == true ? "Lehmer" : " ");

mpfr_clears (Temp1, (mpfr_ptr) 0);

return(1);	
}


// *******************************************************************
// We compute the K value used in Turing's Methods, where K is the 
// smallest positive integer such that:
//
// K >= 0.00313(log Gp)^2 + 0.1039 log Gp,
// 
// where Gp is the high part of the Gram segment [Gn, Gp)
// where Gn = the passed Gram value.
//
// Of course Gp is not exactly knowable, but should always be less than
// Gn + 100.
// by the HardyZWithCount library function.
// *******************************************************************
int ComputeTuringK(mpfr_t Gram)
{
mpfr_t			Gp, logGp, logGpSq, Temp1, Temp2;
unsigned int	K;

mpfr_inits2 (hgt_init.DefaultBits, Gp, logGp, logGpSq, Temp1, Temp2, (mpfr_ptr) 0);

mpfr_add_ui (Gp, Gram, (unsigned long int) 100, MPFR_RNDN);
mpfr_log (logGp, Gp, MPFR_RNDN);
mpfr_sqr (logGpSq, logGp, MPFR_RNDN);
mpfr_mul_d (Temp1, logGp, (double) 0.1039, MPFR_RNDN);
mpfr_mul_d (Temp2, logGpSq, (double) 0.00313, MPFR_RNDN);
mpfr_add (Temp1, Temp1, Temp2, MPFR_RNDN);
mpfr_rint_ceil (Temp1, Temp1, MPFR_RNDU);
K = (unsigned int) mpfr_get_ui (Temp1, MPFR_RNDN);
mpfr_clears (Gp, logGp, logGpSq, Temp1, Temp2, (mpfr_ptr) 0);
return(K);
}


// *******************************************************************
// Following are the two callback functions passed to and then called 
// by the HardyZWithCount library function.
//
// NOTE: To avoid the gcc compiler warning "warning unused parameter"
// we can use the gcc "__attribute__((unused))" specifier.  But
// we prefer the following method:
// 
// #pragma GCC diagnostic push
// #pragma GCC diagnostic ignored "-Wunused-parameter"
// <code with unused parameters here>
// #pragma GCC diagnostic pop
// *******************************************************************

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-parameter"
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

int HardyZCallbackA(mpfr_t t, mpfr_t HardyZ, int i, int CallerID) 
{
int	idx = (CallerID * CountZ) + 1 + i;
hInfo[idx].hzValue = mpfr_get_d (HardyZ, MPFR_RNDN);
return(1);
}


int HardyZCallbackB(mpfr_t t, mpfr_t HardyZ, int i, int CallerID) 
{
hInfo[0].hzValue = mpfr_get_d (HardyZ, MPFR_RNDN);
return(1);
}

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#pragma GCC diagnostic pop
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


/*
per Keith Thompson:
int64_t		Var;
printf("Var = %jd \n", (intmax_t)Var);
*/

// For testing:
// Lehmer pair: 7005.06266  7005.10056 -- also near 10854395965 and 35615956517
// Lehmer triplet: 12125271.88276506   12125272.07355435   12125272.2576189 