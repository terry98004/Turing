// -------------------------------------------------------------------
// Program last modified December 9, 2025. 
// Copyright (c) 2025 Terrence P. Murphy
// MIT License -- see turing.c for details.
// -------------------------------------------------------------------

#include <quadmath.h>
#define MPFR_WANT_FLOAT128 1

#include <time.h>
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

struct GRAMLIST	gList[HGT_TUR_GRAM_PTS_MAX + 1];


// *******************************************************************
// We 
// *******************************************************************
int ComputeTuring(struct TURING tur)
{
mpfr_t				t, nOfGram, Accuracy, HardyZ, Temp1, Temp2;
int					i;
int					MinusOneToN;
uint64_t			ui64N;

// -------------------------------------------------------------------
// Initialize the MPFR system and the variables that will hold the 
// Hardy Z remainder coefficients.
// -------------------------------------------------------------------
InitMPFR(tur.DefaultBits, tur.Threads, tur.DebugFlags, true);
InitCoeffMPFR(hgt_init.DefaultBits);	


mpfr_inits2 (hgt_init.DefaultBits, t, nOfGram, Accuracy, HardyZ, Temp1, Temp2, (mpfr_ptr) 0);

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
HardyZWithCount(gList[tur.CountGram].Gram, Temp2, 1, tur.CountGram, HardyZCallbackC);

// -------------------------------------------------------------------
// We now have the information to determine whether out Gram points are 
// "good" or "bad". We calculate that next.
// -------------------------------------------------------------------
double	GoodTest;
for(i=0; i <= tur.CountGram; i++) {
	GoodTest = gList[i].MinusOneToN * gList[i].HardyZValue[0];
	gList[i].Good = GoodTest > 0 ? 1 : -1;
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
int		i, j;
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
			TUR_HARDY_DECIMALS, gList[i].HardyZValue[0], gList[i].Good == 1 ? "true" : "false" );						
		}
	
	printf("\nThe Gram interval lengths and sub-interval lengths are as follows:\n\n");
	for(i=0; i < tur.CountGram; i++) {
		mpfr_printf("For Gram = %.16Rf, Interval length = %.16Rf, SubInterval len = %.16Rf \n", 
			gList[i].Gram, gList[i].lenInterval, gList[i].lenSubInterval);
		}
	printf("\nThe requested Turing Method data is as follows:\n\n");
	}

// -------------------------------------------------------------------
// We show the 't' and Hardy Z values of the sub-intervals before
// the first Gram point.
// -------------------------------------------------------------------
mpfr_sub (Temp1, gList[0].Gram, gList[0].lenSubInterval, MPFR_RNDN);
mpfr_printf("G( 0) -1, %.*Rf, %*.*f \n", tur.OutputDP, Temp1, TUR_HARDY_WIDTH, 
			TUR_HARDY_DECIMALS, gList[0].HardyZBefore);

// -------------------------------------------------------------------
// For the "CountGram" Gram points, we show the 't' and Hardy Z
// values of the "CountZ" sub-intervals.
// -------------------------------------------------------------------
for(i=0; i < tur.CountGram; i++) {
	for(j = 0; j < tur.CountZ; j++) {
		mpfr_mul_si (Temp1, gList[i].lenSubInterval, j, MPFR_RNDN);
		mpfr_add (Temp1, gList[i].Gram, Temp1, MPFR_RNDN);
		mpfr_printf("G(%2d) %2d, %.*Rf, %*.*f \n", i, j, 
			tur.OutputDP, Temp1, TUR_HARDY_WIDTH, TUR_HARDY_DECIMALS, 
			gList[i].HardyZValue[j]);
		}
	}

// -------------------------------------------------------------------
// We show the 't' and Hardy Z values of the Gram point that immediately
// follows our "CountGram" Gram points.
// -------------------------------------------------------------------	
mpfr_printf("G(%2d)  0, %.*Rf, %*.*f \n", tur.CountGram, tur.OutputDP, 
			gList[tur.CountGram].Gram, TUR_HARDY_WIDTH, TUR_HARDY_DECIMALS, 
			gList[tur.CountGram].HardyZValue[0]);

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
// Following are the three callback functions passed to and then called 
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
gList[CallerID].HardyZValue[i] = mpfr_get_d (HardyZ, MPFR_RNDN);
return(1);
}


int HardyZCallbackB(mpfr_t t, mpfr_t HardyZ, int i, int CallerID) 
{
gList[CallerID].HardyZBefore = mpfr_get_d (HardyZ, MPFR_RNDN);
return(1);
}



int HardyZCallbackC(mpfr_t t, mpfr_t HardyZ, int i, int CallerID) 
{
gList[CallerID].HardyZValue[0] = mpfr_get_d (HardyZ, MPFR_RNDN);
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
