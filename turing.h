// -------------------------------------------------------------------
// Program last modified January 9, 2026. 
// Copyright (c) 2025-2026 Terrence P. Murphy
// MIT License -- see turing.c for details.
// -------------------------------------------------------------------

#define	TUR_HARDY_WIDTH		15
#define	TUR_HARDY_DECIMALS	10

struct TURING {
	char	tBuf[HGT_MAX_CMDLINE_STRLEN + 2];	// holds entered '-t' value
	int 	CountZ;    	 			// number of Z(t) in Gram interval
	int 	CountGram;   			// number of Gram interval	
	bool	Verbose;				// verbose report? true or false
	bool	ShowSeconds;			// report compute second taken T/F
	int		DebugFlags;				// used for debugging
	int		OutputDP;				// digits after '.' in report
	int		DefaultBits;			// Bits for MPFR floating point 
	int		Threads;				// Number of threads to use 
}; 


struct GRAMLIST {
	mpfr_t		Gram; 			// 
	mpfr_t		n;				// g(n) = Gram 
	mpfr_t		lenInterval; 	// Gram(next) - Gram 
	mpfr_t		lenSubInterval;	// lenInterval / CountZ
	int			MinusOneToN;	// 1 if n even, -1 if n odd
	int			Good;			// true if (-1)^n * HardyZ[0] > 0
	double		HardyZBefore;
	double		HardyZValue[HGT_TUR_SUBINTVL_MAX + 2];
	double		IntervalSlope[HGT_TUR_SUBINTVL_MAX + 1];
	int			iFloatBits;		// Bits for MPFR floating point 
	int			DebugFlagsSet;	// Based on user input
}; 


int		ComputeTuring(struct TURING hz);
int 	TuringReport(struct TURING tur);
int 	ComputeTuringK(mpfr_t Gram);

int 	HardyZCallbackA(mpfr_t t, mpfr_t HardyZ, int i, int CallerID);
int 	HardyZCallbackB(mpfr_t t, mpfr_t HardyZ, int i, int CallerID);
int 	HardyZCallbackC(mpfr_t t, mpfr_t HardyZ, int i, int CallerID);

