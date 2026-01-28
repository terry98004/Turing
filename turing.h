// -------------------------------------------------------------------
// Program last modified January 28, 2026. 
// Copyright (c) 2025-2026 Terrence P. Murphy
// MIT License -- see turing.c for details.
// -------------------------------------------------------------------

#define	TUR_HARDY_WIDTH		15
#define	TUR_HARDY_DECIMALS	10

struct TURING {
	char	tBuf[HGT_MAX_CMDLINE_STRLEN + 2];	// holds entered '-t' value
	int 	CountZ;    	 			// number of Z(t) in Gram interval
	int 	CountGram;   			// number of Gram interval	
	bool	Verbose;				// T/F: verbose report
	bool	ShowSeconds;			// T/F: report compute second taken
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
	bool		Good;			// T/F: (-1)^n * HardyZ[0] > 0
	bool		OddZeros;		// T/F: this and next Gram point are both good or both bad
	int			ZerosFound;		// Number of zero crossings between intervals
}; 

struct HARDYINFO {
	double		hzValue;
	double		hzRise;			// Rise (positive) or Fall (negative) of hzValue in interval 
	bool		ZeroCross;		// T/F: zero crossed between previous and current Hardy points
	bool		TowardZero;		// T/F: Value is moving toward zero (no slope is towards)
	bool		Lehmer;			// T/F: possible Lehmer condition 
}; 

int		ComputeTuring(struct TURING hz);
int 	TuringReport(struct TURING tur);
int 	ComputeTuringK(mpfr_t Gram);

int 	HardyZCallbackA(mpfr_t t, mpfr_t HardyZ, int i, int CallerID);
int 	HardyZCallbackB(mpfr_t t, mpfr_t HardyZ, int i, int CallerID);

