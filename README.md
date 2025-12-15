# Turing
Using the libHGT static library, calculates values needed to apply Turing's Method

## Overview

We provide here a C program to calculate the values needed to apply Turing's Method -- used to verify the Riemann Hypothesis up to an ordinate 't' on the critical line.
The actual calculation work is done by calls to functions in the **libHGT** static library.  The code here is a front end to that calculating code that: (1) gathers the user-requested parameters (via the command line) for the calculations, (2) validates those command line parameters (by calls to the library), (3) passes those parameters to the calculating code in the library), (4) saves the calculated values in local arrays, and (5) prints (to stdout) a report of those calculations.

## Building the Executable

For Windows 11 users, an executable is included with any release posted on GitHub.

For other operating systems, you will need to build the executable, as follows.

*  You need the [**gcc**][gcc-gnu-link] C compiler installed on your system.  That installation must include the **libmpfr.a** and **libgmp.a** (floating point) static libraries.

*  From [**libHGT**][libhgt-link], you need to: (1) create the **libhgt.a** static library file, and (2) make that library file plus **hgt.h** visible to the **gcc** compiler.

*  Following the build logic in the **maketuring.bat** file, you need to create the necessary 'makefile', in the form that applies to your operating system and the **gcc** compiler.

You can then build the Turing executable from the provided source files.

## Files

This distribution consists of the following files:

  * [README.md][readme-link]. The file you are currently reading.
  
  * [turing.c][turing-c-link]. The entry pont for our program.  This source code file is 
  quite straightforward. We validate the user’s command line input, save their
  choices and then call the 'ComputeTuring' function, which in turn calls the **libHGT** functions that do the actual calculations.
 
  * [CompTuring.c][CompTuring-c-link]. This source code file includes the 'ComputeTuring', 'ComputeTuringK' and 'TuringReport' functions, which: (1) call the **libHGT** functions that do the actual calculations, (2) compute the number of Gram segments needed, and (3) print the output reports.
  
  * [turing.h][turing-h-link]. The is the only (local) include file for the program.  
  
  * [maketuring.bat][maketuring-bat-link]. The is the "makefile" for the program.  Currently,
  this file is a Windows batch file (**not** an actual makefile), but can be easily converted to 
  a standard makefile.

## Command Line Parameters

*  -t [positive number]	Location of t along the critical line - this parameter is required. (Digits and '.' only).
*  -g [positive integer]	Count of the number of Gram intervals to check - between 1 and 16, defaults to 8.
*  -c [positive integer]	Count of the number of Z(t) values to check in a Gram interval - between 8 and 32, defaults to 8.
*  -p [positive integer]	Decimal point digits of 't values to show in report - between 2 and 60, defaults to 6.
*  -b [positive integer]	Floating point bits: 128 <= b <= 1024 - defaults to 256.
*  -d [positive integer]	Used for debugging only.  Please disregard.
*  -k [positive integer]	Number of threads to use - defaults to 1, maximum of 8.
*  -h			Show command line parameters.  All other parameters will be ignored.
*  -s			Report the total seconds taken to compute the Hardy Z values.
*  -v			Verbose report (otherwise CSV only).

## Terms of use

This **Turing Method Calculator** is free and distributed under the
**MIT License** (MIT). 

We used the [**gcc**][gcc-gnu-link] compiler and the [**MPFR**][mpfr-link] floating point library.
We also used the [**msys2**][msys2-link] software distribution and building platform for windows.
See their respective links for theirs terms of license.  

[website-link]:			https://riemann1859.com
[readme-link]:			https://github.com/terry98004/Turing/blob/master/README.md
[turing-c-link]:		https://github.com/terry98004/Turing/blob/master/turing.c
[CompTuring-c-link]:	https://github.com/terry98004/Turing/blob/master/CompTuring.c
[turing-h-link]:		https://github.com/terry98004/Turing/blob/master/turing.h
[maketuring-bat-link]:	https://github.com/terry98004/Turing/blob/master/maketuring.bat
[mpfr-link]:			https://www.mpfr.org/
[gcc-gnu-link]:			https://gcc.gnu.org/
[msys2-link]:			https://www.msys2.org/
[libhgt-link]:			https://github.com/terry98004/libHGT/
