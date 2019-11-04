// mwanova - Multi-Way Analysis of Variance.
// Copyright (C) 2001  Antonio Santos (amsantos@fc.up.pt)
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
//
// Use this file to mofify some of the constants and limits of the program
// 

#ifndef CONF_H
#define CONF_H 1

#ifdef CGI
#define MAXBUFF 100000
#endif

#define MAXFACTORS 10     	// Maximum number of factors allowed           
#define MAXLEVELS  100     	// Maximum number of levels per factor allowed 
#define MAXNAME	   10     	// Maximum name size (of factors or codes) in chars


// COMBINS is equal to 2^MAXFACTORS

#define COMBINS    1024;

// MAXTERMSIZE is equal to COMBINS/2*MAXNAME and is used to
// store variance components of terms after Cornfield Tukey Rules

#define MAXTERMSIZE 5120	

// Definitions for factor types

#define FIXED      0
#define RANDOM     1

// Some typdefs for the main variables. They are hard-coded because I'm not 
// a professional programmer and, although I understand pointers and dynamic 
// variables, I was not in the mood to strugle with new's and delete's. This 
// way, if the main class (which is dynamic) is created, the program will run 
// for sure. Modify them at your own risk... 

typedef   char 	FACTORNAMES[MAXFACTORS][MAXNAME+1];
typedef   char  DATANAME[MAXNAME+1];
typedef   char	LEVELS[MAXFACTORS];
typedef   char  CODES[MAXFACTORS];
typedef   char  CODENAMES[MAXFACTORS][MAXLEVELS][MAXNAME+1];
typedef   char	FACTORTYPES[MAXFACTORS];
typedef   char	NESTED[MAXFACTORS][MAXFACTORS];

// Some defs for the transformations

#define NOTRANSF	0
#define SQRTRANSF	1
#define LNTRANSF	2
#define LOGTRANSF	3
#define ARCSTRANSF	4
#define ARCSTRANSFRAD	5
#define MULT100         6
#define DIV100		7  

// Some defs for the multiple tests

#define NOMTESTS	0
#define SNK		1
#define TUKEY		2

// Precision of output values

#define PRECISION	3

// Sizes of various variables for output

#define SSSIZE		12
#define	DFSIZE 		5
#define	MSSIZE		12
#define	FRSIZE	 	12
#define	PRSIZE	 	8

// Size of an output row in chars (usually a the size of console)

#define ROWSIZE		78

// Define to get more debug info

//#define DEBUG_DATA   	
//#define DEBUG_MODEL   	
//#define DEBUG_NESTING
//#define DEBUG_GET_PARTIAL_SS
//#define DEBUG_BUILD_ORTHOGONAL_MODEL
//#define DEBUG_BUILD_MODEL
//#define DEBUG_CTRULES

#endif /* !CONF_H */
