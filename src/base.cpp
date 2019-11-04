// base.cpp
//
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
// *************************************************************************
// 


#include <iostream>
#include <string>
#include <iomanip>
#include <cstring>
#include <cstdlib>

#include "base.h"
#include "../config.h"

using namespace std;

//************************************************************************//
//**************************** PUBLIC FUNCTIONS **************************//
//************************************************************************//

base::base()
{
 strcpy(datafilename,"");	// Data file name
 noanova=false;			// Do not compute anova
 verbose=false;			// Be verbose
 orthogonal=false;		// Do not use model
 ctrules=false;			// Display Cornfield Tukey Rules
 mtable=false;			// Show means table
 homogeneity=false;             // Show tests of homogeneity 
 #ifdef CGI
 memset(buffer,0,sizeof(buffer));
 #else
 do_debug=false;		// Debug - very verbose mode
 #endif
 mtests=NOMTESTS;		// Show multiple tests
 pretransf=NOTRANSF;		// Pre-transformation
 transf=NOTRANSF;		// Apply transformation to data 
 alpha=0.05;
}

base::~base()
{
}

bool base::be_verbose()
{
 return verbose;
}

#ifndef CGI
bool base::debug()
{
 return do_debug;
}
#endif
  
bool base::do_anova()
{
 return !noanova;
}
  
bool base::show_orthogonal()
{
 return orthogonal;
}

int  base::pretransform()
{
 return pretransf;
}

int  base::transformation()
{
 return transf;
}

int  base::get_mtest()
{
 return mtests;
}

bool base::show_ctrules()
{
 return ctrules;
}
  
bool base::show_mtable()
{
 return mtable;
} 

bool base::show_mtests()
{
 if(mtests==NOMTESTS) return false;
 else return true;
} 

bool base::show_var_tests()
{
 return homogeneity;
}

double base::get_alpha()
{
 return alpha;
}

 
const char *base::data_file_name()
{
 return datafilename;
}

#ifdef CGI
void base::set_option(const char *option, int type)
{
 if(strstr("SHOWORTHO",option)){
  switch(type){
   case 1: orthogonal=true; break;
   default : orthogonal=false; break;
  }
 }
 if(strstr("SHOWHOMO",option)){
  switch(type){
   case 1: homogeneity=true; break;
   default : homogeneity=false; break;
  }
 }
 if(strstr("CTRULES",option)){
  switch(type){
   case 1: ctrules=true; break;
   default : ctrules=false; break;
  }
 }
 if(strstr("VERBOSE",option)){
  switch(type){
   case 1: verbose=true; break;
   default : verbose=false; break;
  }
 }
 if(strstr("MTABLE",option)){
  switch(type){
   case 1: mtable=true; break;
   default : mtable=false; break;
  }
 }
 if(strstr("PRETRANSF",option)){
  switch(type){
   case 1: pretransf=SQRTRANSF; break;
   case 2: pretransf=LOGTRANSF; break;
   case 3: pretransf=LNTRANSF; break;
   case 4: pretransf=ARCSTRANSF; break;
   case 5: pretransf=MULT100; break;
   case 6: pretransf=DIV100; break;
   default : pretransf=NOTRANSF; break;
  }
 }
 if(strstr("POSTRANSF",option)){
  switch(type){
   case 1: transf=SQRTRANSF; break;
   case 2: transf=LOGTRANSF; break;
   case 3: transf=LNTRANSF; break;
   case 4: transf=ARCSTRANSF; break;
   case 5: transf=MULT100; break;
   case 6: transf=DIV100; break;
   default : transf=NOTRANSF; break;
  }
 }
 if(strstr("MTEST",option)){
  switch(type){
   case 0: mtests=NOMTESTS; break;
   case 1: mtests=SNK; break;
   case 2: mtests=TUKEY; break;
   default: mtests=NOMTESTS; break;
  }
 }
}

void base::set_alpha(double a){
 if((a>=0)&&(a<=1.0)) alpha=a;
 else alpha = 0.05;
}

#else
void base::parse_args(int argc, char *argv[])
{
 int 	i;
 bool 	InputFileExists=false;
 char	token[20];
 
 i=1;
 do{
  if(argv[i][0]=='-'){
   switch(argv[i][1]){
    case 'x': mtable=true; i++; break;         // output means table
    case 'n': noanova=true; i++; break;        // do not ouput anova table
    case 'v': verbose=true; i++; break;        // be verbose
    case 'r': ctrules=true; i++; break;        // output Cornfield-Tukey rules table
    case 'd': do_debug=true; i++; break;       // debuggin info
    case 'h': homogeneity=true; i++; break;    // show homogeneity tests
    case 'o': orthogonal=true; i++; break;     // show orthogonal model
    case 'f': i++;                             // data file name
              if((i<argc)&&(argv[i][0]!='-')){
               strcpy(datafilename,argv[i]);
	       InputFileExists=true;            
	       i++;
	      }
	      else noanova=true;  
	      break;
    case 'p': i++;                             // type of pre-transformation
              if((i<argc)&&(argv[i][0]!='-')){
               strcpy(token,argv[i]);
	       if(strcmp(token,"sqrt")==0) pretransf=SQRTRANSF;
	       if(strcmp(token,"log")==0) pretransf=LOGTRANSF;
	       if(strcmp(token,"ln")==0) pretransf=LNTRANSF;
	       if(strcmp(token,"arcsin")==0) pretransf=ARCSTRANSF;
	       if(strcmp(token,"asin")==0) pretransf=ARCSTRANSFRAD;
	       if(strcmp(token,"mult")==0) pretransf=MULT100;
	       if(strcmp(token,"div")==0) pretransf=DIV100;
	       i++;
	      } 
	      else{
	       pretransf=NOTRANSF;
	       i++;
	      }
    case 't': i++;                             // type of transformation
              if((i<argc)&&(argv[i][0]!='-')){
               strcpy(token,argv[i]);
	       if(strcmp(token,"sqrt")==0) transf=SQRTRANSF;
	       if(strcmp(token,"log")==0) transf=LOGTRANSF;
	       if(strcmp(token,"ln")==0) transf=LNTRANSF;
	       if(strcmp(token,"arcsin")==0) transf=ARCSTRANSF;
	       if(strcmp(token,"asin")==0) transf=ARCSTRANSFRAD;
	       if(strcmp(token,"mult")==0) transf=MULT100;
	       if(strcmp(token,"div")==0) transf=DIV100;
	       i++;
	      } 
	      else{
	       transf=NOTRANSF;
	       i++;
	      } 
	      break;
     case 'm': i++;                             // multiple tests
              if((i<argc)&&(argv[i][0]!='-')){
               strcpy(token,argv[i]);
	       if(strcmp(token,"snk")==0) mtests=SNK;
	       if(strcmp(token,"tukey")==0) mtests=TUKEY;
	       i++;
	      } 
	      else{
	       mtests=NOMTESTS;
	       i++;
	      } 
	      break;	      	      
     case 'a': i++;                             // alpha for multiple tests
              if((i<argc)&&(argv[i][0]!='-')){
               strcpy(token,argv[i]);
	       alpha=atof(token);
	       i++;
	      } 
	      else{
	       alpha=0.05;
	       i++;
	      } 
	      break;	      
    default: i++; break;	      
   } 
  }
  else i++;  
 }while(i<argc);
}


void base::help()
{
 cout << "mwanova - Multi-Way Analysis of Variance  Version " << VERSION << endl;
 cout << "Copyright (C)  2013  Antonio Santos (amsantos@fc.up.pt)" << endl << endl;
 cout << "This program is distributed under the GNU General Public License" << endl;
 cout << "See the COPYING file for details" << endl << endl;
 cout << "Program usage:" << endl;
 cout << "  mwanova -f DataFile [options]" << endl;
 cout << "  cat DataFile | mwanova -f [options]" << endl << endl;
 cout << "OPTIONS:" << endl;
 cout << "  -x  output means table      -n  do not output anova table" << endl;
 cout << "  -o  output orthogonal model -h  homogeneity tests" << endl;
 cout << "  -v  be verbose              -r  output Cornfield-Tukey rules table" << endl;
 cout << endl;
 cout << "  -p sqrt|log|ln|arcsin|asin|mult|div  pre-transform data " << endl;
 cout << "  -t sqrt|log|ln|arcsin|asin|mult|div  transform data " << endl;
 cout << "  -m snk|tukey                         multiple comparison tests" << endl;
 cout << "  -a <alpha> [default 0.05]            alpha for multiple tests" << endl;
 cout << endl;
 cout << "Read the man page for more information" << endl;
}
#endif


//------------------------------------------------------------------------//
// Draws a header with an string 'h' in the middle before outputing any   //
// data. The string is enclosed in '---'.                                 //
//------------------------------------------------------------------------//



void base::header(const char *h)
{
 #ifdef CGI
 cout << "<center><h3>" << h << "</h3></center>" << endl;
 #else
 unsigned int rest,i;
 char s[ROWSIZE+1];
 
 memset(s,0,ROWSIZE+1);
 memset(s,'-',ROWSIZE);
 
 if(strlen(h)<ROWSIZE){
  rest=ROWSIZE-strlen(h);
  rest/=2;
  for(i=0;i<strlen(h);i++) s[rest+i]=h[i];
 } 
 cout << s << endl;
 #endif
}

//------------------------------------------------------------------------//
// Draws a footer after outputing data. The footer is a string filled     //
// with '-----------'.                                                    //
//------------------------------------------------------------------------//

void base::footer()
{
 #ifdef CGI
 cout << "<HR>" << endl;
 #else
 char s[ROWSIZE+1];
 
 memset(s,0,ROWSIZE+1);
 memset(s,'-',ROWSIZE);
 cout << s << endl;
 #endif
}



