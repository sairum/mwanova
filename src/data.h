// mwanova - Multi-Way Analysis of Variance.
// Copyright (C) 2001 Antonio Santos (amsantos@fc.up.pt)
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
// This class holds the data summarized as partial sums, sums of squares and
// number of replicates of all combinations of factor levels. It will check
// if the number of replicates is the same for all levels of factor combinations
// and if there are no missing level combinations. In the first case, missing
// values will be replaced by a dummy value (the average of all other values
// in that combination of factors) and the total degrees of freedom will be
// corrected. If a combination of factor levels is missing the analysis stops.
// The class data will also guess the model being tested by finding which
// factors are nested and which are orthogonal. 

#ifndef DATA_H
#define DATA_H 1

#include "conf.h"
#include "base.h"

// Structure that will hold factor level combinations and respective sums, 
// sums of squares and replicates. A list of these structures will be created
// dynamically for all combinations read from the data file.

struct partial{
 CODES  codes,orig;
 double sum;
 double sum2;
 double var;
 int    n;
 partial  *next; 
 partial  *prev;
};

struct combins{
 CODES codes;
 combins *next;
};

struct group{
 double member1;
 double member2;
 group *next;
};
 
// Main class data

class data: public base{
 private:  
  char factors;			// Number of Factors                   
  LEVELS levels;		// Number Levels per Factor  
  LEVELS origlevels;            // ...the same thing: 'levels' will be
                                //    modified if there are nested factors 
  FACTORNAMES factor_name;	// Names of each Factor     
  DATANAME data_name;		// Name of data valriable
  CODENAMES code_name;		// Names of each Level in each Factor
  FACTORTYPES factor_type;	// Type of factor: 0-Fixed, 1-Random                 
  int  nt;			// Total number of data points      
  int  n;			// Number of replicates  
  int  correct_df;                                           
  
  NESTED nested;		//Bears information about nesting
  
  partial *first,*last;	// Pointers to items of list of 'partial' terms 
  int     npartials;    // Number of partial terms
  CODES   code_line;	// Temporary line to store level codes for each observation
  
  // Private functions
  
  char set_code(int, const char *);
  void set_factor_type(int, char);
  void add_code_line(CODES, double);
  void recode(int , int);
  void multi_comp(CODES , int, const char *, double, int, partial *);
  void compare(int , int , double , int , partial *, partial *);  
    
 public:
  data();
  ~data();
  
  // Functions to store data
  
  bool set_factor(const char *);
  void set_factor_name(int, const char *);
  void set_data_name(const char *);
  bool add_code(int, const char *, int);
  bool add_value(int, double, int);
  
  // Functions to get information about the data
  
  int   get_factors();
  int   get_levels(int);
  int   get_correct_df();
  int   get_n();
  int   get_nt();
  
  char  get_factor_type(int);
  const char  *get_factor_name(int);
  const char  *get_code_name(int, int);
  const char  *get_orig_code_name(int , int);
      
  bool  is_nested_into(int, int);
  int   is_nested(int);
  
  double get_partial_SS(CODES);
  double get_CT();	
  double get_sum_of_squares();
  double get_error_ss();
  int    get_error_df();
  double get_total_ss();
  int    get_total_df();  
  int    get_df(CODES);  
  
  void   equalize();
  void   compute_nesting();
  bool   orthogonalize();  
  void   test_homogeneity();
  double transform(double);
  #ifndef CGI
  bool   read_data();
  #else
  bool   read_data(char *);
  #endif
  void   summary();
  void   get_averages(CODES, const char *, double, int, double);
};

#endif /* !DATA_H */
