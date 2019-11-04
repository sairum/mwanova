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

#ifndef MODEL_H
#define MODEL_H 1

#include "data.h"

struct term{
 CODES  fcode;
 int    ctrow[MAXFACTORS+1];
 double ss;
 double SS;
 double MS;
 int    df;
 char   name[100];
 char   against[100];
 double contrast;
 double F;
 int    df2;
 char   *vars;
 term   *next;
 term   *prev;
};

class model: public data{
 private:
  CODES  fcode;
  char	 fname[200];
  term   *first,*last;
  char   term_name[200];
  
  void   insert_term(int);
  int    get_order(CODES);
  bool   is_included(CODES, CODES);
  double get_SS(term *);
  double get_error_ms();
  void   build_orthogonal_model();  
  void   build_model();
  int    table_entry(int, CODES);
  bool   is_component(CODES, CODES);
  void   ctrules();
    
  const char *set_term_name(CODES, char *);  
  
  void   averages();
  
 public:
  model();
  ~model();
  
  void write_anova();
  void run();  
};

#endif /* !MODEL_H */