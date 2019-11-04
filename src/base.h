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
// This is the base class of program mwanova. It provides some logical
// variables to control program flow as well as a function to parse
// command line arguments. All classes inherit this class and share its
// proprerties.

#ifndef BASE_H
#define BASE_H 1

#include "conf.h"

class base{
 private:
  char datafilename[200];
  #ifdef CGI
  char buffer[MAXBUFF];
  #endif
  
  bool noanova;
  bool verbose;  
  bool orthogonal;
  bool ctrules;
  bool mtable;
  bool homogeneity;
  #ifndef CGI
  bool do_debug;
  #endif
  
  int  transf;
  int  pretransf;  
  int  mtests;
  
  double alpha;
   
 public:
  base();
  ~base();
  
  #ifndef CGI
  bool debug();
  #endif
  
  bool be_verbose();  
  bool do_anova();
  bool show_orthogonal();
  int  pretransform();
  int  transformation();
  int  get_mtest();
  bool show_ctrules();
  bool show_mtable();  
  bool show_mtests();
  bool show_var_tests();
  
  
  const char *data_file_name();
  
  #ifndef CGI
  void parse_args(int argc, char *argv[]);
  void help();
  #else
  void set_option(const char *, int );
  void set_alpha(double);
  #endif
  
  double get_alpha();
  
  void header(const char *);
  void footer();
};

#endif /* !BASE_H */
