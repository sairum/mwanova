// main.cpp
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
// 

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "model.h"
using namespace std;

int main(int argc, char *argv[])
{
 model *d;                                
 #ifdef CGI
 char *buff,*method,*postdata;
 const char *len; 
 long contentlength;
 int check_read = 0;
 #endif 
 
 #ifndef CGI
 d = new model;
 if(d){
  if(argc>1){ 
   d->parse_args(argc,argv);
   if(d->read_data()) d->run();
   
  }else d->help();
  delete d; 
 }
 else cout << "Cannot allocate memory for MWANOVA!... Exiting" << endl;
 #else
 cout << "Content-type: text/html\r\n\r\n";
 cout << "<!DOCTYPE html>\r\n\r\n";
 cout << "<HTML>\r\n\r\n<BODY>\r\n\r\n";
 
 method = getenv("REQUEST_METHOD");
 if(!method){
  fprintf(stdout, "unknown request method\n");
  exit(EXIT_FAILURE);
 }
 if(strcmp(method,"POST") != 0) cout << "This program only accepts POST data" << endl;
 else{
  len = getenv("CONTENT_LENGTH");
  if(!len) contentlength = 0;
  else  contentlength = atoi(len);
  if(contentlength == 0)  cout << "No data was received." << endl;
  else{
   postdata = new char[contentlength  + 1]; 
   if(postdata){
    /* read in POST data */
    check_read = fread(postdata, sizeof(char), contentlength, stdin);
    if(check_read != contentlength){
     cout << "Error reading input" << endl;
     exit(EXIT_FAILURE);
    }
    postdata[contentlength] = '\0';    
    d = new model();
    if(d){
     if(d->read_data(postdata)) d->run();     
     delete [] postdata;  
     delete d; 
    }
   }
  } 
 }
 cout << "</BODY>\r\n\r\n</HTML>\r\n\r\n";
 #endif
 return 0;
}
