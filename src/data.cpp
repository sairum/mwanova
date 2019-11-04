// data.cpp
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
#include <fstream>
#ifdef CGI
#include <sstream>
#endif
#include <cstdlib>
#include <iomanip>
#include <cstring>
#include <cmath>
#include "data.h"
#include "probs.h"
using namespace std;

//************************************************************************//
//**************************** PRVATE STUFF ******************************//
//************************************************************************//

//------------------------------------------------------------------------//
// This function sets the code for a new level 'cname' of factor 'factnum'// 
// First it tests if factor 'factnum' has level 'cname'. If so it returns //
// the char code corresponding to the level. If there is no level with    //
// name 'cname' it adds another level to factor 'ffactnum' and increases  //
// 'levels[]' accordingly, returning the new char code                    //
//------------------------------------------------------------------------//

char data::set_code(int factnum, const char *cname)
{
 int i;
 if((factnum>=0)&&(factnum<factors)){
  for(i=0;i<levels[factnum];i++){
   if(strcmp(code_name[factnum][i],cname)==0) return (char) i;
  } 
  i=levels[factnum];
  levels[factnum]++;
  origlevels[factnum]++;
  strcpy(code_name[factnum][i],cname);
  return (char) i; 
 }
 else return 0; 
}

//------------------------------------------------------------------------//
// This function sets the factor type to 'ftype' (0 - FIXED, 1 - RANDOM)  //
//------------------------------------------------------------------------//

void data::set_factor_type(int factnum, char ftype)
{
 if((factnum>=0)&&(factnum<factors)){
  factor_type[factnum]=ftype;
 }
}

//------------------------------------------------------------------------//
// This function adds an observation to a list of partials. Each partial  //
// has a code line which bears the codes of the levels of each factor in  //
// the analysis. 'cline' is compared with the code lines of all items in  // 
// the list to check if a value should be added to an existent partial or //
// if a new partial should be created. If an item with a similar 'cline'  //
// is present, the value 'val' is added to the sum and sum of squares,    //
// and the 'n' is updated, otherwise a new item is created. Items are     //
// inserted in a sorted way to facilitate the subsequent computations.    //
//------------------------------------------------------------------------//

void data::add_code_line(CODES cline, double val)
{
 bool equal,inserted;
 partial *t1, *t2;
 
 if(!first){	
 		
  // First item is null, so create it! 
    
  first = new partial;
  npartials++;
  memcpy(first->codes,cline,MAXFACTORS);
  memcpy(first->orig,cline,MAXFACTORS);
  first->sum=val;
  first->sum2=pow(val,2);
  first->n=1;
  first->next=NULL;
  last=first;  
 }
 else{				
  
  // Items exist, so compare them with 'cline'
  
  t1=first; 
  
  do{
   equal=true;
   if(memcmp(cline,t1->codes,MAXFACTORS)!=0) equal=false;
   if(!equal) t1=t1->next;
  }while((t1!=NULL)&&(!equal));
  
  if(equal){			
   
   // There is an item similar to 'cline'... add the data
   
   t1->sum+=val;
   t1->sum2+=pow(val,2);
   t1->n++;
   
  }
  else{
  	
   // 'cline' is a new item, create and insert it!
  
   t1=first;
   if(memcmp(cline,t1->codes,MAXFACTORS)<0){   	
   
    // 't1' is inserted before 'first'
    
    t1 = new partial;
    npartials++;
    t1->next=first;
    memcpy(t1->codes,cline,MAXFACTORS);
    memcpy(t1->orig,cline,MAXFACTORS);
    t1->sum=val;
    t1->sum2=pow(val,2);
    t1->n=1;
    first=t1;
   }
   else{
    inserted=false;
    while((t1->next!=NULL)&&(!inserted)){	
    
     // 't1' is inserted into the list
     
     if(memcmp(cline,t1->next->codes,MAXFACTORS)<0){
      t2 = new partial;
      npartials++;
      memcpy(t2->codes,cline,MAXFACTORS);
      memcpy(t2->orig,cline,MAXFACTORS);
      t2->sum=val;
      t2->sum2=pow(val,2);
      t2->n=1;
      t2->next=t1->next;
      t1->next=t2;
      inserted=true;
     }
     t1=t1->next;
    }
    if(!inserted){				
    
     // 't1' is inserted at the end of the list
     
     t1=last;
     last = new partial;
     npartials++;
     memcpy(last->codes,cline,MAXFACTORS);
     memcpy(last->orig,cline,MAXFACTORS);
     last->sum=val;
     last->sum2=pow(val,2);
     last->n=1;
     last->next=NULL;
     t1->next=last;
    } 
   }
  } 
 }
}


//****************************************************************************//
//**************************** PUBLIC FUNCTIONS ******************************//
//****************************************************************************//

data::data()
{
 first=NULL;
 last =NULL;
 
 factors=0;
 n=0;
 nt=0;
 correct_df=0;
 npartials=0;
 
 memset(levels,0,sizeof(levels));
 memset(factor_name,0,sizeof(factor_name));
 memset(factor_type,0,sizeof(factor_type));
 memset(code_name,0,sizeof(code_name));
 memset(code_line,0,sizeof(code_line));
 memset(nested,0,sizeof(nested));
}

//------------------------------------------------------------------------//
// data destructor. Clears up the list of factor level combinations     //
//------------------------------------------------------------------------//

data::~data()
{
 partial *t;
 if(first){
  do{
   t=first->next;
   delete first;
   first=t;
  }while(first); 
 }
 #ifdef DEBUG_DATA
 cout << "Destructing 'data' variable" << endl;
 #endif
}

//------------------------------------------------------------------------//
// This function sets a new factor, increasing the number of factors and  //
// storing its name in 'factor_name'. It also checks if the factor is a   //
// fixed or a random factor (a factor name with a '*' at the end), in     //
// which case it removes the '*' from the name.                           //
//------------------------------------------------------------------------//

bool data::set_factor(const char *fname)
{
 char fn[100];
 
 if((factors>=0)&&(factors<MAXFACTORS)){
  strncpy(fn,fname,100);
  if(strchr(fn,'*')!=NULL){
   factor_type[factors]=RANDOM;
   fn[strlen(fname)-1]=0;
   strncpy(factor_name[factors],fn,MAXNAME);
  }
  else{
   factor_type[factors]=FIXED;
   strncpy(factor_name[factors],fn,MAXNAME);
  } 
  factors++;
 }
 else{
  #ifdef CGI
  cout << "Number of factors exceeds MAXFACTORS (" << MAXFACTORS << ")<p>" << endl;
  #else 
  cerr << "Number of factors exceeds MAXFACTORS (" << MAXFACTORS << ")" << endl;
  #endif
  return false;
 }
 return true;
}

//------------------------------------------------------------------------//
// This functions sets the name of factor 'fact' to 'fname'               //
//------------------------------------------------------------------------//

void data::set_factor_name(int fact, const char *fname)
{
 if((fact>=0)&&(fact<MAXFACTORS)){
  strncpy(factor_name[fact],fname,MAXNAME);
 } 
}


//------------------------------------------------------------------------//
// This functions sets the name of the data vector                        //
//------------------------------------------------------------------------//

void data::set_data_name(const char *dname)
{
 strncpy(data_name,dname,MAXNAME);
}

//------------------------------------------------------------------------//
// This function adds a level code (that has been read from the datafile) //
// to the temporary variable 'code_line', inserting it in the position    //
// corresponding to the factor number. Code names are converted to a      //
// single and unique char (thus, no more than 250 levels can be           //
// represented, but that is fair enough...)                               //
//------------------------------------------------------------------------//

bool data::add_code(int factnum, const char *cname, int l)
{
 char cn[100];
 
 strncpy(cn,cname,100);
 if(strlen(cn)>(MAXNAME)) strncpy(cn,cn,MAXNAME); // Trim long code names 
 if((factnum>=0)&&(factnum<factors)){
  code_line[factnum]=set_code(factnum, cn);
 }
 else{
  #ifdef CGI
  cout << "Number of codes exceeds number of factors (";
  cout << get_factors() << ") in line " << l+1 << "<p>" << endl;
  #else
  cerr << "Number of codes exceeds number of factors (";
  cerr << get_factors() << ") in line " << l+1 << endl;
  #endif
  return false;
 }
 return true;
}

//------------------------------------------------------------------------//
// This function tests if the number of codes in 'code_line' are equal to //          
// the number of factors. If so, it calls add_code_line to add this new   //
// line of codes and the respective data value to a list of factor level  //
// combinations. The variable 'fact' holds the number of codes read before//
// the value 'val' and is set in the function anova::read_data().         //
//------------------------------------------------------------------------//

bool data::add_value(int fact, double val, int l)
{
 // Test if the number of codes is equal to the number of factors
 
 if(fact<get_factors()){
  #ifdef CGI
  cout << "Number of codes in line " << l+1;
  cout << " is less than the number of factors (";
  cout << get_factors() << ")" << "<p>" << endl;
  #else
  cerr << "Number of codes in line " << l+1;
  cerr << " is less than the number of factors (";
  cerr << get_factors() << ")" << endl;
  #endif
  return false;
 } 
 add_code_line(code_line,val);		 // Add 'code_line' to the list
 memset(code_line,0,sizeof(code_line));  // Clear 'code_line'
 nt++;					 // Increment number of replicates
 return true;
}

//------------------------------------------------------------------------//
// This function returns the number of factors in the analysis            //
//------------------------------------------------------------------------//

int data::get_factors()
{
 return (int) factors;
}


//------------------------------------------------------------------------//
// This functions gets the dfs to be subtracted to Error DF (missing data)//
//------------------------------------------------------------------------//

int data::get_correct_df()
{
 return correct_df;
}

//------------------------------------------------------------------------//
// This functions returns the number of replicates per sample             //
//------------------------------------------------------------------------//

int data::get_n()
{
 return n;
}


//------------------------------------------------------------------------//
// This functions returns the total number of replicates                  //
//------------------------------------------------------------------------//

int data::get_nt()
{
 return nt;
}

//------------------------------------------------------------------------//
// This function returns the levels of factor 'factnum'. Levels are       //
// stored in 'levels[]'                                                   //
//------------------------------------------------------------------------//

int data::get_levels(int factnum)
{
 if((factnum>=0)&&(factnum<factors)) return (int) levels[factnum];
 return 0;
}


//------------------------------------------------------------------------//
// This function returns the type of factor 'factnum'. 0 means FIXED and  //
// 1 means RANDOM. Factor types are stored in 'factor_type[]'             //
//------------------------------------------------------------------------//

char data::get_factor_type(int factnum)
{
 if((factnum>=0)&&(factnum<factors)) return factor_type[factnum];
 return 0;
}

//------------------------------------------------------------------------//
// This function returns the name of factor 'factnum'                     //
//------------------------------------------------------------------------//

const char *data::get_factor_name(int factnum)
{
 if((factnum>=0)&&(factnum<factors)) return factor_name[factnum];
 return "";
}

//------------------------------------------------------------------------//
// This function returns the code name for level 'levnum' of factor       //
// 'factnum'. Codes are stored in 'code_name[]'                           //
//------------------------------------------------------------------------//

const char *data::get_code_name(int factnum, int levnum)
{
 if((factnum>=0)&&(factnum<factors)){
  if((levnum>=0)&&(levnum<levels[factnum])){  
   return code_name[factnum][levnum];
  }
 }
 return "";
}

//------------------------------------------------------------------------//
// This function returns the original code name for level 'levnum' of     //
// factor 'factnum'. Remember that 'levels[]' will be modified if there   //
// are nested factors, after orthogonalization.                           //
//------------------------------------------------------------------------//

const char *data::get_orig_code_name(int factnum, int levnum)
{
 if((factnum>=0)&&(factnum<factors)){
  if((levnum>=0)&&(levnum<origlevels[factnum])){  
   return code_name[factnum][levnum];
  }
 }
 return "";
}

//------------------------------------------------------------------------//
// This function returns true if 'fact1' is nested in 'fact2'             //
//------------------------------------------------------------------------//

bool data::is_nested_into(int fact1, int fact2)
{
 if(nested[fact1][fact2]==1) return true;
 else return false;
}

//------------------------------------------------------------------------//
// This function returns the number of factors in which 'fact1' is nested //
//------------------------------------------------------------------------//

int data::is_nested(int fact1)
{
 int i,f;
 f=0;
 for(i=0;i<MAXFACTORS;i++) f+=(int) nested[fact1][i];
 return f;
 //if(memchr(nested[fact1],1,MAXFACTORS)!=NULL) return true;
 //else return false;
}

//------------------------------------------------------------------------//
// This function returns the partial sum of squares for a factor or a     //
// combination of factors. Variable 'cline' is an array with 1s in each   //
// cell that corresponds to a factor being analyzed. Thus [1001000] means //
// that factors 0 and 3 are being considered and the partial SS is an     //
// interaction SS. The list of partials is searched for different         //
// combinations of levels of the factors involved. In each combination a  //
// cumulative sum of values is stored as well as the cumulative sum of    //
// degrees of freedom. In the end, the individual sums are squared and    //
// summed, and finally divided by the cumulative degrees of freedom of the//
// combinations (the df of the first combination is used, but they should //
// be the same for all combinations). During the process a dynamic list of//
// combinations is built and will be discarded in the end.                //
//------------------------------------------------------------------------//

double data::get_partial_SS(CODES cline) 
{
 struct part{
  CODES codes;
  double ss;
  double n;
  part *next;
 };
 
 double  ss,n;
 part    *firstss,*lastss,*s;
 partial *t;
 bool    equal;
 int     i;
 
 firstss=NULL;
 lastss=NULL;
 ss=0;
 n=0;
    
 // Go along list of partials partials
 
 if(first){
  t=first;
  do{
   if(!firstss){     
     
    // No items in the list of partial SS: create!
     
    firstss = new part;
    if(firstss){
     for(i=0;i<get_factors();i++){
      if(cline[i]>0) firstss->codes[i]=t->codes[i];
      else firstss->codes[i]=0;
     }
     firstss->ss=t->sum;
     firstss->n=t->n;
     firstss->next=NULL;
     lastss=firstss;
    }
   }
   else{  
     
    // There are already items in the list of partial SS. 
    // Find out if current 't->codes' is there...
    
    s=firstss;
    do{
     equal=true;
     for(i=0;i<get_factors();i++){
      if(cline[i]>0) if((s->codes[i])!=(t->codes[i])) equal=false;
     }
     if(!equal) s=s->next;
    }while((s)&&(!equal));
     
    if(equal){
     
     // There is one item equal to 't->codes'! Add ss and n values to it
     
     s->ss+=t->sum;
     s->n+=t->n;
    }
    else{
     
     // No item in the list is equal to current 't->codes'. 
     //Add a new one.
     
     s = new part;
     if(s){
      for(i=0;i<get_factors();i++){
       if(cline[i]>0) s->codes[i]=t->codes[i];
       else s->codes[i]=0;
      }
      s->ss=t->sum;
      s->n=t->n;
      s->next=NULL;
      lastss->next=s;
      lastss=s;
     }
    }
   }    
   t=t->next;
  }while(t);
  
  // Square all partials and sum them up, divide the total by the 
  // degrees of freedom of one of the terms (they should be all
  // equal) and destruct the list before returning SS!
  
  if(firstss){
   s=firstss;
   do{
    ss+=pow(s->ss,2);
    s=s->next;
   }while(s);
   
   #ifdef DEBUG_GET_PARTIAL_SS
   for(i=0;i<get_factors();i++) cout << (int) cline[i];
   cout << endl;   
   s=firstss;
   do{
    for(i=0;i<get_factors();i++) cout << (int) s->codes[i];
    cout << ": " << s->ss << endl;
    s=s->next;
   }while(s);
   cout << endl;
   #endif
   
   // Make sure 'n' is bigger than 0
   
   if((firstss->n)>0) ss/=firstss->n;
   else ss=0;
   
   // Destruct list of 'part's
   do{
    s=firstss;
    firstss=s->next;
    delete s;
   }while(firstss);
  }
 }
 return ss;
}

//------------------------------------------------------------------------//
// Compute Correction Term 'CT' which is the squared sum of all values    //
// divided by the total number of observations. This term is essencial    //
// for the computation of partial sums of squares (SS)                    //
//------------------------------------------------------------------------//

double data::get_CT()
{
 partial *t;
 double  sum;
 
 // Make sure a list of partials exists and
 // nt > 0
 
 if(first&&(get_nt()>0)){
  t=first;
  sum=0;
  do{
   sum+=t->sum; 
   t=t->next;
  }while(t);
  
  sum=pow(sum,2);
  sum/=(double)nt;
  return sum;
 }
 else return 0;
}

//------------------------------------------------------------------------//
// This function returns the sum of all squared values. This value is     //
// essential to compute Error or Residual SS.                             //
//------------------------------------------------------------------------//

double data::get_sum_of_squares()
{
 double  sum=0;
 partial *t;
 
 if(first){
  t=first;
  do{
   sum+=t->sum2;
   t=t->next;
  }while(t);
 }
 return sum;
}

//------------------------------------------------------------------------//
// This function returns the sum of squares of the Error or Residual      //    
//------------------------------------------------------------------------//

double data::get_error_ss()
{
 double  sum=0;
 partial *t;
 
 if(first){
  t=first;
  do{
   sum+=pow(t->sum,2);
   t=t->next;
  }while(t);
 }
 #ifdef DEBUG_DATA
 cerr << "Error = " << get_sum_of_squares()-sum/get_error_df() << endl;
 #endif
 return get_sum_of_squares()-sum/get_n();
}

//------------------------------------------------------------------------//
// This function returns the df of the Error or Residual                  //
//------------------------------------------------------------------------//

int data::get_error_df()
{
 int i,df;
 df=1;
 for(i=0;i<get_factors();i++) df*=get_levels(i);
 df*=(get_n()-1);
 df-=correct_df;  // Correct df of error if missing values were inserted
 return df;
}

//------------------------------------------------------------------------//
// This function returns the Total Sum of Sqaures                         //
//------------------------------------------------------------------------//

double data::get_total_ss()
{
 return get_sum_of_squares()-get_CT();
}

//------------------------------------------------------------------------//
// This function returns the dfs of the Total Sum of Squares              //
//------------------------------------------------------------------------//

int data::get_total_df()
{
 return get_nt()-1-correct_df;
}

//------------------------------------------------------------------------//
// This function computes the degrees of freedom of a factor or a         //
// combination of factors (interaction) listed in 'cline'                 //
//------------------------------------------------------------------------//

int data::get_df(CODES cline)
{
 int i,total;
 
 total=1;
 for(i=0;i<get_factors();i++) if(cline[i]>0) total*=(get_levels(i)-1); 
 if(total>0) return total;
 else return 0;
}

//------------------------------------------------------------------------//
// This function tests if all combinations of factor level codes have the //
// same number of replicates. If not, it replaces missing values with the //
// average of values for that particular combinations, updating the       //
// variable 'correct_df' which will be the number of degrees of freedom   //
// to remove from the Error DF in the end of the analysis.                //
//------------------------------------------------------------------------//

void data::equalize()
{
 partial *t;
 int     maxreps=0;   
 int     newreps,i;
 double  average;
 
 // Find out the largest set of replicates
 
 if(first){
  t=first;
  do{
   if((t->n)>maxreps) maxreps=t->n;
   t=t->next;
  }while(t);
  
  // Number of replicates per combination of factor in the analysis 
  // should be equal to the largest number of replicates in a combination 
  // of factor levels, since all other combinations will be 'equalized' 
  // according to that value.
  
  n=maxreps;
  
  // Compute averages for sets with less than 'maxreps' observations and
  // insert as many values as needed, updating 'correct_df'. This value 
  // will be subtracted from Error DF in the final analysis...
  
  t=first;
  do{
   if(((t->n)<maxreps)&&((t->n)>0)){
    newreps=maxreps-(t->n);
    average=(t->sum)/(t->n);
    for(i=0;i<newreps;i++){
     t->sum+=average;
     t->sum2+=pow(average,2);
     t->n++;
     correct_df++;
     nt++;
    } 
   } 
   t=t->next;
  }while(t);
 }
}


//------------------------------------------------------------------------//
// This function tests if factor1 is nested in factor2, updating that     //
// information in a matrix variable 'nested'. This matrix has 0 for each  //
// combination of row/col if factor[row] is orthogonal to factor[col],    //
// and 1 if factor[row] is nested in factor[col].                         //
// The algorithm to find if a factor is nested in another is complex.     //
// The goal is to find whow many combinations of levels of the two        //
// factors exist in the list 'partial'. If the combinations are equal to  //
// the product of the levels of both factors, then the factors are        //
// orthogonal. If the number of combinations is less than the product of  //
// the levels of the two factors, the one with more levels is the nested  //
// factor. This is only true if both factors are not nested into a higher //
// factor or gactors. If so, they might be orthogonal, even if expected   //
// combinations are higher than the real number of combinations. So,      //
// in a first phase all factores are contrasted with all others according //
// to the above method. In a second phase, a check is made to see if      //
// factors were other factors are nested are themselves nested into       //
// another factor. If a factor and a nested factor are nested into a      //
// higher factor the expected number of combinations must be divided by   //
// the number of levels of the higher factor. As a general rule, their    //
// expected combinations should be computed as the product of their levels//
// multiplied by the levels of factors where both are nested.             //
//------------------------------------------------------------------------//

void data::compute_nesting()
{
 int     i,j,l,expected_combins,corrected;
 partial *t,*u;
 char    tc1,tc2,uc1,uc2;
 bool    exists;
 int     combins[MAXFACTORS][MAXFACTORS];
 
 memset(combins,0,sizeof(combins));
 
 // Test if this is a multiway anova and a 'partial' list exists...
 
 if((get_factors()>1)&&(first)){  
  
  #ifdef DEBUG_NESTING
  cout << "DEBUG compute_nesting(): Testing for nested factors" << endl;
  #endif
 
  // Compute the number of combinations of levels for each pair of factors
  
  for(i=0;i<(get_factors()-1);i++){
   for(j=(i+1);j<get_factors();j++){
    //combins=0; 
    t=first;
    do{
     if(t==first){
      combins[i][j]++;  // There is at least this combination
      combins[j][i]++;
     } 
     else{
      exists=false;
      u=first;
      tc1=t->codes[i];
      tc2=t->codes[j];
      do{
       uc1=u->codes[i];
       uc2=u->codes[j];
       if((tc1==uc1)&&(tc2==uc2)) exists=true;	
       u=u->next;
      }while((u!=t)&&(!exists));
      if(!exists){
       combins[i][j]++;
       combins[j][i]++;
      } 
     }
     t=t->next;
    }while(t);   
   }
  }
  
  // Check if number of expected combinations is equal to 
  // the number of observed combinations... if not, assume that
  // the factor with the highest number of levels is the nested
  // one. If both factors have the same number of levels they
  // cannot be nested into each other and might be orthogonal
  // factors nested into a higher factor.
  
  for(i=0;i<(get_factors()-1);i++){
   for(j=(i+1);j<get_factors();j++){
    expected_combins=get_levels(i)*get_levels(j);
    
    #ifdef DEBUG_NESTING
    cout << "Combinations between factor " << i << " and  factor ";
    cout << j << ": " << combins[i][j] << "\t(Expected: ";
    cout << expected_combins << ")" <<  endl;
    #endif
     
    // There is a nested factor  
     
    if(combins[i][j]<expected_combins){   		         
     if(get_levels(i)>get_levels(j)) nested[i][j]=1;
     if(get_levels(i)<get_levels(j)) nested[j][i]=1;
    } 
   } 
  }
  
  // Now check for double nesting...
  
  #ifdef DEBUG_NESTING
  cout << "Retesting for double nesting... " << endl;
  #endif
  
  for(i=0;i<(get_factors()-1);i++){
   for(j=(i+1);j<get_factors();j++){
    corrected=1;
    for(l=0;l<get_factors();l++){
     if((l!=i)&&(l!=j)){
      if(is_nested_into(i,l)&&is_nested_into(j,l)){
       corrected*=get_levels(l);
      }
     }
    }
    expected_combins=get_levels(i)*get_levels(j);
    if((combins[i][j]*corrected)==expected_combins){ 
     nested[i][j]=0;
     nested[j][i]=0;
    }
    else{
     if(get_levels(i)>get_levels(j)) nested[i][j]=1;
     if(get_levels(i)<get_levels(j)) nested[j][i]=1;
    }
   }
  } 
    
  // Finally check if all nested factors are random and
  // if not change their types
  
  for(i=0;i<get_factors();i++){
   if(is_nested(i)){
    if(get_factor_type(i)==FIXED) set_factor_type(i,RANDOM);
   }
  }
 }
}

//------------------------------------------------------------------------//
// This function finds ou all combinations of levels of factors in which  //
// factor 'f' is nested and, for each combination, it recodes the level   //
// codes of factor 'f' so as to orthogonalize it.                         //
//------------------------------------------------------------------------//

void data::recode(int f, int lev)
{
 int i,j;
 partial *t;
 div_t r;
 
 if((f>0)&&(f<get_factors())){
  if(first){
   t=first;
   j=0;
   do{
    j=(int) t->codes[f];
    if(j>=lev){
     r=div(j,lev);
     j=r.rem;
    } 
    t->codes[f]=(char)j;
    t=t->next;
   }while(t); 
  }
 }
}
 
//------------------------------------------------------------------------//
// This function transforms nested factors in the data into orthogonal    //
// factors. Computation of SS and DF is easy for orthogonal data sets.    //
// The SS of an orthogonal data set can then be combined to compute any   //
// other ANOVA model. This function also searches if there are missing    //
// factor level combinations, in which case the analysis stops. When all  //
// factors have been orthogonalized, the total number of combinations in  //
// the list 'partials' must be equal to the product of all factor levels. //
// The variable 'orig' in each partial stores the original levels codes   //
// because in the end it is essential to compute averages of levels and   //
// to assign them the correct names.                                      //
//------------------------------------------------------------------------//

bool data::orthogonalize()
{
 int     i,j,k,nestlev,lev,combins;
 int 	 c,nf;
 int     nested_factors[MAXFACTORS];
 partial *t;
 
 
 memset(nested_factors,0,sizeof(nested_factors));
 
 if((get_factors()>1)&&(first)){
  
  // Find out factors which are nested in others. First those which
  // are nested in one factor, then those which are nested in two, and so on...
  // and put their codes in nested_factors[] in that order.

  nf=0;
  for(i=1;i<get_factors();i++){      
   for(j=0;j<get_factors();j++){
    if(is_nested(j)==i){
     nested_factors[nf]=j;
     nf++;
    }
   }
  }
  
  // Now, for each 'i' factor in nested_factors[] list
  // find out factors where 'i' is nested into and compute the 
  // product of their levels
  
  for(i=0;i<nf;i++){
   memset(code_line,0,sizeof(code_line)); 
   nestlev=1;
   for(j=0;j<get_factors();j++){
    if(is_nested_into(nested_factors[i],j)) nestlev*=get_levels(j);
   }
    
   // If 'i' was to be orthogonal to the factors where it is
   // nested it should have fewer levels. The number of levels
   // of 'i' should be divided by the product of the levels of
   // the factors where it is nested.
    
   lev=get_levels(nested_factors[i])/nestlev;
    
   // Levels of 'i' which are above 'lev' should be recoded
   // using the remainder of their division by 'lev'.

   recode(nested_factors[i],lev);
   levels[nested_factors[i]]=lev;
  }
   
  // Check if all level combinations exist. This is a final check
  // not particularly related with orthogonalization of factors. 
  // First compute number of combinations in the list of partials
  
  t=first;
  combins=0;
  do{
   combins++;
   t=t->next;
  }while(t);
  
  // Then compute the number of expected combinations in an
  // orthogonal model
  
  lev=1;
  for(i=0;i<get_factors();i++) lev*=get_levels(i);
  
  // If 'lev' is not equal to 'combins' there are missing levels...
  // return false and stop the analysis
  
  if(lev!=combins){
   #ifdef CGI
   cerr << "There are missing combinations of factor levels!<p>" << endl;
   cerr << "Bailing out!<p>" << endl;
   #else
   cerr << "There are missing combinations of factor levels!" << endl;
   cerr << "Bailing out!" << endl;
   #endif
   return false;
  }
 }
 return true;
}


//------------------------------------------------------------------------//
// This function tests whether the data is homoscedastic or not by        //
// computing Cochran's C test (Largest Variance/Sum of Variances), and    //
// Fmax (Larges Variance/Smallest Variance)                               //
//------------------------------------------------------------------------//

void data::test_homogeneity()
{
 int df;
 double varmax,varmin,sumvar,var;
 double b1,b2,bc;
 partial *t;
 
 sumvar=varmax=b1=b2=0;
 df=0;
 if(first){
  t=first;
  varmin=(t->sum2-(pow(t->sum,2)/t->n))/(t->n-1);
  do{
   var=(t->sum2-(pow(t->sum,2)/t->n))/(t->n-1);
   sumvar+=var;
   t->var=var;
   if(var>varmax) varmax=var;
   if(varmin>var) varmin=var;
   df+=(t->n-1);
   b2+=(t->n-1)*log(var);
   t=t->next;
  }while(t);
  
  b1=log(sumvar);
  bc=1+1/(3*(npartials-1))*(npartials/(get_n()-1)*1/(df));
  
  if(show_var_tests()){
   #ifdef CGI
   header("HOMOGENEITY TESTS");
   cout << "Cochran's C test: " << endl;
   cout << "\tC for " << npartials << " means and df=" << get_n()-1 << ": ";
   cout << varmax/sumvar;
   cout << " (<I>p=" << cprob(varmax/sumvar,npartials,get_n()-1) << "</I>)" << endl;
   cout << "<P>" << endl;
   cout << "Hartley's Fmax test:" << endl;
   cout << "\tF for " << npartials << " means and " << get_n()-1;
   cout << " dfs: " << varmax/varmin;
   cout << "<P>" << endl;
   cout << "Bartllet's test:" << endl;
   cout << "\tChi-square for df=" << npartials-1 << ": "<< df*b1-b2;
   cout << " (<I>p=" << chiprob((df*b1-b2),npartials-1) << "</I>)";
   cout << endl << "\tChi-square (corrected) for df=" << npartials-1 << ": "<< (df*b1-b2)/bc;
   cout << " (<I>p=" << chiprob((df*b1-b2)/bc,npartials-1) << "</I>)" << endl;
   cout << "<P>" << endl;
   footer();
   #else
   cout << "---------------- HOMOGENEITY TESTS -------------------" << endl << endl;
   cout << "Cochran's C test: " << endl;
   cout << "\tC for " << npartials << " means and df=" << get_n()-1 << ": ";
   cout << varmax/sumvar;
   cout << " (p=" << cprob(varmax/sumvar,npartials,get_n()-1) << ")" << endl;
   cout << endl;
   cout << "Hartley's Fmax test:" << endl;
   cout << "\tF for " << npartials << " means and " << get_n()-1;
   cout << " dfs: " << varmax/varmin;
   cout << endl << endl;
   cout << "Bartllet's test:" << endl;
   cout << "\tChi-square for df=" << npartials-1 << ": "<< df*b1-b2;
   cout << " (p=" << chiprob((df*b1-b2),npartials-1) << ")";
   cout << endl << "\tChi-square (corrected) for df=" << npartials-1 << ": "<< (df*b1-b2)/bc;
   cout << " (p=" << chiprob((df*b1-b2)/bc,npartials-1) << ")" << endl;
   cout << endl << endl;
   #endif
  }
 } 
}

//------------------------------------------------------------------------//
// This function transforms a data value with a specified transformation  //
//------------------------------------------------------------------------//

double data::transform(double v)
{
 switch(pretransform()){
  case SQRTRANSF: v=sqrt(v); break;
  case LNTRANSF:  v=log(v+1); break;
  case LOGTRANSF: v=log10(v+1); break;
  case ARCSTRANSF:  v=asin(sqrt(v))*90/(M_PI/2); break;
  case ARCSTRANSFRAD:  v=asin(sqrt(v)); break; 
  case MULT100: v*=100; break;
  case DIV100: v/=100; break;  
 }
 switch(transformation()){
  case SQRTRANSF: v=sqrt(v); break;
  case LNTRANSF:  v=log(v+1); break;
  case LOGTRANSF: v=log10(v+1); break;
  case ARCSTRANSF:  v=asin(sqrt(v))*90/(M_PI/2); break;
  case ARCSTRANSFRAD:  v=asin(sqrt(v)); break;   
  case MULT100: v*=100; break;
  case DIV100: v/=100; break;  
 }
 return v;
}

#ifndef CGI
//------------------------------------------------------------------------//
// This function reads a anova file. anova files should be in columnar    //
// format, with the last column being the data values. The first line     //
// have the names of the factors, each one with an optional '*' character //
// in the end to be treated as a random factor. Comment lines starting    //
// with an '#' are ignored.                                               //
//------------------------------------------------------------------------//

bool data::read_data()
{
 int		fact,lines;
 ifstream 	df;
 char		string[1000];
 char		*token, temp[100];
 bool		header=true;
 char		delimiters[] = " \t:;,";
 double		v;
 
 lines=0;
 if(strlen(data_file_name())>0){
  df.open(data_file_name());
  if(df){				// Read from data file
   df.getline(string,1000);
   while(!df.eof()){      
    if((strchr(string,'#')==NULL)&&(strlen(string)>0)){
     if(header){			
      token=strtok(string,delimiters);
      do{
       strcpy(temp,token);
       token=strtok(NULL,delimiters);
       if(token!=NULL) if(!set_factor(temp)) return false;
       else set_data_name(temp);
      }while(token!=NULL);
      header=false;				      
     }
     else{
      fact=0;
      token=strtok(string,delimiters);
      do{
       strcpy(temp,token);
       token=strtok(NULL,delimiters);
       if(token!=NULL){
        if(!add_code(fact,temp,lines)) return false; 
	else fact++;
       }
       else{
        v=transform((double) atof(temp));
        if(!add_value(fact,v,lines)) return false;	
       }
      }while(token!=NULL);
     }                
    }
    lines++;
    df.getline(string,1000);
   }                                        
   df.close();
  }
  else {
   cerr << "Error while opening " << data_file_name() << "!... Exiting..." << endl;
   return false;
  }
 }
 else{ 			            // Read data from standard input
  cin.getline(string,1000);
  while(!cin.eof()){      
   if((strchr(string,'#')==NULL)&&(strlen(string)>0)){
    if(header){			
     token=strtok(string,delimiters);
     do{
      strcpy(temp,token);
      token=strtok(NULL,delimiters);
      if(token!=NULL) if(!set_factor(temp)) return false;
      else set_data_name(temp);
     }while(token!=NULL);
     header=false;				      
    }
    else{
     fact=0;
     token=strtok(string,delimiters);
     do{
      strcpy(temp,token);
      token=strtok(NULL,delimiters);
      if(token!=NULL){
       if(!add_code(fact,temp,lines)) return false; 
       else fact++;
      }
      else{
       v=transform((double) atof(temp));
       if(!add_value(fact,v,lines)) return false;	
      }
     }while(token!=NULL);
    }                
   }
   lines++;
   cin.getline(string,1000);
  }
 } 
 return true; 
}
#else

bool data::read_data(char *buff)
{
 istringstream ins(buff);
 char a[200];
 bool datafile=false;
 bool showortho=false;
 bool showhomo=false;
 bool mtable=false;
 bool ctrules=false;
 bool verbose=false;
 bool pretransf=false;
 bool postransf=false;
 bool mtest=false;
 bool mtalpha=false;
 bool header=true;
 bool isname;
 
 int f=0,lines=0;
 double v;

 
 do{
  ins >> a;
  //cout << a << "<p>" << endl;
  if(!strstr(a,"----")){ 		// ignore tokens with "----timestamp"  
   if(datafile){			// Reading data values
    if(header){
     if(!strstr(a,"DATA")){
      if(!set_factor(a)) return false;
     }
     else{
      header=false;
      set_data_name("DATA");
      f=0;
      lines++;
     } 
    }
    else{
     if(f<get_factors()){
      if(!add_code(f,a,lines)) return false;
      else f++;
     }
     else{
      v=transform((double) atof(a));
      if(!add_value(f,v,lines)) return false;
      f=0;
      lines++;
     } 
    }
   }
   else{ 
    if(showortho){
     set_option("SHOWORTHO",atoi(a));
     showortho=false;
    }
    if(showhomo){
     set_option("SHOWHOMO",atoi(a));
     showhomo=false;
    }
    if(ctrules){
     set_option("CTRULES",atoi(a));
     ctrules=false;
    }
    if(verbose){
     set_option("VERBOSE",atoi(a));
     verbose=false;
    }
    if(mtable){
     set_option("MTABLE",atoi(a));
     mtable=false;
    }  
    if(pretransf){
     set_option("PRETRANSF",atoi(a));
     pretransf=false;
    }
    if(postransf){
     set_option("POSTRANSF",atoi(a));
     postransf=false;
    }   
    if(mtest){
     set_option("MTEST",atoi(a));
     mtest=false;
    }
    if(mtalpha){
     set_alpha(atof(a));
     mtalpha=false;
    }
    if(strstr(a,"octet-stream")) datafile=true; 
    if(strstr(a,"text/plain")) datafile=true;  
    if(strstr(a,"SHOWORTHO")) showortho=true; 
    if(strstr(a,"SHOWHOMO")) showhomo=true;
    if(strstr(a,"CTRULES")) ctrules=true;
    if(strstr(a,"MTABLE")) mtable=true;
    if(strstr(a,"VERBOSE")) verbose=true;
    if(strstr(a,"PRETRANSF")) pretransf=true;
    if(strstr(a,"POSTRANSF")) postransf=true;
    if(strstr(a,"MTEST")) mtest=true;
    if(strstr(a,"ALPHA")) mtalpha=true;
   }
  }
 }while(!ins.eof());
 return true; 
}
#endif

//------------------------------------------------------------------------//
// This function displays a summary of the data (factors, names, levels,  //
// types, and so on)                                                      //
//------------------------------------------------------------------------//

void data::summary()
{
 int fact,lev,i;
 
 #ifdef DEBUG_DATA
 partial *t;
 #endif
 
 if(be_verbose()){
  cout << "---------------- SUMMARY OF FACTORS -----------------" << endl;
  for(fact=0;fact<get_factors();fact++){  
   cout << "Factor " << setw(MAXNAME+1) << get_factor_name(fact) << " (Type: ";
   if(get_factor_type(fact)==RANDOM) cout << "Random) ";
   else cout << "Fixed ) ";
   cout << get_levels(fact) << " levels (";
   for(lev=0;lev<get_levels(fact);lev++){
    if(lev==0) cout << get_code_name(fact,lev);
    else cout << "," << get_code_name(fact,lev); 
   } 
   cout << ") ";
   if(is_nested(fact)){
    cout << " Nested in ";
    for(i=0;i<get_factors();i++) 
    if(is_nested_into(fact,i)) cout << get_factor_name(i) << " ";
   }
   cout << endl;
  }
  cout << "-----------------------------------------------------" << endl;  
 }
 
 #ifdef DEBUG_DATA 
 if(first){
  cout << "Table of partials..." << endl;
  t=first;
  do{
   for(i=0;i<get_factors();i++) cout << (int) t->codes[i];
   cout << " ";
   for(i=0;i<get_factors();i++) cout << (int) t->orig[i];
   cout << " Sum: " << t->sum << "\tSum2: " << t->sum2;
   t->var=(t->sum2-(pow(t->sum,2)/t->n))/(t->n-1);
   cout << "\tVar: " << t->var << "\t n: " << t->n << endl;
   t=t->next;
  }while(t);
 } 
 #endif
}

