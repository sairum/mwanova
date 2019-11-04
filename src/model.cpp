// model.cpp
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
#include <cstring>
#include <cstdio>
#include <iomanip>
#include "model.h"
#include "probs.h"

using namespace std;

//************************************************************************//
//**************************** PRVATE STUFF ******************************//
//************************************************************************//


//------------------------------------------------------------------------//
// This function inserts 'terms' in an ordered list of ANOVA terms. The   //
// variable 'fact' holds the factor number being inserted. It is added to //
// 'fcode' (a 1 is inserted in fcode[fact]). 'fcode' is an external       //
// variable because since 'insert_terms' is recursive it must know what   //
// other factors were being analysed when it was called. If 'fact' is     //
// bigger than 0 (first factor), say 2 (factor 3), a term is inserted for //
// factor 2. Then the funciton calls itself for factor 0 and 1. But in    //
// both cases 'fname' will have already a 1 in fname[2]. This means that  //
// when it is called for fact=0, 'fname' will bear 1 in fname[2] and      //
// another 1 will be added in fname[0]. With two 1s this variable is      //
// the interaction between factor 0 and factor 2.                         //
//------------------------------------------------------------------------//

void model::insert_term(int fact)
{
 term *t,*c;
 int  i;
 
 fcode[fact]=1;
 
 if(!first){					
  
  // No terms, insert first one
  
  first = new term;
  memcpy(first->fcode,fcode,MAXFACTORS);
  memset(first->name,0,sizeof(first->name));
  memset(first->ctrow,0,sizeof(first->ctrow));
  strcpy(first->against,"No Test");
  first->SS=0;
  first->MS=0;
  first->ss=0;
  first->df=0;
  first->F=0;
  first->df2=0;
  first->contrast=0;
  first->vars=NULL;
  first->next=NULL;
  first->prev=NULL;
  last=first;
 }  
 else{
  c = new term;
  if(c){
   memcpy(c->fcode,fcode,MAXFACTORS);
   memset(c->name,0,sizeof(c->name));
   memset(c->ctrow,0,sizeof(c->ctrow));
   strcpy(c->against,"No Test");
   c->SS=0;
   c->MS=0;
   c->ss=0;
   c->df=0;
   c->F=0;
   c->df2=0;
   c->contrast=0;
   c->vars=NULL;
   c->next=NULL;
   c->prev=last;
   last->next=c;
   last=c;
  }
 }
 //if(fact>0) for(i=(fact-1);i>=0;i--) insert_term(i); 
 if(fact>0) for(i=0;i<fact;i++) insert_term(i); 
 fcode[fact]=0;
}

//------------------------------------------------------------------------//
// This function returns the order of a term described by 'cline'. The    //
// order is the number of factors involved. Main factors have order 1,    //
// first order interactions have order 2, etc.                            //
//------------------------------------------------------------------------//

int model::get_order(CODES cline)
{
 int i,order=0;
 for(i=0;i<MAXFACTORS;i++) if(cline[i]>0) order++;
 return order;
}

//------------------------------------------------------------------------//
// This function tests if code line 1 has all the items listed in code    //
// line 2. Thus [1001000] is included in [10011010] since the items 0 and //
// 3 in code line 1 are also present in code line 2.                      //
//------------------------------------------------------------------------//

bool model::is_included(CODES c1, CODES c2)
{
 int i;
 for(i=0;i<MAXFACTORS;i++) if((c1[i]>0)&&(c2[i]==0)) return false;
 return true;
}

//------------------------------------------------------------------------//
// This function returns the SS for 'f' term in the ANOVA. SSs for a      //  
// give term are computed in the following way:                           //
//                                                                        //
// a) for main factors (e.g., A)                                          //
//    SS(A)= pSS(A)-CT  (where pSS(A) is the partial SS of factor A)      //
//                                                                        //
// b) for the first order interactions (e.g., AxB)                        //
//    SS(AB)=pSS(AB)-pSS(A)-pSS(B)+CT                                     //
//                                                                        //
// c) for second order interactions (e.g., AxBxC)                         //
//    SS(ABC)=pSS(ABC)-pSS(AB)-pSS(BC)-pSS(AC)+pSS(A)+pSS(B)+pSS(C)-CT    //
//                                                                        //
// and so on... Note the changes in the signs. A given factor only has    //
// a pSS of a factor of a higher order, thus the function 'get_order()'   //
// to rule ou all factors of the same or lower orders.                    // 
//------------------------------------------------------------------------//

double model::get_SS(term *f)
{
 int    ord1,ord2;
 double ss=0;
 term   *t;
 
 ss=f->ss;
 
 if(first){
  t=first;
  do{
   if(get_order(t->fcode)<get_order(f->fcode)){
    if(is_included(t->fcode,f->fcode)){
     ord1=get_order(f->fcode);
     ord2=get_order(t->fcode);
     if(((ord1-ord2)%2)>0) ss-=t->ss;
     else ss+=t->ss; 
    }
   }
   t=t->next;
  }while(t);
 }
 if((get_order(f->fcode)%2)>0) ss-=get_CT();
 else ss+=get_CT();
 return  ss;
}


//------------------------------------------------------------------------//
// Returns the Error MS                                                   //
//------------------------------------------------------------------------//

double model::get_error_ms()
{
 if(get_error_df()>0) return get_error_ss()/get_error_df();
 return 0;
}


//------------------------------------------------------------------------//
// This function builds an orthogonal model from the factors read from the//
// data file. An ANOVA model is made of terms. Each factor has its own    //
// term. Each pair of factors have a first order interaction. Three       //
// factors have 3 first order interactions and 1 second order interaction,//
// and so on...                                                           //
// The algorithm is, again, complex and uses recursive functions to find  //
// all the terms involved in the model. For a fully orthogonal model, the //
// number of terms will be 2^n where n is the number of factors.          //
// Function 'insert_item' is recursive. It will compute the terms for the //
// item being passed ('i') plus all the interactions between factor 'i'   //
// and all other factors considered previously. If i=2 (factor 2) it will //
// insert terms for factor 2, the interaction between factor 1 x factor 2 //
// and the interaction between factor 0 and factor 2.                     //
//------------------------------------------------------------------------//

void model::build_orthogonal_model()
{
 int    i;
 term   *t,*q,s;
 
 // First create a list of all possible combinations (terms) of factors
 
 //memset(fcode,0,sizeof(fcode));
 for(i=0;i<get_factors();i++){
  memset(fcode,0,sizeof(fcode));
  insert_term(i);
 }
 
 // Set their names
  
 if(first){
  t=first;
  do{
   strcpy(t->name,"");
   for(i=0;i<MAXFACTORS;i++){
    if(t->fcode[i]>0){
     if(strlen(t->name)>0) strcat(t->name,"*");
     strcat(t->name,get_factor_name(i));
    }
   }
   t=t->next;
  }while(t);
 }
 
 // Sort terms (factors, first order interactions, ...)

 
 if(first&&(get_factors()>1)){
  t=first;
  do{
   q=t->next;
   do{
    if(get_order(t->fcode)>get_order(q->fcode)){
     memcpy(s.fcode,t->fcode,MAXFACTORS);
     memcpy(s.ctrow,t->ctrow,MAXFACTORS+1);
     s.ss=t->ss;
     s.SS=t->SS;
     s.MS=t->MS;
     s.df=t->df;
     memcpy(s.name,t->name,100);
     memcpy(s.against,t->against,100);
     s.contrast=t->contrast;
     s.F=t->F;
     s.df2=t->df2;
     s.vars=t->vars;
     
     memcpy(t->fcode,q->fcode,MAXFACTORS);
     memcpy(t->ctrow,q->ctrow,MAXFACTORS+1);
     t->ss=q->ss;
     t->SS=q->SS;
     t->MS=q->MS;
     t->df=q->df;
     memcpy(t->name,q->name,100);
     memcpy(t->against,q->against,100);
     t->contrast=q->contrast;
     t->F=q->F;
     t->df2=q->df2;
     t->vars=q->vars;
     
     memcpy(q->fcode,s.fcode,MAXFACTORS);
     memcpy(q->ctrow,s.ctrow,MAXFACTORS+1);
     q->ss=s.ss;
     q->SS=s.SS;
     q->MS=s.MS;
     q->df=s.df;
     memcpy(q->name,s.name,100);
     memcpy(q->against,s.against,100);
     q->contrast=s.contrast;
     q->F=s.F;
     q->df2=s.df2;
     q->vars=s.vars;
    }
    q=q->next;
   }while(q);
   t=t->next;
  }while(t->next);
 }
 
 // Compute the partial sums of squares for each term, and insert 
 // the correspondent degrees of freedom.
 
 #ifdef DEBUG_GET_PARTIAL_SS
 cerr << "DEBUG get_partial_SS(): Table of partial SS " << endl;
 #endif
 
 if(first){
  t=first;
  do{
   t->ss=get_partial_SS(t->fcode); 
   t->df=get_df(t->fcode);
   #ifdef DEBUG_GET_PARTIAL_SS
   cerr << t->name << "\t" << t->ss << "\t" << t->df << endl;
   #endif
   t=t->next;
  }while(t);
 }
 
 #ifdef DEBUG_BUILD_ORTHOGONAL_MODEL
 cerr << "DEBUG build_orthogonal_model(): Partial SS of ANOVA" << endl;
 if(first){
  t=first;
  do{
   for(i=0;i<MAXFACTORS;i++) cerr << (int) t->fcode[i];
   cerr << "\t" << t->name;
   cerr << "\tSS: " << t->ss << "\tdf: " << t->df << endl;
   t=t->next;
  }while(t);
 }
 cerr << endl;
 #endif
 
 // Now compute the real sums of squares fora all terms in the list
 
 if(first){
  t=first;
  do{
   t->SS=get_SS(t); 
   t=t->next;
  }while(t);
 }
 
 #ifdef DEBUG_BUILD_ORTHOGONAL_MODEL
 cerr << "DEBUG build_orthogonal_model(): SS and dfs of ANOVA" << endl;
 if(first){
  t=first;
  do{
   for(i=0;i<MAXFACTORS;i++) cerr << (int) t->fcode[i];
   cerr << "\t" << t->name << "\tSS: " << t->SS << "\tdf: " << t->df << endl;
   t=t->next;
  }while(t);
  cerr << "Error term: " << get_error_ss() << " df: " << get_error_df() << endl;
  cerr << "Total: " << get_total_ss() << " df: " << get_total_df() << endl;
 }
 #endif
}

//------------------------------------------------------------------------//
// This function sets the name of a term in the model. If the term is an  //
// interaction it adds the character '*' between factors. If a factor is  //
// nested in another one it puts the latter between parenthesis.          //
//------------------------------------------------------------------------//

const char *model::set_term_name(CODES cline, char *name)
{

 int i,j;

 strcpy(name,"");
 for(i=0;i<MAXFACTORS;i++){
  if(cline[i]==1){
   if(is_nested(i)){
    if(strlen(name)>0) strcat(name,"*");
    strcat(name,get_factor_name(i));
    strcat(name,"(");
    for(j=0;j<MAXFACTORS;j++){
     if(is_nested_into(i,j)&&(cline[j]==2)){
      strcat(name,get_factor_name(j));
      strcat(name,"*");
     } 
    }
    // The last character of 'name' has an '*' so replace it with ')'
    name[strlen(name)-1]=')';
   }
   else{
    if(strlen(name)>0) strcat(name,"*");
    strcat(name,get_factor_name(i));
   }
  }
 }
 return name;
}


//------------------------------------------------------------------------//
// This function builds the ANOVA model specified by the data file. The   //
// algorithm finds out the SSs and dfs for nested terms and interactions  //
// by combining the necessary terms of the orthogonal model. First, for   //
// each term in the list of terms the program checks if there is or are   //
// nested terms. If so, 'term->fcode' is updated, inserting a '2' in each //
// factor that "nests" a nested factor. For example if 3th factor is      //
// nested in the 2nd factor, its 'fcode' is [0010000] and will be updated //
// to look like [0210000]. If 2nd factor is nested in the 3th, and the 1st//
// factor is nested in the 2nd and 3th, the 'fcode' of the term for the   //
// first factor [1000000] will be updated to [1220000]... complicated, ah?//
// In a second stage, all terms with the same 'fcode' will be clumped     //
// into a single term, summing the SS and dfs, and deleting the rest.     //
//------------------------------------------------------------------------//

void model::build_model()
{
 term  *t,*s,q;
 int   i,j;
 
 if(!show_orthogonal()){   	// If option -o ignore nesting!
  if(first){
   t=first;
   do{
    for(i=0;i<get_factors();i++){
     if(t->fcode[i]>0){
      if(is_nested(i)){
       for(j=0;j<get_factors();j++) 
          if(is_nested_into(i,j)) t->fcode[j]=2; 
      }
     }
    } 
    t=t->next;
   }while(t);
  }
 }
 
 #ifdef DEBUG_BUILD_MODEL
 cerr << "DEBUG build_model(): Recoding 'fcodes'..." << endl;
 if(first){
  t=first;
  do{
   for(i=0;i<MAXFACTORS;i++) cout << (int) t->fcode[i];
   cerr << "\t" << t->name << "\tSS: " << t->SS << "\tdf: " << t->df << endl;
   t=t->next;
  }while(t);
 }
 cerr << endl;
 #endif
 
 // Now add up the SS and df of all the terms with the same codes to the 
 // highest order term, zeroing SS and df of the terms that are added
 
 if(first){
  t=first;
  do{
   s=t->next;
   if(s){       
    do{
     #ifdef DEBUG_BUILD_MODEL
     cerr << "Comparing ";
     for(i=0;i<MAXFACTORS;i++) cerr << (int) t->fcode[i];
     cerr << " with ";
     for(i=0;i<MAXFACTORS;i++) cerr << (int) s->fcode[i];
     cerr << endl;
     #endif
     if(memcmp(t->fcode,s->fcode,MAXFACTORS)==0){
      #ifdef DEBUG_BUILD_MODEL
      cerr << "Adding them " << endl;
      #endif
      t->SS+=s->SS;
      t->df+=s->df;
      strcpy(t->name,set_term_name(t->fcode,t->name));
      s->SS=0;
      s->df=0;
      memset(s->fcode,0,sizeof(fcode));
      strcpy(s->name,"");
     }
     s=s->next;
    }while(s);
   }
   t=t->next;
  }while(t);
 }
 
 // Now delete entries with SS and df equal to zero
 
 if(first){
  t=first;
  do{
   if(get_order(t->fcode)==0){
    if(t==first){
     s=t;
     t=t->next;
     t->prev=NULL;
     first=t;
     delete t; 
    }
    else{
     if(t==last){
      s=t;
      last=t->prev;
      last->next=NULL;
      t=NULL;
      delete s;
     }
     else{
      s=t;
      t->prev->next=t->next;
      t->next->prev=t->prev;
      t=t->next;
      delete s;
     }
    }     
   }
   else{
    t->MS=t->SS/t->df;
    t=t->next;
   } 
  }while(t);
 }
 
 #ifdef DEBUG_MODEL
 if(first){
  t=first;
  do{
   for(i=0;i<MAXFACTORS;i++) cout << (int) t->fcode[i];
   cout << "\t" << t->name << "\tSS: " << t->SS << "\tdf: " << t->df << endl;
   t=t->next;
  }while(t);
 }
 cerr << "Error term: " << get_error_ss() << " df: " << get_error_df() << endl;
 cerr << "Total: " << get_total_ss() << " df: " << get_total_df() << endl;
 #endif
 
 // Now reorder terms, firts those with order=0, then order=1, and so on

 if(first){
  t=first;
  while(t->next){
   s=t->next; 
   if(s){
    do{
     if(get_order(t->fcode)>get_order(s->fcode)){
      memcpy(q.fcode,s->fcode,MAXFACTORS);
      memcpy(q.ctrow,s->ctrow,MAXFACTORS+1);
      q.ss=s->ss;
      q.SS=s->SS;
      q.MS=s->MS;
      q.df=s->df;
      memcpy(q.name,s->name,100);
      memcpy(q.against,s->against,100);
      q.contrast=s->contrast;
      q.F=s->F;
      q.df2=s->df2;
      q.vars=s->vars;  
      
      memcpy(s->fcode,t->fcode,MAXFACTORS);
      memcpy(s->ctrow,t->ctrow,MAXFACTORS+1);
      s->ss=t->ss;
      s->SS=t->SS;
      s->MS=t->MS;
      s->df=t->df;
      memcpy(s->name,t->name,100);
      memcpy(s->against,t->against,100);
      s->contrast=t->contrast;
      s->F=t->F;
      s->df2=t->df2;
      s->vars=t->vars;
      
      memcpy(t->fcode,q.fcode,MAXFACTORS);
      memcpy(t->ctrow,q.ctrow,MAXFACTORS+1);
      t->ss=q.ss;
      t->SS=q.SS;
      t->MS=q.MS;
      t->df=q.df;
      memcpy(t->name,q.name,100);
      memcpy(t->against,q.against,100);
      t->contrast=q.contrast;
      t->F=q.F;
      t->df2=q.df2;
      t->vars=q.vars;
     }
     s=s->next; 
    }while(s);
   } 
   t=t->next; 
  }
 }
 
 #ifdef DEBUG_MODEL
 if(first){
  t=first;
  do{
   for(i=0;i<MAXFACTORS;i++) cout << (int) t->fcode[i];
   cout << "\t" << t->name << "\tSS: " << t->SS << "\tdf: " << t->df << endl;
   t=t->next;
  }while(t);
 }
 cerr << "Error term: " << get_error_ss() << " df: " << get_error_df() << endl;
 cerr << "Total: " << get_total_ss() << " df: " << get_total_df() << endl;
 #endif
}


//------------------------------------------------------------------------//
// Computes the entry in the multiplier table for subscript j and term i  //
// It's used in ctrules()                                                 //
//------------------------------------------------------------------------//

int model::table_entry(int fact, CODES f)
{ 
 if(f[fact]>0){		// Factor is present in the term  
  if(f[fact]==2) // Factor is in parenthesis (the term factor is nested in it)
         return 1;       
  else{
   if(get_factor_type(fact)==RANDOM)  // Factor is random
         return 1;
   else return 0;                 // Factor is fixed
  }
 }
 else return (int) get_levels(fact);
 return 0;
}

//------------------------------------------------------------------------//
// This function checks if variance term with code c1 is component of     //
// variance term with code c2. It is used in ctrules()                    //
//------------------------------------------------------------------------//

bool model::is_component(CODES c1, CODES c2)
{
 int i;
 
 for(i=0;i<MAXFACTORS;i++)
  if((c2[i]>0)&&(c1[i]==0)) return false;
 return true; 
}

//------------------------------------------------------------------------//
// Computes Cornfield-Tukey Rules to find ou what are the F ratios to be  //
// computed.                                                              //
//------------------------------------------------------------------------//

void model::ctrules()
{
 int  i,l,coef;
 term *t,*s;				   
 char tempa[MAXTERMSIZE],tempb[MAXTERMSIZE],tempc[100];
 char *p;
 
 // Build Cornfield-Tukey table of multipliers. This table is necessary 
 // to find out the variance components that are present in each term   
 // so as to find the appropriate F ratios to be computed. The table,
 // itself, is separated by terms, that is, each row of the table is
 // in one of the terms of the analysis (each row is called 'ctrow')               

 if(first){
  t=first;
  do{
   memset(t->ctrow,0,sizeof(t->ctrow));
   for(i=0;i<MAXFACTORS;i++) t->ctrow[i]=table_entry(i,t->fcode);

   // 'ctrow' has one more cell for the Error. For all terms this
   // cell is filled with the number of replicates
   
   t->ctrow[MAXFACTORS+1]=get_n();

   t=t->next;
  }while(t);
  
  #ifdef DEBUG_CTRULES
  cerr << "DEBUG ctrules(): Table of multipliers" << endl;
  cerr << "-------- CORNFIELD-TUKEY TABLE OF MULTIPLIERS --------" << endl;
  t=first;
  do{
   cerr << setw(20) << setiosflags(ios::left);
   cerr << t->name;
   for(i=0;i<get_factors();i++) 
     cerr << setw(5) << setiosflags(ios::left) << t->ctrow[i];   
   cerr << setw(5) << setiosflags(ios::left) << t->ctrow[MAXFACTORS+1] << endl;
   t=t->next;
  }while(t);
  cerr << "------------------------------------------------------" << endl;
  #endif

 
  // For every term in the list check what are the variance components 
  // measured (or included) in it. To do so contrast the term with all
  // terms in the analysis starting from below. If a variance component 
  // is present in the term being considered add it to 'tempa' multiplied 
  // by a coeficient computed from the multiplier table. Store the variance 
  // attributable to the term itself in a separate variable to append in 
  // the end (this must be done since as factors are inserted in the list in
  // an unordered way, sometimes the term being analysed involves 
  // components which are above in the list, which means that they would be 
  // appended after the variance attributable to the term itself; the F tests 
  // are computed by comparing the 'vars' string for one term *minus the last
  // variance component (attributable to himself)* with all other terms).
 
  t=first;
  do{       
   s=last;
   #ifndef CGI
   //sprintf(tempa,"%de",get_n());
   sprintf(tempa,"e");
   #else
   strcpy(tempa,"&sigma;&sup2;<SUB><I>e</I></SUB>");
   #endif
   do{
    if(is_component(s->fcode,t->fcode)){
     coef=get_n();
     for(i=0;i<get_factors();i++){      
      if(t->fcode[i]==0){
       coef*=s->ctrow[i];
      }
     }
     if(coef>0){
      #ifndef CGI
      sprintf(tempb,"+%d%s",coef,s->name);
      #else
      sprintf(tempb,"+%d&sigma;&sup2;<SUB><I>%s</I></SUB>",coef,s->name);
      #endif
      if(t==s) strcpy(tempc,tempb); // Store self component
      else strcat(tempa,tempb);     // ... or add to term list
     } 
    }        
    s=s->prev;
   }while(s);
   
   // Add the self component variance
   
   strcat(tempa,tempc);
   
   l=strlen(tempa);
   t->vars = new char[l+1];  	
   strcpy(t->vars,tempa);
   
   t=t->next;
  }while(t);
  
  // Now insert Error Term in the list
  
  s = new term;
  
  strcpy(s->name,"Error");
  memset(s->fcode,0,sizeof(s->fcode));
  memset(s->ctrow,0,sizeof(s->ctrow));
  s->SS=get_error_ss();
  s->df=get_error_df();
  s->MS=get_error_ms();
  #ifndef CGI
  //sprintf(tempa,"%de",get_n()); 
  sprintf(tempa,"e");
  #else
  strcpy(tempa,"&sigma;&sup2;<SUB><I>e</I></SUB>");  
  #endif 
  l=strlen(tempa);
  s->vars = new char[l+1];  	
  strcpy(s->vars,tempa);
   
  s->next=NULL;
  s->prev=last;
  last->next=s;
  last=s;  
  
  // Check what are the tests to be done and compute F ratios

  t=first;
  do{
   strcpy(tempa,t->vars);
   p=rindex(tempa,'+');      
   p[0]=(char) 0;
   s=last;
   do{  
    if(strcmp(tempa,s->vars)==0){
     strcpy(t->against,s->name);
     t->contrast=s->MS;
     if((s->MS>0)&&(s->df>0)){
      t->F=(t->MS/s->MS);
      t->df2=s->df;
     } 
    }      
    s=s->prev;
   }while(s);
   
   t=t->next;
   
   // It is 'while(t->next)' because one doesn't want to contrast the last 
   // term (Error or Residual) with any other...
   
  }while(t->next);
  if(show_ctrules()){
   #ifndef CGI
   header(" Cornfield-Tukey Rules ");
   t=first;
   do{
    cout << setiosflags(ios::left);
    cout << setw(20) << t->name << " " << t->vars << endl;    
    t=t->next;
   }while(t);
   #else
   header(" Cornfield-Tukey Rules ");
   cout << "<TABLE>";
   t=first;
   do{
    cout << "<TR><TD>" << t->name << "</TD><TD>" << t->vars << "</TD></TR>" << endl;    
    t=t->next;
   }while(t);
   cout << "</TABLE>" << endl;
   #endif
  }
  footer();
 }
}

//------------------------------------------------------------------------//
// Computes averages for main factors and interactions and preforms       //
// multiple comparison tests as needed.                                   //
//------------------------------------------------------------------------//

void model::averages()
{
 term *t;
 if((show_mtable()||show_mtests())&&first){
  t=first;
  do{
   get_averages(t->fcode, t->name, t->contrast, t->df2, fprob(t->F,t->df,t->df2));
   t=t->next;
   
   // It is 'while(t->next)' because one doesn't want to compute 
   // averages for the Error or Residual term...
  }while(t->next);
 }
}

//************************************************************************//
//**************************** PUBLIC FUNCTIONS **************************//
//************************************************************************//

model::model()
{
 first=NULL;
 last=NULL;
 memset(fcode,0,sizeof(fcode));
 memset(fname,0,sizeof(fname));
}

model::~model()
{
 term *t;
 if(first){
  do{
   t=first->next;
   if(first->vars) delete [] first->vars;
   delete first;
   first=t;
  }while(first); 
 }
 #ifdef DEBUG_MODEL
 cout << "Destructing 'model' variable" << endl;
 #endif
}

void model::write_anova()
{
 unsigned int   NS;
 term *t;
 
 // Check for the size of the longest term name
 
 NS=0;
 if(first){
  t=first;
  do{
   if(strlen(t->name)>NS) NS=strlen(t->name);
   t=t->next;
  }while(t);
 }
 if(NS<20) NS=20;
 
 // Print the ANOVA table 
 #ifndef CGI
 header(" ANOVA Results ");
 cout << setiosflags(ios::left);
 cout << setw(NS) << "Source of Variation";
 cout << resetiosflags(ios::left) << setiosflags(ios::right);
 cout << setw(SSSIZE) << "SS";
 cout << setw(DFSIZE) << "DF";
 cout << setw(MSSIZE) << "MS";
 cout << setw(FRSIZE) << "FR";
 cout << setw(PRSIZE) << "P";
 cout << " Against" << endl;
 footer();
 
 if(first){
  t=first;
  do{
   cout << setiosflags(ios::left);
   cout << setw(NS) << t->name;   
   cout << resetiosflags(ios::left);
   cout << setiosflags(ios::right|ios::fixed);
   cout << setprecision(PRECISION);
   cout << setw(SSSIZE) << t->SS;
   cout << setw(DFSIZE) << t->df;
   cout << setw(MSSIZE) << t->MS;
   if(t!=last){
    if(t->df2>0){ 
     cout << setw(FRSIZE) << t->F; 
     cout << setw(PRSIZE) << fprob(t->F,t->df,t->df2);
     cout << " " << t->against << endl;
    }
    else{
     cout << setw(FRSIZE) << "-";
     cout << setw(PRSIZE) << "-";
     cout << " " << t->against << endl;
    }
   }   
   else cout << endl;
   t=t->next;
  }while(t);
 } 
 footer();
 cout << setiosflags(ios::left);
 cout << setw(NS) << "Total";
 cout << resetiosflags(ios::left);
 cout << setiosflags(ios::right|ios::fixed);
 cout << setprecision(PRECISION);
 cout << setw(SSSIZE) << get_total_ss();
 cout << setw(DFSIZE) << get_total_df() << endl << endl;
 #else
 header("ANOVA Results");
 cout << "<TABLE RULES='groups' CELLPADDING=8>" << endl;
 cout << "<THEAD>" << endl;
 cout << "<TR>" << endl;
 cout << "<TD>Source of Variation</TD>";
 cout << "<TD>SS</TD>";
 cout << "<TD>DF</TD>";
 cout << "<TD>MS</TD>";
 cout << "<TD>FR</TD>";
 cout << "<TD>P</TD>";
 cout << "<TD>Against</TD></TR>" << endl;
 cout << "</THEAD>" << endl;
 cout << "<TBODY>" << endl; 
 if(first){
  t=first;
  do{
   cout << "<TR><TD>" << t->name << "</TD>";   
   cout << resetiosflags(ios::left);
   cout << setiosflags(ios::right|ios::fixed);
   cout << setprecision(PRECISION);
   cout << "<TD>" << t->SS << "</TD>";
   cout << "<TD>" << t->df << "</TD>";
   cout << "<TD>" << t->MS << "</TD>";
   if(t!=last){
    if(t->df2>0){ 
     cout << "<TD>" << t->F << "</TD>"; 
     cout << "<TD>" << fprob(t->F,t->df,t->df2) << "</TD>";
     cout << "<TD>" << t->against << "</TD>" << endl;
    }
    else{
     cout << "<TD>" << "-" << "</TD>";
     cout << "<TD>" << "-" << "</TD>";
     cout << "<TD>" << t->against << "</TD>" << endl;
    }
   }   
   else cout << "</TR>" <<  endl;
   t=t->next;
  }while(t);
 } 
 cout << "</TBODY>" << endl;
 cout << "<TFOOT>" << endl;
 cout << setiosflags(ios::left);
 cout << setw(NS) << "<TR><TD>Total</TD>";
 cout << resetiosflags(ios::left);
 cout << setiosflags(ios::right|ios::fixed);
 cout << setprecision(PRECISION);
 cout << "<TD>" << get_total_ss() << "</TD>";
 cout << "<TD>" << get_total_df() << "</TD>";
 cout << "<TD></TD><TD></TD><TD></TD><TD></TD></TR>" << endl;
 cout << "</TFOOT>" << endl;
 cout << "</TABLE>" << endl;
 footer();
 #endif
}

void model:: run()
{ 
 equalize();
 compute_nesting();
 if(orthogonalize()){
  summary();
  test_homogeneity();
  build_orthogonal_model();
  build_model(); 
  ctrules();
  write_anova();
  averages();
 } 
}
