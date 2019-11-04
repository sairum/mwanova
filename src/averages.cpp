// averages.cpp
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
#include <cmath>
#include "data.h"
#include "probs.h"

using namespace std;

void data::compare(int fact, int total, double err, int dferr, partial *pfirst, partial *plast){
 int 	i,j,k,range;
 double a1,a2,m1,m2,t,p;
 group *q,*qfirst,*qlast;
 bool included;
 partial *t1,*t2;
 
 qfirst=qlast=NULL;
 
 if(pfirst&&plast&&(total>1)){ 
  t1=pfirst;
  j=0;
  do{
   a1=t1->sum/t1->n;
   t2=plast;
   k=total;
   range=k-j;
   do{
    a2=t2->sum/t2->n;
    included=false; 
    // There are homogeneous groups already!
    if(qfirst){ 
     q=qfirst;     
     do{
      m1=q->member1;
      m2=q->member2;
      if((a1>=m1)&&(a1<=m2)&&(a2>=m1)&&(a2<=m2)) included=true;
      q=q->next;
     }while(q&&(!included));
    }
    if(!included){
     // Preform multiple tests....
     if(t1==t2) p=1;
     else{
      t=fabs(a1-a2);
      t=t/sqrt(err/t1->n);
      switch(get_mtest()){
       case SNK: p=qprob(t,range,dferr); break;
       case TUKEY: p=qprob(t,total,dferr);  break;
      }
     } 
     if(p>get_alpha()){
      #ifdef DEBUG_DATA
      cout << a1 << "===" << a2 << "(t=" << t << ",k=" << range << ",err=" << err << ",dferr=" << dferr << ",p=" << p << ")" << endl;
      #endif
      if(!qfirst){
       q= new group;
       q->member1=a1;
       q->member2=a2;
       qfirst=q;
       qlast=qfirst;
       q->next=NULL;
      }
      else{
       q= new group;
       q->member1=a1;
       q->member2=a2;
       qlast->next=q;
       qlast=q;
       q->next=NULL;
      } 
     }
     #ifdef DEBUG_DATA
     else cout << a1 << "!=" << a2 << "(t=" << t << ",k=" << range << ",err=" << err << ",dferr=" << dferr << ",p=" << p << ")" << endl;
     #endif
    }
    t2=t2->prev;
    k--;
   //}while((t2!=t1)&&(t2->prev));  
   }while(t2!=(t1->prev)); 
   t1=t1->next;
   j++;
  }while(t1);

  #ifdef DEBUG_DATA
  if(qfirst){
   q=qfirst;
   do{
    cerr << "(" << q->member1 << "," << q->member2 << ") ";
    q=q->next;
   }while(q);
   cerr << endl;
  }
  #endif
    
  #ifdef CGI
  cout << "<table><theader><tr><td>Level</td><td>Average</td><td>n</td></theader>" << endl;
  #endif
  
  t1=pfirst;
  do{
   #ifndef CGI
   cout << get_code_name(fact,t1->codes[fact]) << "\t" << t1->sum/t1->n << "\t" << t1->n;
   #else
   cout << "<tr><td>" << get_code_name(fact,t1->codes[fact]) << "</td><td>" << t1->sum/t1->n << "</td><td>" << t1->n << "</td>";
   #endif
   if(qfirst){
    q=qfirst;
    do{
     a1=t1->sum/t1->n;
     m1=q->member1;
     m2=q->member2;
     //The 0.000000001 is necessary because float comparison is never precise...
     if((a1>(m1-0.000000001))&&(a1<(m2+0.000000001))){
      #ifndef CGI
      cout << "\t*";
      #else
      cout << "<td>&bull;</td>";
      #endif
     }
     else{
      #ifndef CGI
      cout << "\t ";
      #else
      cout << "<td>&nbsp;</td>";
      #endif
     }     
     q=q->next;
    }while(q);
    #ifndef CGI
    cout << endl;
    #else
    cout << "</tr>";
    #endif
   } 
   else{
    #ifndef CGI
    cout << "\t " << endl;
    #else
    cout << "<td>&nbsp;</td></tr>";
    #endif
   }
   t1=t1->next;
  }while(t1);  
  
  #ifdef CGI
  cout << "</table>" << endl;
  #endif
  
  if(qfirst){
   do{
    q=qfirst->next;
    delete qfirst;
    qfirst=q;
   }while(qfirst); 
  }
 }
}

void data::multi_comp(CODES cline, int fact, const char *fname, double err, int dferr, partial *firstp){
 partial *t,*p,*pfirst,*plast;
 CODES cl;
 combins *firstc,*lastc,*c;
 bool    equal,found;
 int 	 i,npartials;
 
 firstc=NULL;
 lastc=NULL;
 
 pfirst=NULL;
 plast=NULL;
 
 if(firstp){
  t=firstp;
  // Find out all combinations of levels of other factors (excluding 'fact' and factors nested
  // in 'fact') to analyse differences between levels of 'fact'
  do{ 
   memset(cl,0,sizeof(cl));
   for(i=0;i<get_factors();i++){
    if((cline[i]>0)&&(i!=fact)&&(!is_nested_into(i,fact))) cl[i]=t->codes[i]; 
   } 
   //if(found){
    if(!firstc){ // First combination of factors
     c = new combins;      
     memcpy(c->codes,cl,MAXFACTORS);     
     c->next=NULL;
     firstc=c;
     lastc=c;     
    }
    else{
     // Find out if this combination of factors already exists
     c=firstc;
     equal=false;
     do{   
      if(memcmp(cl,c->codes,MAXFACTORS)==0) equal=true;
      c=c->next;
     }while((c)&&(!equal));
      
     if(!equal){
      c=new combins;
      memcpy(c->codes,cl,MAXFACTORS);
      lastc->next=c;
      lastc=c;   
      c->next=NULL;
     }    
    } 
   //}   
   t=t->next;
  }while(t);
  
  
  // Go along the list of combinations of factors, build a list of averages
  // of factor 'fact' for each combination of levels of other factors
  
  if(firstc&&first){    
   c=firstc;
   do{
    t=firstp;
    do{
     memcpy(cl,t->codes,MAXFACTORS);
     cl[fact]=0;
     if(memcmp(cl,c->codes,MAXFACTORS)==0){
      // This is a partial that should be added
      if(!pfirst){
       p = new partial;
       npartials=1;
       memcpy(p->codes,t->codes,MAXFACTORS);
       p->sum=t->sum;
       p->n=t->n;
       pfirst=p;
       plast=p;
       p->next=NULL;
       p->prev=NULL;
      }
      else{
       p = new partial;
       npartials++;
       memcpy(p->codes,t->codes,MAXFACTORS);
       p->sum=t->sum;
       p->n=t->n;
       plast->next=p;
       p->prev=plast;
       plast=p;
       p->next=NULL;
      } 
     }
     t=t->next;
    }while(t);
    
    if(pfirst){
     //#ifdef DEBUG_DATA
     #ifndef CGI
     cout << "For levels of Factor " << get_factor_name(fact);
     for(i=0;i<get_factors();i++){
      if((cline[i]>0)&&(i!=fact)){
       cout << " in level " << get_code_name(i,c->codes[i]) << " of Factor " << get_factor_name(i);
      } 
     } 
     cout << endl;
     #else
     cout << "For levels of Factor <b>" << get_factor_name(fact) << "</b>";
     for(i=0;i<get_factors();i++){
      if((cline[i]>0)&&(i!=fact)){
       cout << " in level " << get_code_name(i,c->codes[i]) << " of Factor <b>" << get_factor_name(i) << "</b>";
      } 
     } 
     cout << endl;
     #endif
     
     #ifdef DEBUG_DATA
     p=pfirst;
     do{
      for(i=0;i<MAXFACTORS;i++) cout << (int) p->codes[i];
      cout << " " << p->sum/p->n << " " << p->n << endl;
      p=p->next;
     }while(p); 
     #endif
     
     compare(fact,npartials,err,dferr,pfirst,plast);
     
     // Destroy partial averages for this combination of levels of factors    
     do{
      p=pfirst;
      pfirst=p->next;
      delete p;
     }while(pfirst); 
    }
    
    c=c->next;
   }while(c);
  }
  
  // Destroy combinations list
  if(firstc){    
   do{
    c=firstc;
    firstc=c->next;
    delete c;
   }while(firstc);    
  }

 }   
}

void data::get_averages(CODES cline, const char *fname, double err, int dferr, double f)
{ 
 partial *t,*firstp,*lastp,*p,q;
 bool    equal,found;
 int     i,j;
 
 firstp=NULL;
 lastp=NULL;
  
 if((!show_mtable())&&(!show_mtests())) return;
 
 if(first){
  if(show_mtable()){
   #ifndef CGI
   cout << "Averages for " << fname << endl;
   #else
   cout << "<CENTER><H3>Averages for " << fname << "</H3></CENTER>";
   #endif
   footer();
  }

  // Create a list of averages
  
  t=first;
  do{   
   if(!firstp){           // No partials in the list, add first
    firstp = new partial;
    if(firstp){
     for(i=0;i<MAXFACTORS;i++){
      if(cline[i]>0) firstp->codes[i]=t->orig[i];
      else firstp->codes[i]=0;
     }
     firstp->sum=t->sum;
     firstp->sum2=t->sum2;
     firstp->n=t->n;
     firstp->next=NULL;
     lastp=firstp;
    }
   }
   else{                  // There are partials in the list...
    p=firstp;    
    do{
     equal=true;
     for(i=0;i<MAXFACTORS;i++){
      if(cline[i]>0) if((p->codes[i])!=(t->orig[i])) {equal=false; break;}
     }
     if(!equal) p=p->next;
    }while((p)&&(!equal));
    
    if(equal){
     
     // There is one item equal to 't->codes'! Add sums and n values to it
     
     p->sum+=t->sum;
     p->sum2+=t->sum2;
     p->n+=t->n;
    }
    else{
     // No item in the list is equal to current 't->codes'. 
     //Add a new one.
     
     p = new partial;
     if(p){
      for(i=0;i<MAXFACTORS;i++){
       if(cline[i]>0) p->codes[i]=t->orig[i];
       else p->codes[i]=0;
      }
      p->sum=t->sum;
      p->sum2=t->sum2;
      p->n=t->n;
      p->next=NULL;
      lastp->next=p;
      lastp=p;
     }
    } 
   }
   t=t->next;
  }while(t);
  
  if(firstp&&show_mtable()){
   // Compute variances
   
   p=firstp;
   do{
    p->var=(p->sum2-pow(p->sum,2)/p->n)/(p->n-1);
    p=p->next;
   }while(p);
   
   #ifndef CGI
   // Output averages
   for(i=0;i<get_factors();i++){
    if(cline[i]>0) cout << setw(4) << get_factor_name(i);
   }
   cout << setw(6) << "n";
   cout << setw(17) << "Average" << setw(17) << "Variance" << endl;
   footer();
   p=firstp;
   do{
    for(i=0;i<get_factors();i++){
     if(cline[i]>0) cout << setw(4) << get_orig_code_name(i,p->codes[i]);
    } 
    cout << setw(6) << p->n;
    cout << setw(17) << p->sum/p->n << setw(17) << p->var << endl;
    p=p->next;
   }while(p); 
   #else
   cout << "<PRE>" << endl;
   for(i=0;i<get_factors();i++){
    if(cline[i]>0) cout << get_factor_name(i) << "\t";
   }
   cout << "n" << "\t";
   cout << "Average" << "\t" << "Variance" << endl;
   
   p=firstp;
   do{
    for(i=0;i<get_factors();i++){
     if(cline[i]>0) cout << get_orig_code_name(i,p->codes[i]) << "\t";
    } 
    cout << p->n << "\t";
    cout <<  p->sum/p->n << "\t" << p->var << endl;
    p=p->next;
   }while(p); 
   cout << "</PRE>";
   #endif
   footer();
  } 
  
  if(firstp&&show_mtests()&&(f<get_alpha())){
   //Are there any fixed factors in cline?
   found=false;
   for(i=0;i<get_factors();i++){
    if((cline[i]>0)&&(get_factor_type(i)==FIXED)) found=true;
   } 
   // No... bail out
   // Yes, proceed with mtests
   if(found){
    #ifndef CGI
    cout << "Multiple comparison tests for " << fname << endl;
    #else
    cout << "<CENTER><H3>Multiple comparison tests for " << fname << "</H3></CENTER>";
    #endif
    footer();
   
    // Sort averages
   
    p=firstp;
    do{
     t=p->next;
     do{
      if((p->sum/p->n)>(t->sum/t->n)){
       memcpy(q.codes,p->codes,MAXFACTORS);
       memcpy(q.orig,p->orig,MAXFACTORS);
       q.sum=p->sum;
       //q.sum2=p->sum2;
       q.n=p->n;
       memcpy(p->codes,t->codes,MAXFACTORS);
       memcpy(p->orig,t->orig,MAXFACTORS);
       p->sum=t->sum;
       //p->sum2=t->sum2;
       p->n=t->n;
       memcpy(t->codes,q.codes,MAXFACTORS);
       memcpy(t->orig,q.orig,MAXFACTORS);
       t->sum=q.sum;
       //t->sum2=q.sum2;
       t->n=q.n;
      }
      t=t->next;
     }while(t);
     p=p->next;
    }while(p->next);
  
    for(i=0;i<get_factors();i++){
     if((cline[i]>0)&&(get_factor_type(i)==FIXED)){
      multi_comp(cline,i,fname,err,dferr,firstp);
     }
    }
   }     
  }
  
  // Destruct list of averages
  if(firstp){ 
   do{
    p=firstp;
    firstp=p->next;
    delete p;
   }while(firstp); 
  }
 }
}



