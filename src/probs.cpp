// probs.cpp
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
// This file contains code to compute probabilities of several statistics
// It is mainly a collection of code of the excellent lib CCMATH 2.2.0
// developed by Daniel A. Atkinson (DanAtk@aol.com). Applied Statistics
// algorithms were used for the computation of studentized range
// probabilities. Igor Baskir (baskir@univer.kharkov.ua) is the author
// of the code that computes Cochran's C probabilities.

#include <cmath>
#include "conf.h"
using namespace std;

//------------------------------------------------------------------------//
// Algorithm AS66 Applied Statistics (1973) vol22 no.3                    //
// Evaluates the tail area of the standardised normal curve               //
// from x to infinity if upper is .true. or                               //
// from minus infinity to x if upper is .false.                           //
//------------------------------------------------------------------------//
#include "probs.h"
using namespace std;

double alnorm(double x)
{
 const double ltone  = 7.0;
 const double utzero = 18.66;
 const double con    = 1.28;
 const double p      = 0.398942280444;
 const double q      = 0.39990348504;
 const double r      = 0.398942280385;
 const double a1     = 5.75885480458;
 const double a2     = 2.62433121679;
 const double a3     = 5.92885724438;
 const double b1     =-29.8213557807;
 const double b2     = 48.6959930692;
 const double c1     =-3.8052E-8;
 const double c2     = 3.98064794E-4;
 const double c3     =-0.151679116635;
 const double c4     = 4.8385912808;
 const double c5     = 0.742380924027;
 const double c6     = 3.99019417011;
 const double d1     = 1.00000615302;
 const double d2     = 1.98615381364;
 const double d3     = 5.29330324926;
 const double d4     =-15.1508972451;
 const double d5     = 30.789933034;

 double z,y,res;
 bool up;

 up = true;
 z  = x;
 if(!(z >= 0)){
   up = !(up);
   z = -z;
 }
 if((z <= ltone) || up && (z <=utzero)){
   y = 0.5*z*z;
   if(z > con)  res = r*exp(-y)/(z+c1+d1/(z+c2+d2/(z+c3+d3/(z+c4+d4/(z+c5+d5/(z+c6))))));
   else res = 0.5-z*(p-q*y/(y+a1+b1/(y+a2+b2/(y+a3))));
 } 
 else res = 0;
 if(!(up)) res = 1-res;
 return res;
}


//------------------------------------------------------------------------//
// Algorithm AS 190  Appl. Statist. (1983) Vol.32, No. 2                  //
// Incorporates corrections from Appl. Statist. (1985) Vol.34 (1)         //
// Evaluates the probability from 0 to q for a studentized                //  
// range having v degrees of freedom and r samples.                       // 
// Uses subroutine ALNORM = algorithm AS66.                               //
//                                                                        //
// Arrays vw and qw store transient values used in the                    //
// quadrature summation.  Node spacing is controlled by                   //
// step.  pcutj and pcutk control truncation.                             //
// Minimum and maximum number of steps are controlled by                  //
// jmin, jmax, kmin and kmax.  Accuracy can be increased                  //
// by use of a finer grid - Increase sizes of arrays vw                   // 
// and qw, and jmin, jmax, kmin, kmax and 1/step proportionally.          //
// of freedom. Uses AS prtrng and alnorm                                  //
//------------------------------------------------------------------------//

double prtrng(double q, int v, int r) 
{
 short  int ifault;
 double retprtrng;
 double vw[31], qw[31];
 double vmax=120.0;
 double cv1=0.193064705;
 double cv2=0.293525326;
 double cvmax=0.39894228;
 double cv[5]={0, 0.318309886, -0.268132716e-2, 0.347222222e-2, 0.833333333e-1};
 double g, gmid, r1, c, h, v2, gstep, pk1, pk2, gk, pk;
 double w0, pz, x, hj, ehj, pj;
 double pcutj=0.00003, pcutk=0.0001, step=0.45;
 int  jmin=3, jmax=15, kmin=7, kmax=15;
 int  k, j, jj, jump;

  

  /* check initial values */
  retprtrng=0.0;
  ifault=0;
  if ((v < 1.0) || (r < 2.0)) ifault=1;
  if ((q <= 0.0) || (ifault == 1)) return retprtrng;

  /* computing constants, local midpoint, adjusting steps */
  g=step * pow(r, -0.2);
  gmid=0.5 * log((double)r);
  r1=r-1.0;
  c=log(r*g*cvmax);
  if (v <= vmax){ 
  
   h=step * pow(v, -0.5);
   v2=v * 0.5;
   if (v == 1.0) c=cv1;
   if (v == 2.0) c=cv2;
   if (!((v == 1.0) || (v == 2.0)) ) c=sqrt(v2) * cv[1] /(1.0 + ((cv[2]/v2 + cv[3])/ v2 + cv[4]) / v2);
   c=log(c*r*g*h);  
  }
  
  /* computing integral */
at20:
  gstep=g;
  qw[1]=(-1.0);
  qw[jmax+1]=(-1.0);
  pk1=1.0;
  pk2=1.0;

  for (k=1; k<=kmax; k++) {
    gstep=gstep-g;
  at21:
    gstep=(-1.0)*gstep;
    gk=gmid + gstep;
    pk=0.0;
    if ((pk2 <= pcutk) && (k > kmin)) goto at26;
    w0=c-gk*gk*0.5;
    pz=alnorm(gk);
    x=alnorm(gk-q) - pz;
    if (x > 0.0) pk=exp(w0+r1*log(x));
    if (v > vmax) goto at26;

    jump=(-1)*jmax;
  at22:
    jump=jump + jmax;
    for (j=1; j<=jmax; j++) {
      jj=j+jump;
      if (qw[jj] > 0.0) goto at23;
      hj=h*j;
      if (j < jmax) qw[jj+1]=(-1.0);
      ehj=exp(hj);
      qw[jj]=q*ehj;
      vw[jj]=v*(hj+0.5-ehj*ehj*0.5);
    at23:
      pj=0.0;
      x=alnorm(gk-qw[jj]) - pz;
      if (x > 0.0) pj=exp(w0+vw[jj]+r1*log(x));
      pk=pk+pj;
      if (pj > pcutj) goto at24;
      if ((jj > jmin) || (k> kmin)) goto at25;
    at24:
      if (1==2) x=0; // dummy statement;
    }
  at25:
    h=(-1.0)*h;
    if (h < 0.0) goto at22;
at26:
    retprtrng=retprtrng + pk;
    if ((k > kmin) && (pk <= pcutk) && (pk1 <= pcutk)) goto ende;
    pk2=pk1;
    pk1=pk;
    if (gstep > 0.0) goto at21;
  }
  
ende:
  return(retprtrng);;
}

//------------------------------------------------------------------------//
//  Extracted from gaml.c    CCMATH mathematics library source code.      //
//                                                                        //
//  Copyright (C)  2000   Daniel A. Atkinson    All rights reserved.      //
//  This code may be redistributed under the terms of the GNU library     //
//  public license (LGPL). ( See the lgpl.license file for details.)      //
//------------------------------------------------------------------------//

double gaml(double x)
{ 
 double g,h;
 for(g=1.; x<30. ;g*=x,x+=1.); 
 h=x*x;
 g=(x-.5)*log(x)-x+.918938533204672-log(g);
 g+=(1.-(1./6.-(1./3.-1./(4.*h))/(7.*h))/(5.*h))/(12.*x);
 return g;
}

//------------------------------------------------------------------------//
//  Extracted from qgama.c    CCMATH mathematics library source code.     //
//                                                                        //
//  Copyright (C)  2000   Daniel A. Atkinson    All rights reserved.      //
//  This code may be redistributed under the terms of the GNU library     //
//  public license (LGPL). ( See the lgpl.license file for details.)      //
//------------------------------------------------------------------------//


double qgama(double x,double a)
{ 
 double ap,ro,f,t,gaml(double); int k;
 ro=a*log(x)-x-gaml(a);
 ap=6.25; if(a>ap) ap=a; t=(x-ap)/(f=sqrt(2.*ap));
 if(x<4.5 || t< -1. || (t< -.5 && a<10.)){
  for(f=t=1.,ap=a; t>1.e-14 ;){ t*=x/(ap+=1.); f+=t;}
  return 1.-exp(ro)*f/a;
 }
 else{
  if(t<0. && a<10.) k=18;
  else{ 
   if(t>3.){ k=(int) ceil(19.-3.*t); if(k<4) k=4;}
   else k=6+(int)ceil(f*(2.05-.8*t+.091*t*t));
  } 
  for(f=x; k>0 ;){ t=k--; f=x+(t-a)/(1.+t/f);}
  return exp(ro)/f;
 }
}


//------------------------------------------------------------------------//
//  Extracted from qbeta.c    CCMATH mathematics library source code.     //
//                                                                        //
//  Copyright (C)  2000   Daniel A. Atkinson    All rights reserved.      //
//  This code may be redistributed under the terms of the GNU library     //
//  public license (LGPL). ( See the lgpl.license file for details.)      //
//------------------------------------------------------------------------//

double qbeta(double x, double a, double b)
{ 
 double ro,t,ts,f,gaml(double); int nf;
 ro=gaml(a)+gaml(b)-gaml(a+b);
 if(x<.5) nf=1; else{ x=1.-x; t=a; a=b; b=t; nf=0;}
 f=t=exp(a*log(x)+b*log(1.-x)-ro)/a;
 for(ts=0.,b+=a-1.; t>1.e-12 || t>ts ;){
   b+=1.; a+=1.; ts=t; t*=x*b/a; f+=t; 
 }
 if(nf) return f; else return 1.-f;
}


//*********************************************************************//
//                           MAIN FUNCTIONS                            //
//*********************************************************************//

//------------------------------------------------------------------------//
// This function computes F exact probability for v1 and v2  degrees of   //
// freedom. Uses CCMATH's qbeta.                                          //
//------------------------------------------------------------------------//

double fprob(double f, int n1, int n2)
{
 return qbeta(n2/(n2+n1*f),0.5*n2,0.5*n1);
}

//------------------------------------------------------------------------//
// This function computes the inverse of the F distribution for v1 and v2 //
//  degrees of freedom. Uses CCMATH's qbeta.                              //
//------------------------------------------------------------------------//

double fiprob(double f, int n1, int n2)
{
 double z;
 z=qbeta(0.95,(double)0.5*n2,(double) 0.5*n1);
 return (n2*(1-z))/(n1*z);
}

//------------------------------------------------------------------------//
// This function computes Q exact probability for k means and df degrees  // 
// of freedom. Uses AS prtrng and alnorm.                                 //
//------------------------------------------------------------------------//

double qprob(double q, int k, int df)
{
 return 1-prtrng(q,df,k);
}

//------------------------------------------------------------------------//
// This function computes Cochran's C exact probability for k means and   //
// df degrees of freedom.                                                 //
// I must thank Igor Baskir <baskir@univer.kharkov.ua> for this algorithm.//
//------------------------------------------------------------------------//

double cprob(double c, int k, int df)
{
 double P=0;
 if((c>0)&&(k>1)){
  P=fprob((1.0/c-1.0)/(double)(k-1.0),((k-1)*df),df)*k;
  if(P>1) P=fabs(ceil(P)-P);
 }
 return P;
}


//------------------------------------------------------------------------//
// This function computes Chi-square exact probability of X^2 with        // 
// df degrees of freedom.                                                 //
//------------------------------------------------------------------------//


double chiprob(double c, int df)
{ 
 if((c<0.0)||(df<1.0)) return 0;
 return qgama(c*0.5,(double) df*0.5);
}




