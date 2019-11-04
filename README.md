# mwanova - Multi-Way Analysis of Variance.

**mwanova** is a simple program to compute analysis of variance. It is distributed under the GNU General Public License. You should read the file COPYING before installing and running mwanova.

Before compiling take a look at the Makefile and *conf.h* and redefine any variables you need. Most of them control the sizes of static variables for data and computed values. Be careful if you need to change them.

This version of mwanova is a revision of the former mwanova-1.1, including an almost complete rewrite of the code which abolished some of its previous limitations (and introduced some bugs). mwanova has still some limitations as any other piece of software. First of all, it  can only deal with as much as 10 factors. I can hardly imagine an experiment with 10 or more factors, but if you need one, just change the definitions MAXFACTORS in the conf.h file (hope you have enough RAM for that...).

The maximum number of data values (MAXDATA), levels per factor (MAXLEVELS) and combinations (COMBINS) have been eliminated. COMBINS is still used but it does not limit the analysis up to 10 factors. Previously, COMBINS was the size of an array of all possible factor combinations for a 10 factor analysis, and was set to 1024. An analysis with 11 factors would have 2048 combinations and was only
possible if this constant was set to 2048. The current version of mwanova  creates a dynamically linked list of terms and the only limitation is RAM. However, COMBINS is still used to compute some string sizes to accommodate variance estimates of terms used in Cornfield-Tukey Rules. If you want to test a model with more than 10 *random* factors you should consider reading *model.cpp* and changing the value of COMBINS in  *conf.h*. Note that this is only true for an analysis with RANDOM factors. An analysis of 20 FIXED factors is easily accomplished by mwanova (well, it takes a while) without any change of COMBINS (but you must change MAXFACTORS).

Another thing that was in the TO DO list is now fully supported in **mwanova**: to automatically handle missing data! However, missing level combinations are not allowed...

**mwanova** does not depend anymore on external libraries for computation  of probabilities. This was possible thanks to Daniel A. Atkinson who developed CCMATH, an excellent library from where portions of code were grabbed. These reside in "probs.cpp". Igor Baskir also contributed with an algorithm to
compute Cochran's C probabilities.

## Installation

Installation of mwanova is straightforward. First unpack the source using

`$ tar zxf mwanova-version.tar.gz`
	
This will make a directory called mwanova-<version> with the source code inside. Just cd into mwanova-<version>

`$ cd mwanova-version`
   
To compile and install mwanova type

`> autoreconf --install`
`> ./configure`
`> make`
`> make install`

Optionally you can configure mwanova to compile as a CGI executable with 

`> ./configure --enable-gci`

and the manual in /usr/local/man/man1 typing

`> make install-man`

(you should do this as root).


This is the second public release of mwanova, and since I have several minor versions (I started it in 1996), I have named it version 2.0. This version is a beta release. I've fixed all the bugs I was able to detect, but probably there are many more. If you find a bug, don't hesitate, and e-mail me. Suggestions are also welcome.

I've only tested mwanova in a Linux system, because I couldn't find a C++ compiler in our HP-UX... if you port the program to other architectures please tell me.

The algorithm is a bit complex. I'm no professional programmer or statistician. I don't have a clue of matrix algebra, but probably most things should have been done using matrixes. Anyway, mwanova works, and it works fine for me. I've compared its results with the results from some commercial packages and everything was OK (I've only tested fully orthogonal  and fully nested models, but as mixed models can be derived from orthogonal models, I did also some tests for mixed models and they were all OK). Now mwanova must pass through all other possible combinations of data and models, and this is only possible if several different people use it.

In the directory data, an example of a data file and model file will produce the following results:

`> mwanova -f Data/data.dat`

```
------------------------ ANOVA RESULTS ---------------------------------
Source of Variation        SS  DF          MS           F       P  Against 
------------------------------------------------------------------------
F1                     46.482   1      46.482       7.345   0.019  Residual
F2                      0.082   1       0.082       0.013   0.911  Residual
F1xF2                  11.207   1      11.207       1.771   0.208  Residual
F3                      4.148   2       2.074       0.328   0.727  Residual
F1xF3                   8.651   2       4.325       0.684   0.523  Residual
F2xF3                   0.766   2       0.383       0.061   0.942  Residual
F1xF2xF3                1.486   2       0.743       0.117   0.890  Residual
Residual               75.940  12       6.328
------------------------------------------------------------------------
TOTAL                 148.760  23   
```         

If you find the program useful, please e-mail me telling so. Don't forget to cite it if you use mwanova in any published paper... thanks, and enjoy it. 

## DONE TO DO's

- [x] Handling of missing data automatically!

- [x] Multiple comparison tests (Tukey, SNK, Schaffe) included in the program!

- [x] Better ANOVA table with terms ordered by importance (main factors, first order interactions, etc).
    
- [x] Incorporated probability functions in the program!
        
## STILL TO DO (if I can find some time to...)

- [ ] Handling of "Non-Testable" F-Ratios automatically

- [ ] Compute asymmetrical designs (not unbalanced data-sets!). Several experiments demand crossed factors and a single external control, not crossed or nested in the other factors. I must understand this first, because I only know how to compute this for a one-way ANOVA (but it is extensible to more factors).

- [ ] Compute BACI (Before-After, Control versus Impact) designs... Very useful when dealing with impact assessment data. Guess this is not easy...

## License
Copyright (C) 2001  Antonio Santos (amsantos@fc.up.pt)

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA



