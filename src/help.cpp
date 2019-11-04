// help.cpp
//
// mwanova - Multi-Way Analysis of Variance.
// Copyright (C) 1999  Antonio Santos (amsantos@fc.up.pt)
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
#include <iostream>
using namespace std;

void help()
{
 cout << "mwanova - Multi-Way Analysis of Variance  (Version 2.1.1)" << endl;
 cout << "Copyright (C)  2008  Antonio Santos (amsantos@fc.up.pt)" << endl << endl;
 cout << "This program is distributed under the GNU General Public License" << endl;
 cout << "See the COPYING file for details" << endl << endl;
 cout << "Program usage:" << endl;
 cout << "  mwanova -f DataFile [options]" << endl;
 cout << "  cat DataFile | mwanova -f [options]" << endl << endl;
 cout << "OPTIONS:" << endl;
 cout << "  -a  output means table   -n  do not output anova table" << endl;
 cout << "  -s  no-check mode        -t |sqrt|log|ln|arcsin|  transform data " << endl;
 cout << "  -v  be verbose           -r  output Cornfield-Tukey rules table" << endl;
 cout << endl;
 cout << "Read the man page for more information" << endl;
}
