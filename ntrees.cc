// ntrees.cc   File to demonstrate to test random tree generator
// W.B.Langdon cs.bham.ac.uk
// See ETL-TR-95-35 14 Nov 1995 Hitoshi Iba Random tree generator for GP
// Or Random Generation of Trees, Laurent Alonso and Rene Schott, Kulwer, 1995

#define main_version "19 October 1997 $Revision: 1.5 $"

//WBL 29 Dec 21  Ok with gcc version 4.8.5 and gcc version 9.3.1

//WBL  6 May 01  Tweak so compatible with gcc 2.95.2

//WBL 21 Jan 98  Add display of numberdnz

//WBL 19 Oct 97  Cut down version of rand_tree-test.cc with only tree counting

//WBL 26 Aug 97  Extend to take into account number of each functions

//WBL 22 Aug 97  Use for generating plots of Cn(s) for ant

//WBL 12 Aug 97  Use log of factorial, add remove i[4] loop

//WBL 11 Aug 97  Extend range of size before overflow
//               BUGFIX limit (efficiency saving only), 
//               BUGFIX display of done time

//compiling on Linux g++ -o ntrees -O2 -Wno-write-strings ntrees.cc
//compiling in 2020  g++ -o ntrees -O2 ntrees.cc -lm
//compiling on SUN   g++ -o ntrees -O2 ntrees.cc -lg++ -lm
//compiling on Alpha cxx -o ntrees -O5 ntrees.cc -lm
// Includes main function

#include <iostream>
#include <float.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
using namespace std;

//dont change!
#define max_arity 4

double* log_fact;
int log_fact_size = 0;
void new_log_fact(const int x) {
  log_fact = new double[x+1];
  log_fact_size = x+1;
  double f = 0;
  log_fact[0] = 0;
  for(int i=1;i<log_fact_size;i++) log_fact[i] = (f += log(double(i)));
}
inline
double log_factorial(const int x) {
  assert (x<log_fact_size);
  return log_fact[x];
}

double tree_count(const int size, int arity_count[], const int arity[]) {

  //actual count of distinct tree shapes with a particular choice of
  //node arities is given by Cn(s) (8). However we are interested in
  //the proportion of the total rather than the absolute number. This
  //espacially important since Cn(s) overflows.

  // value returned is log(actual number of trees)
  // arity_count[0] set to number of terminals

  int residue = size;
  {for(int i=1;i<=max_arity;i++) residue -= arity_count[i];}
  arity_count[0] = residue;

  double trees = log_factorial(size)-log(size);
  {for(int i=0;i<=max_arity;i++) trees -= log_factorial(arity_count[i]);}

  //now include number of ways of labeling tree with given terminal+functions
  {for(int i=0;i<=max_arity;i++) 
    if(arity[i] != 0)  trees += arity_count[i] * log(arity[i]);}

  return trees;
}

int limit(const int size, const int arity[], const int index) {
  return (arity[index]==0)? 0 : size / index;
}


double all_arity(const int size, const int arity[]) {
//Given tree size and Arities of functions
//return all possible combinations of arities 

  //max_arity = 4;
  int number = 0;
  double max = -FLT_MAX;
  int best[max_arity+1];
  cout<<size<<flush;
int i[5];
for(i[1]=0;i[1]<=limit((size-1),                   arity,1);i[1]++)
for(i[2]=0;i[2]<=limit((size-1)-i[1],              arity,2);i[2]++)
for(i[3]=0;i[3]<=limit((size-1)-i[1]-i[2]*2,       arity,3);i[3]++) {
  const int r3 = (size-1)-i[1]-i[2]*2-i[3]*3;
  i[4] = r3/4; 
  if(r3%4==0 && (arity[4]>0 || i[4]==0)) {
    number++;
    double tree = tree_count(size,i,arity);
    if(tree>max) {max = tree; memcpy(best,i,sizeof(i)); }
  }
}

//now run it all again but report those some fraction of best

  double total = 0;
  int numbernz = 0;
  int numberdnz= 0;
//cout<<"All possible combinations of arities for tree of size "<<size<<endl;
for(i[1]=0;i[1]<=limit(size,                   arity,1);i[1]++)
for(i[2]=0;i[2]<=limit(size-i[1],              arity,2);i[2]++)
for(i[3]=0;i[3]<=limit(size-i[1]-i[2]*2,       arity,3);i[3]++) {
  const int r3 = (size-1)-i[1]-i[2]*2-i[3]*3;
  i[4] = r3/4; 
  if(r3%4==0 && (arity[4]>0 || i[4]==0)) {
    const double tree = tree_count(size,i,arity);
    const double ratio = exp(tree-max);
    total += ratio;
    if(ratio>FLT_EPSILON) numbernz++;
    if(ratio>DBL_EPSILON) numberdnz++;
}
}
 cout<<"\t"<<number<<" "<<numbernz<<"("<<numberdnz<<")";
 cout<<" "<<total<<flush;
 if(total>0) {
 cout<<" "<<(max+log(total))<<flush;
 cout<<" "<<(max+log(total))/log(10)<<flush; 
 if((max+log(total))<85) {
   cout<<" "<<exp(max+log(total))<<endl;
   return exp(max+log(total));
  }
 else {
   cout<<" -1\n";
   return -1;
   }
 }
 else {
   cout<<" 0 0 0\n";
   return 0;
 }
}//end all_arity

///////////////////////////////  Help
void help() {
  cerr<<" Program to count trees. Version "<<main_version<<endl;
  cerr<<" Command line:  "
      <<" treesize Nleaftypes Narity1 [Narity2] [Narity3] [Narity4]\n";
  cerr<<" \n";
  cerr<<" treesize\t= number of nodes (internal and external) in the tree \n";
  cerr<<"         \t  start,end      \t range of sizes \n";
  cerr<<"         \t  start,end,incr \t ditto but calculate only incrth \n";
  cerr<<" \n";
  cerr<<" Nleaftypes\t= number of different leaves (terminals, arity zero \n";
  cerr<<"           \t  functions) which the trees may contain \n";
  cerr<<" Narity1   \t= number of different internal nodes with \n";
  cerr<<"           \t  one branch (ie functions with one argument) \n";
  cerr<<" Narity2   \t= number of different internal nodes (functions) with\n";
  cerr<<"           \t  two branches (ie two arguments) \n";
  cerr<<" Narity3   \t= three branches \n";
  cerr<<" Narity4   \t= four branches \n";
  cerr<<" \n";
  cerr<<" Output:  7 numbers per line \n";
  cerr<<" \tSize of tree \n";
  cerr<<" \tNumber of arity combinations which yeild trees of this size \n";
  cerr<<" \tNumber of such combinations with >FLT_EPSILON chance of being "
      <<" picked \n";
  cerr<<" \tTotal probability before normalising \n";
  cerr<<" \tNatural Logrithm (Number of trees of this size) \n";
  cerr<<" \tLogrithm Base10(Number of trees of this size) \n";
  cerr<<" \tNumber of trees of this size \n";
  cerr<<" \n";
  cerr<<" \tTotal for all displayed sizes (if available)\n";
  cerr<<" \n";
  exit(EXIT_FAILURE);
}//end help

void error(char* message, char* message2)
{
  cerr<<message<<message2<<endl;
  help();
}

///////////////////////////////  MAIN

int main(int argc, char * argv[])
{
    if(argc<=3||argc>7) error("Need 3, 4, 5 or 6 arguments","");
    int start, end, incr;
    const int t = sscanf(argv[1],"%d,%d,%d",&start,&end,&incr);
    switch(t) {
    case 1: end=start; incr=1; break;
    case 2: incr=1; break;
    case 3: break;
    otherwise: error("Syntax error on ",argv[1]); break;
    }
    if(start<1 || end<1 || incr<1 ) {
      cerr<<start<<","<<end<<","<<incr<<" "<<flush;
      error("Error on ",argv[1]);
    }

    const int arityc[5] = {0,0,0,0,0};
    {for(int i=2; i<argc;i++) {
      if((sscanf(argv[i],"%d",&arityc[i-2])!=1) || (arityc[i-2] < 0) )
	error("Syntax error on ", argv[i]);
    }}
    new_log_fact(end);
    double total = 0;
    {for(int i=start;i<=end;i+=incr) {
      const double t = all_arity(i,arityc);
      if(t<0) total = -1;
      else    total+= t;
    }}
    if(start!=end && total>0) cout<<total<<endl;
    return(EXIT_SUCCESS);
}
