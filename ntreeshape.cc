// ntrees.cc   File to demonstrate to calculate number of trees of each depth
// W.B.Langdon cwi.nl
// based on ntrees.cc r1.2

#define main_version "12 March 1999 $Revision: 1.22 $"

//WBL 20 Aug 02  Add code from ntrees.cc r1.3 to make self contained 
//WBL  5 Apr 99  Add command line choice and make cache a proper class
//               Add calculation of mean shape
//WBL 16 Mar 99  re-definition of depth so root depth=1 (not inside count)
//WBL 14 Mar 99  Add gradient
//WBL 12 Mar 99  Add calculation of percentiles
//WBL 14 Dec 98  

//compile g++ -o ntreeshape -O3 ntreeshape.cc
//Running on   SGI   ./ntreeshape 1,141,2 1,100

// Includes main function

#include <iostream.h>
#include <float.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

#define BOOL int
#define TRUE 1
#define FALSE 0

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

//Binary trees only

class Cache {
private:
  enum {INVALID = -1 };
  enum {invalid = '\xff' };
  union datatype {double value; int flag;};
  datatype* data;
  int size, depth, stride;
  inline int addressS(const int s) const { return s/2; }; //only binary trees
  inline int address(const int s, const int d, int& address) const {
    if((s%2)!=0 && s>=0 && s<=size && d>=0 && d<=depth) { //only binary trees
      address = addressS(s)*stride+d;
      return TRUE;
    }
    return FALSE;
  };
public:
  Cache(const int s, const int d): size(s), depth(d), stride(d+1) {
    const int datasize = (addressS(size)+1)*stride;
    data = new datatype[datasize];
    memset(data,invalid,sizeof(datatype)*datasize);
  };
  ~Cache() { delete[] data; };
  inline BOOL cached(const int s, const int d, double& out) const {
    int a;
    if(address(s,d,a) && data[a].flag != INVALID ) {
      out = data[a].value;
      return TRUE;
    }
    else return FALSE;
  }
  inline void set(const int s, const int d, const double in) {
    int a;
    if(address(s,d,a)) data[a].value = in;
  }
};

Cache* cache;

//number of binary trees of size and max depth (root node is depth=0,size=1)
double count(const int size, const int depth /*, const char* ds*/) {
  double ans;
  if(size<=0||depth<0)//||size<=2*depth)
    ans = 0;
  else if(depth==0)
    ans = (size==1)? 1 : 0;
  else {
    double d;
    if(cache->cached(size,depth,d)) return d;

    double sum1 = 0; double sum2 = 0;
    {double cmd1=1;//dummy initial value
    for(int m=2*depth-1;(m<size-1)&&(cmd1>0);m+=2) {
      cmd1 = count(m,depth-1 /*,"1a"*/);
      for(int s=0;s<depth-1;s++)
	sum1 += cmd1*count(size-m-1,s /*,"1b"*/);
    }}
    {double c1=1; double c2=1;//dummy initial values
    for(int m=2*depth-1;(m<size)&&((size-m-1)>=(2*depth-1))&&c1>0;m+=2){
      c1 = count(m,       depth-1 /*,"2a"*/);
      c2 = count(size-m-1,depth-1 /*,"2b"*/);
      sum2 += c1*c2;
    }}
    ans = 2*sum1+sum2;
  }
  //cout<<"count("<<size<<","<<depth<<")="<<ans<<endl;
  cache->set(size,depth,ans);
  return ans;
}//end count

//as count but do quick sanity check first (root node is depth=1,size=1)
double count2(const int i, const int d /*, const char* ds*/) {
  return ((i%2==1)&&(i<pow(2,d)))? count(i,d-1 /*,"x"*/) : 0;
}


///////////////////////////////  Help
void help() {
  cerr<<" Program to count Binary tree shapes. Version "<<main_version<<endl;
  cerr<<" Command line:  "
      <<" treesize treedepth format\n";
  cerr<<" \n";
  cerr<<" treesize\t= number of nodes (internal and external) in the tree \n";
  cerr<<"         \t  start,end      \t range of sizes \n";
  cerr<<"         \t  start,end,incr \t ditto but calculate only incrth \n";
  cerr<<" \n";
  cerr<<" depth   \t= maxdistance from root to leaf in the tree \n";
  cerr<<"         \t  start,end      \t range of sizes \n";
  cerr<<"         \t  start,end,incr \t ditto but calculate only incrth \n";
  cerr<<"         \t  (depth is defined so the root is size=1 and depth=1) \n";
  cerr<<" format: \n";
  cerr<<" \tS\tsize first \n";
  cerr<<" \tD\tdepth first \n";
  cerr<<" \tP\tpercentages \n";
  cerr<<" \tG\tgradient arrows (for gnuplot) \n";
  cerr<<"  omitted\tall four options \n";
  cerr<<" \n";
  cerr<<" Output: (S and D options) 4 numbers per line \n";
  cerr<<" \tSize of tree \n";
  cerr<<" \tDepth of tree \n";
  cerr<<" \tNumber of binary trees with this combination of size and shape \n";
  cerr<<" \tNumber of trees of this size \n";
  cerr<<" \tLog10 this number (0 represented by 0.01)\n";
  cerr<<" \n";
  cerr<<" \tTotal for all displayed sizes (if available)\n";
  cerr<<" \n";
  exit(EXIT_FAILURE);
}//end help

void error(const char* message, const char* message2)
{
  cerr<<message<<message2<<endl;
  help();
}

int parsearg(const char* s, int* start, int* end, int* incr, const int min) {
  //cout<<"parsearg `"<<s<<"'\n";
  int dummy;
  const int t = sscanf(s,"%d,%d,%d,%d",start,end,incr,&dummy);
  switch(t) {
  case 1: *end=*start; *incr=1; break;
  case 2: *incr=1; break;
  case 3: break;
  default: error("Syntax error on ",s); break;
  }
  if(*start<min || *end<min || *incr<1 ) {
    cerr<<*start<<","<<*end<<","<<*incr<<" "<<flush;
    error("Error on ",s);
  } 
  return t;
}
void pr1(const double t, const char* sep)
{
  cout<<t<<sep<<flush;
  if(t>0) cout<<log10(t);
  else    cout<<"0.01";
}
void pr2(const int d, const double t,const char* sep)
{
  cout<<d<<sep<<flush;
  pr1(t,sep);
}
void pr(const int i, const int d, const double t, const char* s)
{
  cout<<i<<"\t"<<flush;
  pr2(d,t,"\t");
  if(s!=0)cout<<"\t"<<s;
  cout<<endl;
}
///////////////////////////////  MAIN

class pair {
  void set(const int d, const double v) { loc = d; value = v;}
public:
  int loc;
  double value, threshold;
  pair(): loc(0), value(0), threshold(0) {;};
  void thresh(const double v) {threshold=v;}
  void update(const int d, const double t, const double v) {
    if(value==0 && t>threshold) set(d,v);
  }
  void peakset(const int d, const double v) {
    if(v>value) set(d,v);
  }
  void cset(const int d, const double v) {
    if(loc==d) value = v;
  }
  void print() { pr2(loc,value," "); cout<<"\t"; }
};//endclass pair
class mm {
  pair min, pl, peak, pu, max;
  double sum;
  double estimatedmean; //for improved numerical stability
  double msum, msum2;
public:
  mm(const int size) {
          int binary[5] = {0,0,size/2,0,0};
    const int arityc[5] = {1,0,1,0,0};
    const double total = exp(tree_count(size,binary,arityc));
    pl.thresh(total*0.05);
    pu.thresh(total*0.95);
    min.loc  = 1+int(log(size)/log(2)); //(root node is depth=1,size=1)
    pu.loc   = min.loc;
    peak.loc = min.loc;
    pl.loc   = min.loc;
    max.loc  = (size+1)/2;
    max.value= exp((max.loc-1)*log(2));
    sum = 0;
    const double pi = 3.1415927;
    estimatedmean = 2*sqrt(pi*(size-1)/2);//sedgewick p256 T5.8 ignore O(n**0.25)
    msum = 0;
    msum2 = 0;
  }
  void update(const int d, const double t) {
    const double dev = d-estimatedmean;
    msum   += dev*t;
    msum2  += dev*dev*t;
    sum += t;
    min.cset(d,t);
    pl.update(d,sum,t);
    peak.peakset(d,t);
    max.cset(d,t);
    pu.update(d,sum,t);
  }
  void print(const int i, const char* s) {
    cout<<i<<"\t"<<flush;
    pr1(sum," ");cout<<" "<<flush;
    min.print();
    pl.print();
    peak.print();
    pu.print();
    max.print();
    if(s!=0)cout<<"\t"<<s;
    const double m    = msum/sum;
    const double mean = m+estimatedmean;
    const double var  = msum2/sum-m*m;
    const double sd   = (var>0)? sqrt(var) : 0;
    cout<<"\t"<<mean;
    cout<<"\t"<<sd;
    cout<<"\t"<<estimatedmean;
    cout<<endl;
  }
};//endclass mm;

int main(int argc, char * argv[])
{
    //cout<<argc<<endl;
    if(argc<3) error("Need 2 or 3 arguments","");
    int  start,  end,  incr;
    int dstart, dend, dincr;
    const int t  = parsearg(argv[1], &start, &end, &incr,1);
    const int t2 = parsearg(argv[2],&dstart,&dend,&dincr,0);

    cout<<"#Program to count binary tree shapes Version "<<main_version<<endl;
    cout<<"# "<<argv[1]<<" t "<<t<<" t2 "<<t2<<endl;
    cout<<"#size "<< start<<";"<< end<<";"<< incr<<"\t"<<flush;
    cout<<"depth "<<dstart<<";"<<dend<<";"<<dincr<<"\t"<<endl;
    cache = new Cache(end,dend);

    BOOL sizefirst  = TRUE; // S
    BOOL depthfirst = TRUE; // D
    BOOL percent    = TRUE; // P
    BOOL gradient   = TRUE; // G
    if(argc==4) {
      sizefirst  = FALSE; // S
      depthfirst = FALSE; // D
      percent    = FALSE; // P
      gradient   = FALSE; // G
      switch (toupper(argv[3][0])) {
      case 'S': sizefirst = TRUE; break;
      case 'D': depthfirst= TRUE; break;
      case 'P': percent   = TRUE; break;
      case 'G': gradient  = TRUE; break;
      default: error("Unknown option",argv[3]);
      }
    }

    double total = 0;
    if(sizefirst) {
    cout<<"#size first\n";
    for(int i=start;i<=end;i+=incr) {
      for(int d=dstart;d<=dend;d+=dincr) {
	const double t = count2(i,d);
	pr(i,d,t,0);
      }
      cout<<endl;
    }
    cout<<endl;
    }

    if(depthfirst) {
    cout<<"#depth first\n";
    for(int d=dstart;d<=dend;d+=dincr) {
      double lastt=-1;
      for(int i=start;i<=end;i+=incr) {
	const double t = count2(i,d);
	if(t>0&&lastt==0)     pr(i-incr,d,lastt,"df");
	if(t>0||lastt>0)      pr(i  ,   d,    t,"df");
	lastt=t;
      }
      cout<<endl;
    }
    cout<<endl;
    }
    if(total>0) cout<<total<<endl;
    
    if(percent||gradient) new_log_fact(end);

    if(percent) {
    cout<<"#size percentiles\n";
    for(int i=start;i<=end;i+=incr) {
      mm data(i);
      for(int d=dstart;d<=dend;d+=dincr) {
	data.update(d,count2(i,d));
      }
      data.print(i,"ps");
    }}

/**********************************************************************
    NOT USED -- was too difficult to read scaled arrows in gnuplot
    cout<<"#gradient, scaling prepass "<<flush;
    {double max_z = 0;
     for(int i=start;i<=end;i+=incr) {
      for(int d=dstart;d<=dend;d+=dincr) {
	const double f0  = count2(i,  d);
	const double dfd = count2(i,  d+1) - f0;
	const double dfs = count2(i+2,d)   - f0;
	if(dfd>max_z) max_z = dfd; //acurate enough for log,log scaling
	if(dfs>max_z) max_z = dfs;
      }
    }
    double llmaxz = log10(log10(max_z));
    cout<<max_z<<" "<<llmaxz<<endl;
**********************************************************************/
    if(gradient) {
    cout<<"#gradient\n";
    for(int i=start;i<=end;i+=incr) {
      mm data(i);
      for(int d=dstart;d<=dend;d+=dincr) {
	const double f0 = count2(i,d);
	if(f0>0) {
	  const double dfd = count2(i,  d+1) - f0;
	  const double dfs = count2(i+2,d)   - f0;
	  const double theta = atan2(dfs,dfd);
	  const double sint = sin(theta);
	  const double cost = cos(theta);
	  //const double z1 = fabs(dfd*sint+dfs*cost);
	  //const double z = (z1<max_z/60.017709)? .25 : log10(log10(z1))/llmaxz;
	  const double z = 0.67;
	  cout<<"set arrow from ";
	  cout<<d<<","<<i<<" to "<<d+z*cost<<","<<i+2*z*sint<<flush;
	  //cout<<"#\t"<<f0<<"\t"<<dfs<<"\t"<<dfd<<"\t"<<theta<<" "<<(theta*180.0/3.1415927)<<"\t"<<z1<<" "<<z;
	  cout<<endl;
	  }
    }}}
    return(EXIT_SUCCESS);
}
