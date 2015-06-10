#ifndef NW
#define NW


#include "minhash.h"


//#define DEBUG 0

using namespace std;


extern int  nw(const string&, const string&,string&, string&, bool);

extern int  nw_align(int **, char **,const string&, const string&,string&, string&,int);

extern void  dpm_init        ( int **, char **, int, int, int );
extern void  print_al        ( string&, string& );
extern void  print_matrix    ( int ** const, const string&, const string& );
extern void  print_traceback ( char ** const, const string&, const string& );
extern int   max             ( int, int, int, char * );

#endif
