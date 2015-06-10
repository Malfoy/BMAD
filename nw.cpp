//
//  nw.cpp
//  PBMOG
//
//  Created by malfoy on 11/03/2015.
//  Copyright (c) 2015 malfoy. All rights reserved.
//

#include "nw.h"

using namespace std;


int nw(const string& seq_1, const string& seq_2, string& seq_1_al, string& seq_2_al, bool prm){
	int  d = 2 ;  /* gap penalty */

	int  L1 = (int)seq_1.length();
	int  L2 = (int)seq_2.length();

	// Dynamic programming matrix
	int ** F = new int * [ L2+1 ];
	for( int i = 0; i <= L2; i++ ){
		F[ i ] = new int [ L1+1 ];
	}

	// Traceback matrix
	char ** traceback = new char * [ L2+1 ];
	for( int i = 0; i <= L2; i++ ){
		traceback[ i ] = new char [ L1+1];
	}

	// Initialize traceback and F matrix (fill in first row and column)
	dpm_init( F, traceback, L1, L2, d );

	// Create alignment
//	cout<<"nwalign"<<endl;
	nw_align( F, traceback, seq_1, seq_2, seq_1_al, seq_2_al, d );
//	cout<<"end"<<endl;
#if 0
	int  L_al = seq_1_al.length();
	cout << "Length after alignment: " << L_al << endl;
#endif

	if( prm )
	{
		cout << "\nDynamic programming matrix: " << "\n\n";
		print_matrix( F, seq_1, seq_2 );

		cout << "\nTraceback matrix: " << "\n\n";
		print_traceback( traceback, seq_1, seq_2 );

		cout << endl;
	}

	for( int i = 0; i <= L2; i++ )  delete F[ i ];
	delete[] F;
	for( int i = 0; i <= L2; i++ )  delete traceback[ i ];
	delete[] traceback;

	return  0 ;
}


void  dpm_init( int ** F, char ** traceback, int L1, int L2, int d )
{
	F[ 0 ][ 0 ] =  0 ;
	traceback[ 0 ][ 0 ] = 'n' ;

	int i=0, j=0;

	for( j = 1; j <= L1; j++ )
	{
//		F[ 0 ][ j ] =  -j * d ;
		F[ 0 ][ j ] =  0 ;
		traceback[ 0 ][ j ] =  '-' ;
	}
	for( i = 1; i <= L2; i++ )
	{
//		F[ i ][ 0 ] =  -i * d ;
		F[ i ][ 0 ] = 0 ;
		traceback[ i ][ 0 ] =  '|' ;
	}
}


int nw_align(int ** F,char ** traceback,const string& seq_1, const string& seq_2,string& seq_1_al, string& seq_2_al, int d/* Gap penalty*/){
	int        k = 0, x = 0, y = 0;
	int        fU, fD, fL ;
	char       ptr, nuc ;
	int        i = 0, j = 0;

	const int  a =  1;   // Match
	const int  b = -1;   // Mismatch

	const int  s[ 4 ][ 4 ] =
	  { { a, b, b, b },    /* substitution matrix */
		{ b, a, b, b },
		{ b, b, a, b },
		{ b, b, b, a } } ;

	int  L1 = (int)seq_1.length();
	int  L2 = (int)seq_2.length();

	for( i = 1; i <= L2; i++ )
	{
		for( j = 1; j <= L1; j++ )
		{
			nuc = seq_1[ j-1 ] ;

			switch( nuc )
			{
				case 'A':  x = 0 ;  break ;
				case 'C':  x = 1 ;  break ;
				case 'G':  x = 2 ;  break ;
				case 'T':  x = 3 ;
			}

			nuc = seq_2[ i-1 ] ;

			switch( nuc )
			{
				case 'A':  y = 0 ;  break ;
				case 'C':  y = 1 ;  break ;
				case 'G':  y = 2 ;  break ;
				case 'T':  y = 3 ;
			}

			fU = F[ i-1 ][ j ] - d ;
			fD = F[ i-1 ][ j-1 ] + s[ x ][ y ] ;
			fL = F[ i ][ j-1 ] - d ;

			F[ i ][ j ] = max( fU, fD, fL, &ptr ) ;

			traceback[ i ][ j ] =  ptr ;
		}
	}
	i-- ; j-- ;
//	cout<<"while"<<endl;
	while( i > 0 || j > 0 )
	{
		switch( traceback[ i ][ j ] )
		{
			case '|' :      seq_1_al += '-' ;
				seq_2_al += seq_2[ i-1 ] ;
				i-- ;
				break ;

			case '\\':      seq_1_al += seq_1[ j-1 ] ;
				seq_2_al += seq_2[ i-1 ] ;
				i-- ;  j-- ;
				break ;

			case '-' :      seq_1_al += seq_1[ j-1 ] ;
				seq_2_al += '-' ;
				j-- ;
		}
//		cout<<traceback[ i ][ j ]<<endl;
//		cin.get();
		k++ ;
	}
//	cout<<"end"<<endl;

//	reverse( seq_1_al.begin(), seq_1_al.end() );
//	reverse( seq_2_al.begin(), seq_2_al.end() );

	return  0 ;
}


int  max( int f1, int f2, int f3, char * ptr )
{
	int  max = 0 ;

	if( f1 >= f2 && f1 >= f3 )
	{
		max = f1 ;
		*ptr = '|' ;
	}
	else if( f2 > f3 )
	{
		max = f2 ;
		*ptr = '\\' ;
	}
	else
	{
		max = f3 ;
		*ptr = '-' ;
	}

	return  max ;
}


void  print_matrix( int ** F, const string& seq_1, const string& seq_2 )
{
	int  L1 = (int)seq_1.length();
	int  L2 = (int)seq_2.length();

	cout << "        ";
	for( int j = 0; j < L1; j++ )
	{
		cout << seq_1[ j ] << "   ";
	}
	cout << "\n  ";

	for( int i = 0; i <= L2; i++ )
	{
		if( i > 0 )
		{
			cout << seq_2[ i-1 ] << " ";
		}
		for( int j = 0; j <= L1; j++ )
		{
			cout.width( 3 );
			cout << F[ i ][ j ] << " ";
		}
		cout << endl;
	}
}


void  print_traceback( char ** traceback, const string& seq_1, const string& seq_2 )
{
	int  L1 = (int)seq_1.length();
	int  L2 = (int)seq_2.length();

	cout << "    ";
	for( int j = 0; j < L1; j++ )
	{
		cout << seq_1[ j ] << " ";
	}
	cout << "\n  ";

	for( int i = 0; i <= L2; i++ )
	{
		if( i > 0 )
		{
			cout << seq_2[ i-1 ] << " ";
		}
		for( int j = 0; j <= L1; j++ )
		{
			cout << traceback[ i ][ j ] << " ";
		}
		cout << endl;
	}
}


void  print_al( string& seq_1_al, string& seq_2_al )
{
	cout << seq_1_al << endl;
	cout << seq_2_al << endl;
}


