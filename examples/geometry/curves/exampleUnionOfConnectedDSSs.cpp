/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

/**
 * @file exampleUnionOfConnectedDSSs.cpp
 * @ingroup Examples
 * @author siviigni (login) (\c Unknown )
 * Unknown
 *
 * @date 2014/01/31
 *
 * An example file named exampleUnionOfConnectedDSSs.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "ConfigExamples.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/base/Common.h"

//! [ArithmeticalDSSHeader]
#include "DGtal/geometry/curves/ArithmeticalDSS.h"
//! [ArithmeticalDSSHeader]

///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
using namespace Z2i;

/**
 * @brief Function that illustrates the basic usage of
 * a naive DSS. 
 */

Integer determinant(Point u, Point v)
{
  return (u[0]*v[1]-u[1]*v[0]);
}


void computeUnionOfDSSs(NaiveDSS8<Integer> S1, NaiveDSS8<Integer> S2)
{
  //! [ArithmeticalDSSNaiveCtor]
  // Construct a naive DSS
  // NaiveDSS8<Integer> S1( 3, 7,                   //slope
  // 			      Point(0,0), Point(14,6), //ending points 
  // 			      Point(0,0), Point(14,6), //upper points
  // 			      Point(2,0), Point(9,3)  //lower points
  // 			      );
  // //! [ArithmeticalDSSNaiveCtor]
  // NaiveDSS8<Integer> S2( 2, 5,                   //slope
  // 			      Point(15,6), Point(23,9), //ending points 
  // 			      Point(19,8), Point(19,8), //upper points
  // 			      Point(16,6), Point(21,8)  //lower points
  // 			      );
 
  
  // NaiveDSS8<Integer> S1( 2, 5,                   //slope
  // 			 Point(0,0), Point(8,3), //ending points 
  // 			 Point(0,0), Point(5,2), //upper points
  // 			 Point(2,0), Point(7,2)  //lower points
  // 			 );
  // //! [ArithmeticalDSSNaiveCtor]
  // NaiveDSS8<Integer> S2( 3, 7,                   //slope
  // 			 Point(9,3), Point(17,6), //ending points 
  // 			 Point(15,6), Point(15,6), //upper points
  // 			 Point(10,3), Point(17,6)  //lower points
  // 			 );
 
  // NaiveDSS8<Integer> S1( 1, 2,                   //slope
  // 			 Point(0,0), Point(4,2), //ending points 
  // 			 Point(1,1), Point(3,2), //upper points
  // 			 Point(0,0), Point(4,2)  //lower points
  // 			 );
  // //! [ArithmeticalDSSNaiveCtor]
  // NaiveDSS8<Integer> S2( 3, 8,                   //slope
  // 			 Point(5,2), Point(16,7), //ending points 
  // 			 Point(5,2), Point(8,4), //upper points
  // 			 Point(13,5), Point(16,7)  //lower points
  // 			 );
  

  
  // Trace to the standard output
  //trace.info() << S1 << " " << S2 << std::endl; 

  
  // Check if one DSS can be directly added to the other one
  bool easyUnion = false;
  if(S1.b() < S2.b())
    { 
      if(S2.isInDSL(S1.Uf()) && S2.isInDSL(S1.Ul()) && S2.isInDSL(S1.Lf()) && S2.isInDSL(S1.Ll()))
	easyUnion = true;
    }
  else
    {
      if(S1.isInDSL(S2.Uf()) && S1.isInDSL(S2.Ul()) && S1.isInDSL(S2.Lf()) && S1.isInDSL(S2.Ll()))
	easyUnion = true;
    }
  
  Point Uf, Ul, Lf, Ll;
  Integer a, b, mu;
  
  if(!easyUnion)
    {
      // if slope of S1 is greater than slope of S2
      if(determinant(Point(S1.b(), S1.a()),Point(S2.b(),S2.a()))<0)
	{
	  // First lower leaning point
	  Lf = S1.Lf();
	  // Second lower leaning point
	  Integer d = determinant(S2.Lf()-Lf, S2.Ll()-Lf);
	  if(d<=0)
	    Ll = S2.Ll();
	  else
	    Ll = S2.Lf();
	  
	  // Compute characteristics
	  b = Ll[0]-Lf[0];
	  a = Ll[1]-Lf[1];
	  mu = b-1-a*Lf[0]+b*Lf[1];
	  std::cout << a << " " << b << " " << mu << " " << std::endl;
	  std::cout << Lf << " " << Ll << std::endl;
	  
	  // Compute the upper leaning point 
	  Integer r;
	  bool notADSS = false;
	  r = a*(S1.Uf())[0]-b*(S1.Uf())[1]+mu;
	  if(r !=0)
	    {
	      if(r<0 || r>=b)
		notADSS = true;
	      
	      r = a*(S1.Ul())[0]-b*(S1.Ul())[1]+mu;
	      if(r!=0)
		{
		  if(r<0 || r>=b)
		    notADSS = true;
		  
		  r = a*(S2.Uf())[0]-b*(S2.Uf())[1]+mu;
		  if(r!=0)
		    {
		      if(r<0 || r>=b)
			notADSS = true;
		      
		      r = a*(S2.Ul())[0]-b*(S2.Ul())[1]+mu;
		      
		      if(r<0 || r>=b)
			notADSS = true;
		      
		      assert(r==0);
		      Uf = S2.Ul();
		    }
		  else
		    Uf = S2.Uf();
		}
	      else
		Uf = S1.Ul();
	    }
	  else
	    Uf = S1.Uf();
	  
	  Ul = Uf;
	  
	  std::cout << Uf << " " << Ul << std::endl;
	 
	  if(notADSS)
	    std::cout << "the union is not a DSS.\n";
	}
      
      else // the slope of S1 is lower than the slope of S2
	{
	  Uf = S1.Uf();
	  Integer d = determinant(S2.Uf()-Uf, S2.Ul()-Uf);
	  if(d<=0)
	    Ul = S2.Ul();
	  else
	    Ul = S2.Uf();
	  
	  b = Ul[0]-Uf[0];
	  a = Ul[1]-Uf[1];
	  mu = a*Uf[0]+b*Uf[1];
	  std::cout << a << " " << b << " " << mu << " " << std::endl;
	  std::cout << Uf << " " << Ul << std::endl;
	  
	  // Compute the lower leaning point 
	  bool notADSS = false;
	  Integer r;
	  r = a*(S1.Lf())[0]-b*(S1.Lf())[1]+mu;
	  if(r !=b-1)
	    {
	      if(r<0 || r>=b)
		notADSS = true;
	      
	      r = a*(S1.Ll())[0]-b*(S1.Ll())[1]+mu;
	      if(r!=b-1)
		{
		  if(r<0 || r>=b)
		    notADSS = true;
		  
		  r = a*(S2.Lf())[0]-b*(S2.Lf())[1]+mu;
		  if(r!=b-1)
		    {
		      if(r<0 || r>=b)
			notADSS = true;
		      
		      r = a*(S2.Ll())[0]-b*(S2.Ll())[1]+mu;
		      
		      if(r<0 || r>=b)
			notADSS = true;
		      assert(r==b-1);
		      Lf = S2.Ll();
		    }
		  else
		    Lf = S2.Lf();
		}
	      else
		Lf = S1.Ll();
	    }
	  else
	    Lf = S1.Lf();
	  
	  Ll = Lf;
	  
	  std::cout << Lf << " " << Ll << std::endl;
	  if(notADSS)
	    std::cout << "the union is not a DSS.\n";
	}
      
    }
  
}

void testUnionOfConnectedDSSs()
{

}


///////////////////////////////////////////////////////////////////////////////

int main( int argc, char** argv )
{
  trace.beginBlock ( "Example exampleUnionOfConnectedDSSs" );
  
  computeUnionOfDSSs();
  
  trace.info() << "Args:";
  for ( int i = 0; i < argc; ++i )
    trace.info() << " " << argv[ i ];
  trace.info() << endl;

  trace.endBlock();
  return 0;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
