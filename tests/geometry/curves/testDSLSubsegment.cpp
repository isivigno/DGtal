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
 * @file testDSLSubsegment.cpp
 * @ingroup Tests
 * @author Isabelle Sivignon (\c isabelle.sivignon@gipsa-lab.grenoble-inp.fr )
 * gipsa-lab Grenoble Images Parole Signal Automatique (CNRS, UMR 5216), CNRS, France
 *
 * @date 2012/07/17
 *
 * Functions for testing class DSLSubsegment.
 *
  */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "DGtal/base/Common.h"
///////////////////////////////////////////////////////////////////////////////

#include <map>
#include "DGtal/geometry/curves/DSLSubsegment.h"
#include "DGtal/arithmetic/StandardDSLQ0.h"
#include "DGtal/kernel/CPointPredicate.h"
#include "DGtal/arithmetic/IntegerComputer.h"
#include "DGtal/arithmetic/SternBrocot.h"
#include "DGtal/arithmetic/LighterSternBrocot.h"
#include "DGtal/arithmetic/LightSternBrocot.h"
#include "DGtal/arithmetic/Pattern.h"
#include "DGtal/geometry/curves/ArithDSSIterator.h"
#include "DGtal/geometry/curves/ArithmeticalDSS.h"

//#include <gmpxx.h>

using namespace std;
using namespace DGtal;

///////////////////////////////////////////////////////////////////////////////
// Functions for testing class DSLSubsegment.
///////////////////////////////////////////////////////////////////////////////

//#define CHECK_RES

template <typename Integer>
bool singleTest(Integer a, Integer b, Integer mu, Integer x1, Integer x2)
{
  typedef double Number;
  typedef DGtal::DSLSubsegment<Integer,Integer> DSLSubseg;
  typedef DGtal::DSLSubsegment<Integer,Number> DSLSubsegD;
  
  typedef typename DSLSubseg::Point Point;
  DGtal::IntegerComputer<Integer> ic;

  Integer y1 = ic.floorDiv(a*x1+mu,b);
  Integer y2 = ic.floorDiv(a*x2+mu,b);
  Point A = Point(x1,y1);
  Point B = Point(x2,y2);
  
  DSLSubseg D(a,b,mu,A,B);
  
  Number alpha = (Number) a/(Number) b;
  Number beta = (Number) mu/(Number) b;
  Number precision = (Number) 1/(2*b);

  DSLSubsegD DD(alpha,beta,A,B, precision);

  std::cout << "(" << a << "," << b << "," << mu << ") (" << alpha << "," << beta << "," << precision << ")" << std::endl;
  std::cout  << A << " " << B << std::endl;

  std::cout << "res = " << "(" << D.aa << "," << D.bb << "," << D.Nu << ")" << std::endl;
  std::cout << "res float = " << "(" << DD.aa << "," << DD.bb << "," << DD.Nu << ")" << std::endl;

  assert(D.aa == DD.aa && D.bb == DD.bb && D.Nu == DD.Nu);
  //if(!(D.aa == DD.aa && D.bb == DD.bb && D.Nu == DD.Nu))
  //nberrors++;
		

}



template <typename Integer,typename SmallInteger>
bool testDSLSubsegment( unsigned int nbtries, Integer bb, Integer modx)
{
  //typedef mpf_class Number;
  typedef long double Number;
  typedef DGtal::DSLSubsegment<Integer,Integer> DSLSubseg;
  typedef DGtal::DSLSubsegment<Integer,Number> DSLSubsegD;
  

  typedef ArithDSSIterator<Integer,8> DSSIterator;
  typedef ArithmeticalDSS<DSSIterator,Integer,8> ArithDSS;
  

  typedef typename DSLSubseg::Point Point;
  
  DGtal::IntegerComputer<Integer> ic; 
  
  // Point A(1,5);
  // Point B(6,9);
  // DSLSubseg DD(2,3,15,A,B);

  // std::cout << "aa=" << DD.aa << " bb=" << DD.bb << " Nu=" << DD.Nu << std::endl;
  
  Integer b;
  
  // std::cout << "# a b mu a1 b1 mu1 Ax Ay Bx By" << std::endl;
  
  long double timeTotalSubseg=0,timeTotalSubsegD=0;
  
  clock_t timeBeginSubseg, timeEndSubseg;
  clock_t timeBeginSubsegD, timeEndSubsegD;
  
  int nb = 0;
  int nberrors = 0;
  for ( unsigned int i = 0; i < nbtries; ++i )
    {
      // generate b as a power of 10
      // the parameters of the DSL can be expressed as (a,b,mu) with a,b,mu integers or (a/b,mu/b) as decimal numbers 
      // SmallInteger p(random() % m);
      
      // Integer b = pow(10.0,p);
      
      Integer a( random() % bb);
      
      //std::cout << " a=" <<  a << " b=" << b << std::endl; 
      
      //Number alpha = (Number) a/(Number) bb;
      
      Integer g = ic.gcd(a,bb);
      a = a/g;
      b = bb/g;

      
      Number alpha = (Number) a/b;
      //Number alpha((double) a/b,500);

      if ( ic.gcd( a, b ) == 1 )
        {
	  
          for ( unsigned int j = 0; j < 5; ++j )
            {
              //Integer mu = random() % (2*(Integer) pow(10.0,m));

	      Integer mu = random() % (2*b);
	      
	      Number beta = (Number) mu/(Number) b;
	      //Number beta((double) mu/b,500);
	      
              for (Integer x = 0; x < 10; ++x )
                {
                  nb ++;
		  Integer x1 = random() % modx;
                  Integer x2 = x1 + 1+ random()%modx;
		  //Integer x2 = x1 + 1 + ( random() % modx );
                  
                  std::cout << "(" << a << "," << b << "," << mu << ") (" << alpha << "," << beta << ")" << std::endl;
		  
                  Integer y1 = ic.floorDiv(a*x1+mu,b);
                  Integer y2 = ic.floorDiv(a*x2+mu,b);
                  Point A = Point(x1,y1);
                  Point B = Point(x2,y2);
		  
		  std::cout  << A << " " << B << std::endl;
		  
		  // DSLSubsegment algorithm
		  
		  timeBeginSubseg = clock();
		  DSLSubseg D(a,b,mu,A,B);
		  timeEndSubseg = clock();
		  timeTotalSubseg += ((double)timeEndSubseg-(double)timeBeginSubseg)/(((double)CLOCKS_PER_SEC)/1000);
		  
		  std::cout << "res = " << "(" << D.aa << "," << D.bb << "," << D.Nu << ")" << std::endl;
		  // // DSLSubsegment algorithm using floating points
		  
		  timeBeginSubsegD = clock();
		  //Number precision((double) 1/(5*b),500);
		  Number precision = (double) 1/(2*b);
		  DSLSubsegD DD(alpha,beta,A,B,precision);
		  timeEndSubsegD = clock();
		  //timeTotalSubsegD += ((double)timeEndSubsegD-(double)timeBeginSubsegD)/(((double)CLOCKS_PER_SEC)/1000);
		  std::cout << "res float = " << "(" << DD.aa << "," << DD.bb << "," << DD.Nu << ")" << std::endl;
		  
		  // Compare both results
		  assert(D.aa == DD.aa && D.bb == DD.bb && D.Nu == DD.Nu);
		  if(D.aa == DD.aa && D.bb == DD.bb && D.Nu == DD.Nu)
		    timeTotalSubsegD += ((double)timeEndSubsegD-(double)timeBeginSubsegD)/(((double)CLOCKS_PER_SEC)/1000);
		  else
		    nberrors++;
		  
#ifdef CHECK_RES
		  // Check if the result is ok comparing with ArithmeticalDSS recognition algorithm
		  DSSIterator  it(a,b,-mu,A);
                  ArithDSS myDSS(it);
                  
		  //  timeBeginDSS = clock();
                  while ( (*(myDSS.end()))[0] <=x2 && myDSS.extendForward())
                    {}
		  //timeEndDSS = clock();
		  
                  //std::cout << "a =" << myDSS.getA() << " b =" << myDSS.getB() << " mu =" << myDSS.getMu() << std::endl << std::endl;
		  
		  		  
                  if(D.aa != myDSS.getA() || D.bb != myDSS.getB() || D.Nu != - myDSS.getMu())
		    {
		      std::cout << "ERROR " << std::endl;
		      std::cout << a << " " << b << " " << mu << " " << x1 << " " << x2 << std::endl;    
		      std::cout << "aa=" << D.aa << " bb=" << D.bb << " Nu=" << D.Nu << std::endl;
		      std::cout << "a =" << myDSS.getA() << " b =" << myDSS.getB() << " mu =" << myDSS.getMu() << std::endl << std::endl;
		      assert(D.aa == myDSS.getA() && D.bb == myDSS.getB() && D.Nu == - myDSS.getMu());
		    }
		  
		  //timeTotalDSS += ((double)timeEndDSS-(double)timeBeginDSS)/((double)CLOCKS_PER_SEC)*1000;
#endif CHECK_RES
		  
		  
		  
		  
		}
	      
	    }
	}
    }
  
  std::cout << nb << " " << nberrors ;
  std::cout << " " << (long double) timeTotalSubseg/(nb);
  std::cout << " " << (long double) timeTotalSubsegD/((nb-nberrors));
  
  
  return true;
}






///////////////////////////////////////////////////////////////////////////////
// Standard services - public :

int main( int argc, char** argv )
{
  typedef DGtal::int64_t Integer;
  typedef DGtal::int32_t SmallInteger;
  SmallInteger m = 3; // b = 10^p with p <= m
  
  unsigned int nbtries = ( argc > 1 ) ? atoi( argv[ 1 ] ) :5;
  
  SmallInteger p = 1;
  Integer b;
  
  srand(time(NULL));
  
  // for(p=3;p<=m;p++)
  //   {
  //     b = (Integer) pow(10.0,p);
  //     //std::cout << b << std::endl;
  //     for(Integer modx = 10; modx <= b;modx+=modx/8)
  // 	{
  // 	  std::setprecision(15);
  // 	  std::cout << b << " " << modx << " ";
  // 	  testDSLSubsegment<Integer,SmallInteger>( nbtries, b, modx);
  // 	  std::cout << std::endl;
  // 	}
  //   }
  
  
  // for(p=3;p<=m;p++)
  //   {
  //     b = (Integer) pow(10.0,p);
  //     //std::cout << b << std::endl;
  //     for(Integer modx = 10; modx <= b;modx+=5)
  // 	{
  // 	  std::setprecision(15);
  // 	  std::cout << b << " " << modx << " ";
  // 	  testDSLSubsegment<Integer,SmallInteger>( nbtries, b, modx);
  // 	  std::cout << std::endl;
  // 	}
  //   }
  
  singleTest<Integer>(41,1000,301,722,732);
  singleTest<Integer>(547,1000,930,467,471);

 
  return 1;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
