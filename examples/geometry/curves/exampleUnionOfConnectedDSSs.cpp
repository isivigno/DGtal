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
#include "DGtal/arithmetic/IntegerComputer.h"
#include "DGtal/geometry/curves/ArithmeticalDSSComputer.h"
#include "DGtal/geometry/curves/ArithDSSIterator.h"
#include "DGtal/geometry/curves/StabbingLineComputer.h"
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
using namespace Z2i;

//#define TRACE

/**
 * @brief Function that illustrates the basic usage of
 * a naive DSS. 
 */

Integer determinant(Point u, Point v)
{
  return (u[0]*v[1]-u[1]*v[0]);
}


NaiveDSS8<Integer> computeUnionOfDSSs(NaiveDSS8<Integer> S1, NaiveDSS8<Integer> S2)
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
 
 
  // S1 = NaiveDSS8<Integer>( 1, 1,                   //slope
  // 			 Point(5,6), Point(10,11), //ending points 
  // 			 Point(5,6), Point(10,11), //upper points
  // 			 Point(5,6), Point(10,11)  //lower points
  // 			 );
  // S2 = NaiveDSS8<Integer>( 1, 1,                   //slope
  // 			 Point(11,11), Point(25,25), //ending points 
  // 			 Point(11,11), Point(25,25), //upper points
  // 			 Point(11,11), Point(25,25)  //lower points
  // 			 );
 
 

  
  // Trace to the standard output
  //trace.info() << S1 << "\n" << S2 << std::endl << std::endl; 

  NaiveDSS8<Integer> DSSres(S1);
  // Check if one DSS can be directly added to the other one
  bool easyUnion = false;
  if(S1.b() < S2.b())
    { 
      if(S2.isInDSL(S1.Uf()) && S2.isInDSL(S1.Ul()) && S2.isInDSL(S1.Lf()) && S2.isInDSL(S1.Ll()))
	{
	  easyUnion = true;
	  DSSres = S2;
	}
    }
  else
    {
      if(S1.isInDSL(S2.Uf()) && S1.isInDSL(S2.Ul()) && S1.isInDSL(S2.Lf()) && S1.isInDSL(S2.Ll()))
	{
	  easyUnion = true;
	  DSSres = S1;
	}
    }
  
  Point Uf, Ul, Lf, Ll;
  Integer a, b, mu;
  
  bool notADSS = false;
  if(!easyUnion)
    {
      // if slope of S1 is greater than slope of S2
      // then the slope of the result will be less than the slope of S1
      Integer d = determinant(Point(S1.b(), S1.a()),Point(S2.b(),S2.a())); 
      if(d<0 )
	{
	  // Lower leaning points
	  // S1.Lf may be lower leaning point, S1.Ll is definitely not
	  
	  // First potential lower leaning point
	  Lf = S1.Lf();
	  
	  // S2.Ll may be lower leaning point, S2.Lf is not
	  // Second potential lower leaning point
	  Ll = S2.Ll();
	  
	  // Upper leaning points
	  // S1.Ul may be an upper leaning point, S1.Uf is not
	  Uf = S1.Ul();
	  
	  // S2.Uf may be an upper leaning point, S2.Ul is not
	  Ul = S2.Uf();
	  
	}
      else // the slope of S1 is lower than the slope of S2
	if(d>0)
	  {
	    // then the slope of the result will be greater than the slope of S1
	    
	    
	    // Lower leaning points
	    // S1.Ll may be lower leaning point, S1.Lf is definitely not
	    
	    // First potential lower leaning point
	    Lf = S1.Ll();
	    
	    // S2.Lf may be lower leaning point, S2.Ll is not
	    // Second potential lower leaning point
	    Ll = S2.Lf();
	    
	    // Upper leaning points
	    // S1.Uf may be an upper leaning point, S1.Ul is not
	    Uf = S1.Uf();
	    
	    // S2.Ul may be an upper leaning point, S2.Uf is not
	    Ul = S2.Ul();
	  
	  }
	else
	  if(d==0 && S1.mu() < S2.mu())
	    {
	      Lf = S1.Lf();
	      
	      Ll = S2.Lf();
	  
	      Uf = S1.Ul();
	      
	      Ul = S2.Ul();
	    }
	  else // d==0 && S1.mu() > S2.mu()
	    {
	      Lf = S1.Ll();
	      Ll = S2.Ll();
	      Uf = S1.Uf();
	      Ul = S2.Uf();
	    }
      
      // Test whether Lf and Lf are both lower leaning points. Then
      // Uf and Ul must belong to the DSS defined by these lower
      // leaning points, and at least one of them if an upper
      // leaning point.  
      b = Ll[0]-Lf[0];
      a = Ll[1]-Lf[1];
      mu = b-1-a*Lf[0]+b*Lf[1];
      
      Integer r1 = a*Uf[0]-b*Uf[1]+mu;
      Integer r2 = a*Ul[0]-b*Ul[1]+mu;
      if(r1 >= 0 && r1<b && r2>=0 && r2<b)
	{
	  assert(r1==0 || r2==0);
	  if(r1!=0)
	    Uf = Ul;
	  if(r2!=0)
	    Ul = Uf;
	}
      else
	{
	  // Test whether Uf and Uf are both lower leaning points. Then
	  // Lf and Ll must belong to the DSS defined by these upper
	  // leaning points, and at least one of them is a lower
	  // leaning point.  
	  b = Ul[0]-Uf[0];
	  a = Ul[1]-Uf[1];
	  mu = -a*Uf[0]+b*Uf[1];
	  
	  r1 = a*Lf[0]-b*Lf[1]+mu;
	  r2 = a*Ll[0]-b*Ll[1]+mu;
	  
	  if(r1 >= 0 && r1<b && r2>=0 && r2<b)
	    {
	      assert(r1==b-1 || r2==b-1);
	      if(r1!=b-1)
		Lf = Ll;
	      if(r2!=b-1)
		Ll = Lf;
	    }
	  else
	    {
	      notADSS = true;
	      std::cout << "union is not possible\n-------------\n\n";
	    }
	}
    }
  
  //std::cout << Uf << " " << Ul << " " << Lf << " " << Ll << std::endl;
  if(!notADSS && !easyUnion)
    DSSres = NaiveDSS8<Integer>(a,b,S1.back(),S2.front(),Uf,Ul,Lf,Ll);
  
  
  return DSSres;
  
}	

void testUnionOfConnectedDSSs(int modb, int modx, int nbtries)
{
  // Container of digital points
  typedef std::vector<Z2::Point> Container;
  // Iterator on the container
  typedef Container::const_iterator ConstIterator;
  // Input points
  Container contour1, contour2;
  
  typedef NaiveDSS8Computer<ConstIterator> DSSComputer;  
  // Construction of the computer
  DSSComputer theDSSComputer;    
  
  typedef ArithDSSIterator<Integer,8> DSSIterator;
  

  IntegerComputer<Integer> ic;
  
  long double CPUTimeUnion = 0, CPUTimeReco = 0, CPUTimeStabbing = 0;
  
  clock_t timeBegin, timeEnd;
  int nb = 0;
  //timeBegin = clock();
  for ( unsigned int i = 0; i < nbtries; ++i )
    {
      // Pick up the DSL slope
      Integer b( random() % modb + 1 );
      Integer a( random() % b +1);
      
      if ( ic.gcd( a, b ) == 1 )
        {
	  nb++;
          for ( unsigned int j = 0; j < 5; ++j )
            {
	      // Pick up the DSL intercept
              Integer mu = random() % (2*modb);
              //DSL D( a, b, mu );
              
              for (Integer x = 0; x < 10; ++x )
                {
                  // modx modulates the length of the subsegments
		  
		  // Pick up the beginning of the first subsegment
		  Integer x1 = random() % modx;
		  // Pick up the end of the first subsegment
                  Integer x2 = x1 + 1 + ( random() % modx );
                  
		  /************************************************/
		  /* Case of connected subsegments */
		  // The beginning of the second subsegment is set to x2. 
		  Integer x3 = x2+1;
		  
		  // Pick up the end of the second subsegment
		  //Integer x4 = x3 + 1 + ( random() % modx );
		  Integer x4 = x3 + 1 + modx ;
		  
                  //std::cout << a << " " << b << " " << mu << " " << x1 << " " << x2 << std::endl;
                  
                  Integer y1 = ic.floorDiv(a*x1+mu,b);
                  Integer y2 = ic.floorDiv(a*x2+mu,b);
		  Integer y3 = ic.floorDiv(a*x3+mu,b);
		  Integer y4 = ic.floorDiv(a*x4+mu,b);
                  
		  Point A = Point(x1,y1);
                  Point B = Point(x2,y2);
                  Point C = Point(x3,y3);
		  Point D = Point(x4,y4);

// #ifdef TRACE		  
// 		  std::cout << "\n" << A  << " " << B << " " << C << " " << D << std::endl;
// #endif
		  /**************************************************/
		  
		  
		  // Computation of the parameters of the two segments using the recognition algorithm -> could use the DSLsubsegment algorithm
		  DSSIterator  it(a,b,-mu,A);
		  contour1.clear();
		  while((*it)[0] <=x2)
		    {
		      contour1.push_back(*it);
		      ++it;
		    }
		  
		  NaiveDSS8<Integer> DSS1(contour1.begin(),contour1.end());

		  contour2.clear();
		  it = DSSIterator(a,b,-mu,C);
		  while((*it)[0] <=x4)
		    {
		      contour2.push_back(*it);
		      ++it;
		    }
		  NaiveDSS8<Integer> DSS2(contour2.begin(),contour2.end());
		  		 
		  
#ifdef TRACE	  
		  std::cout << "DSS1: " << A << " " << B << " a =" << DSS1.a() << " b =" << DSS1.b() << " mu =" << DSS1.mu() << std::endl << std::endl;
		  std::cout << "DSS2: " << C << " " << D << " a =" << DSS2.a() << " b =" << DSS2.b() << " mu =" << DSS2.mu() << std::endl << std::endl;
#endif		  
                  
		  // Computation of the union of the two segments
		  
		  timeBegin = clock();
		  NaiveDSS8<Integer> DSS = computeUnionOfDSSs(DSS1,DSS2);
		  timeEnd = clock();
		  CPUTimeUnion += ((double)timeEnd-(double)timeBegin)/((double)CLOCKS_PER_SEC);
		  
#ifdef TRACE
		  std::cout << "DSS Union " << DSS << std::endl << std::endl;
		  //std::cout << "---------------------------" << std::endl << std::endl;
#endif
		  // Verification and comparison using the arithmetical recognition algorithm
		  
		  ConstIterator itt = contour2.begin();
		  NaiveDSS8<Integer> DSSGroundTruth = DSS1;
		  // the points of DSS2 are added one by one
		  timeBegin = clock();
		  while(itt != contour2.end())
		    {
		      DSSGroundTruth.extendFront(*itt);
		      ++itt;
		    }
		  timeEnd = clock();
		  
		  assert(DSS == DSSGroundTruth);
		  //std::cout << "ArithDSS " << DSSGroundTruth << std::endl << std::endl;
		  CPUTimeReco += ((double)timeEnd-(double)timeBegin)/((double)CLOCKS_PER_SEC);
		  
		  // Comparison with stabbingLineComputer adding only
		  // DSS2 leaning points
		  
		  typedef pair<Z2::Point,Z2::Point> Pair;
		  typedef std::vector< Pair > PairContainer;
		  // Iterator on the container
		  typedef PairContainer::const_iterator ConstPairIterator;
		  
		  // Build the "contour" as the concatenation of DSS1 and
		  // the leaning points of DSS2
		  PairContainer contour;
		  ConstIterator myIt;
		  myIt = contour1.begin();

		  while(myIt!=contour1.end())
		    {
		      Pair p = make_pair(*myIt+Point(0,1),*myIt);
		      contour.push_back(p);
		      ++myIt;
		    }
		  myIt = contour2.begin();
		  while(myIt!= contour2.end())
		    {
		      if(*myIt == DSS2.Uf() || *myIt == DSS2.Ul() || *myIt == DSS2.Lf() || *myIt == DSS2.Ll())
			{
			  Pair p = make_pair(*myIt+Point(0,1),*myIt);
			  contour.push_back(p);
			}
		      ++myIt;
		    }
		  
		  // Display the "contour"
		  // ConstPairIterator ii;
		  // ii = contour.begin();
		  // while(ii!=contour.end())
		  //   {
		  //     std::cout << *ii << " ";
		  //     ++ii;
		  //   }
		  
		  // std::cout << std::endl << std::endl;
		  
		  StabbingLineComputer<ConstPairIterator> s; 
		  // //extension
		  s.init( contour.begin() );
		  bool ok=true; bool noClock =true;
		  while ( s.end() != contour.end() && ok)
		    {
		      if(((s.end())->first)[0] > x2 && noClock )
			{
			  noClock = false;
			  timeBegin=clock();
			}
		      ok = s.extendFront();
		    }		  
		  
		  timeEnd = clock();
		  // std::cout << s.Uf() << " " << s.Ul() << " " << s.Lf() << " " << s.Ll() << std::endl;
		  //std::cout << "----------------\n";
		  CPUTimeStabbing += ((double)timeEnd-(double)timeBegin)/((double)CLOCKS_PER_SEC);
		  
		  
		  
	
		 

		}
	    }
	}
    }
  std::cout << modx << " " << CPUTimeUnion/(double) nb << " " << CPUTimeReco/(double) nb << " " << CPUTimeStabbing/(double) nb << std::endl;
}

///////////////////////////////////////////////////////////////////////////////

int main( int argc, char** argv )
{
  trace.beginBlock ( "Example exampleUnionOfConnectedDSSs" );
  
  // computeUnionOfDSSs();
  int x;
  for(x = 10;x<10000;x=x*2)
    testUnionOfConnectedDSSs(200,x,1000);
  
  //testUnionOfConnectedDSSs(20,10,1);

  trace.info() << "Args:";
  for ( int i = 0; i < argc; ++i )
    trace.info() << " " << argv[ i ];
  trace.info() << endl;
  
  trace.endBlock();
  return 0;
  
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
