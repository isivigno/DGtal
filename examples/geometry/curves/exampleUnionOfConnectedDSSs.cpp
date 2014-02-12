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

//#define DEBUG

Integer determinant(Point u, Point v)
{
  return (u[0]*v[1]-u[1]*v[0]);
}


NaiveDSS8<Integer> fastUnionOfDSSs(NaiveDSS8<Integer> S1, NaiveDSS8<Integer> S2)
{

  Point Uf, Ul, Lf, Ll;
  bool inDSL = true;
  
  vector<Point> LPoints;
  

  /*************************************************/
  /*************************************************/
  /***** Simple tests ******************************/
  
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
  // NaiveDSS8<Integer> S2( 2, 5,                   //slope
  // 			 Point(9,3), Point(17,7), //ending points 
  // 			 Point(12,5), Point(17,7), //upper points
  // 			 Point(9,3), Point(14,5)  //lower points
  // 			 );
 
  // NaiveDSS8<Integer> S1( 1, 1,                   //slope
  // 			 Point(5,6), Point(10,11), //ending points 
  // 			 Point(5,6), Point(10,11), //upper points
  // 			 Point(5,6), Point(10,11)  //lower points
  // 			 );
  // NaiveDSS8<Integer> S2( 1, 1,                   //slope
  // 			 Point(11,11), Point(18,18), //ending points 
  // 			 Point(11,11), Point(18,18), //upper points
  // 			 Point(11,11), Point(18,18)  //lower points
  // 			 );
  ///5 6 16 [PointVector] {10, 11} [PointVector] {18, 18} [PointVector] {11, 11} [PointVector] {11, 11}
  
  
  // NaiveDSS8<Integer> S1( 2, 5,                   //slope
  // 			 Point(0,0), Point(9,3), //ending points 
  // 			 Point(0,0), Point(5,2), //upper points
  // 			 Point(2,0), Point(7,2)  //lower points
  // 			 );
  // NaiveDSS8<Integer> S2( 1, 2,                   //slope
  // 			 Point(9,3), Point(14,5), //ending points 
  // 			 Point(11,4), Point(13,5), //upper points
  // 			 Point(10,3), Point(12,4)  //lower points
  // 			 );
  /// res = (3,8,1); [PointVector] {5, 2} [PointVector] {13, 5} [PointVector] {2, 0} [PointVector] {10, 3}
  /**************************************************/
  /**************************************************/

  bool easyUnion = false;
  Point Uff, Ull, Lff, Lll;

  NaiveDSS8<Integer> DSSres(S1);


  if((S2.Uf())[0]< (S2.Lf())[0])
    {
      LPoints.push_back(S2.Uf());
      LPoints.push_back(S2.Lf());
      LPoints.push_back(S2.Ul());
      LPoints.push_back(S2.Ll());
    }
  else
    {
      LPoints.push_back(S2.Lf());
      LPoints.push_back(S2.Uf());
      LPoints.push_back(S2.Ll());
      LPoints.push_back(S2.Ul());
    }
  vector<Point>::const_iterator it;
  Integer r;
  it = LPoints.begin();
  Ull = S1.Ul();
  Lll = S1.Ll();
  while(it != LPoints.end() && inDSL)
    {
      
      r = S1.remainder(*it)-S1.mu();
      if(r>= S1.b() || r< 0)
	inDSL = false;
      else
	{
	  if(r==0)
	    Ull = *it;
	  
	  if(r==(S1.b()-1))
	    Lll = *it;
	  ++it;
	}
    }
  //std::cout << *it << std::endl;
  if(!inDSL)
    {
      //Integer r = S1.remainder(*it);
      //std::cout << r << std::endl;
      assert(r<0 || r>=S1.b());
      if(r>=S1.b()) // lower exterior point -> the slope decreases
	{
	  Uf = S1.Ul();
	  Lf = S1.Lf();
	}
      else // upper exterior point -> the slope increases
	{
	  Uf = S1.Uf();
	  Lf = S1.Ll();
	}
    }
  else
    {
      easyUnion = true;
      //std::cout << "easy" << std::endl;
      //DSSres = S1; // !!!! amend with the computation of the new leaning points: the characteristics do not change, but the leaning points may change. 
      DSSres = NaiveDSS8<Integer>(S1.a(),S1.b(),S1.back(),S2.front(),S1.Uf(),Ull,S1.Lf(),Lll);
    }
  
  /****************************/
  
  if(!easyUnion)
    {
      LPoints.clear();
      if((S1.Ul())[0]< (S1.Ll())[0])
	{
	  LPoints.push_back(S1.Ll());
	  LPoints.push_back(S1.Ul());
	  LPoints.push_back(S1.Lf());
	  LPoints.push_back(S1.Uf());
	}
      else
	{
	  LPoints.push_back(S1.Ul());
	  LPoints.push_back(S1.Ll());
	  LPoints.push_back(S1.Uf());
	  LPoints.push_back(S1.Lf());
	}
      
      inDSL = true;
      Uff = S2.Uf();
      Lff = S2.Lf();
      it = LPoints.begin();
      while(inDSL && it != LPoints.end())
	{
	  r = S2.remainder(*it) - S2.mu();
	  if(r<0 || r>=S2.b())
	    inDSL = false;
	  else
	    {
	 
	      if(r==0)
		Uff = *it;
	      
	      if(r==S2.b()-1)
		Lff = *it;
	      ++it;
	    }
	}
      
      if(!inDSL)
	{
	  // Integer r = S2.remainder(*it);
	  // std::cout << r << std::endl;
	  //assert(r<0 || r>=S2.b());
	  if(r>=S2.b()) // lower exterior point -> the slope is greater than DSS2
	    {
	      Ul = S2.Uf();
	      Ll = S2.Ll();
	    }
	  else // upper exterior point -> the slope is lower than DSS2
	    {
	      Ul = S2.Ul();
	      Ll = S2.Lf();
	    }
	}
      else
	{
	  easyUnion = true;
#ifdef DEBUG
	  std::cout << "easy" << std::endl;
#endif
	  //DSSres = S2; // !!!! amend with the computation of the new leaning points: the characteristics do not change, but the leaning points may change.
	  DSSres = NaiveDSS8<Integer>(S2.a(),S2.b(),S1.back(),S2.front(),Uff,S2.Ul(),Lff,S2.Ll());
	}
    }
  
  bool notADSS;
  Integer a,b,mu;
  Integer aa,bb,muu;
  if(!easyUnion)
    {
#ifdef DEBUG
      std::cout << "debug " << Uf << " " << Ul << " " << Lf << " " << Ll << std::endl;
#endif 
      
      b = Ll[0]-Lf[0];
      a = Ll[1]-Lf[1];

      bb = Ul[0]-Uf[0];
      aa = Ul[1]-Uf[1];
	      
      
      
      if(b == bb && a == aa)
	DSSres = NaiveDSS8<Integer>(a,b,S1.back(),S2.front(),Uf,Ul,Lf,Ll);
      
      // if(a == aa && b == bb) // then it's done, Uf, Ul, Lf, Ll are all leaning points
      // 	{
      // 	  DSSres = NaiveDSS8<Integer>(a,b,S1.back(),S2.front(),Uf,Ul,Lf,Ll);
      // 	}
      else
	{
	  
	  // Test whether Lf and Lf are both lower leaning points. Then
	  // Uf and Ul must belong to the DSS defined by these lower
	  // leaning points, and at least one of them if an upper
	  // leaning point.  
	  
	  mu = b-1-a*Lf[0]+b*Lf[1];
	  Integer r1 = a*Uf[0]-b*Uf[1]+mu;
	  Integer r2 = a*Ul[0]-b*Ul[1]+mu;
	  
	  notADSS = false;
	  if(r1 >= 0 && r1<b && r2>=0 && r2<b) // Uf and Ul belong to the DSS defined by Lf and Ll
	    {
	      assert(r1==0 || r2==0);
	      // if(r1!=0 && r2 !=0)
	      //std::cout << "ni r1 ni r2 !! " << a << " " << b << " " << Uf << " " << r1 << " " << Ul << " " << r2 << std::endl;
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
		  
		  std::cout << S1 << "\n" << S2 << std::endl;
		  std::cout << "union is not possible\n-------------\n\n";
		}
	    }
	  
	  
      // std::cout << a << " " << b << " " << mu << " " << r1 << " " << r2 << std::endl;
      // std::cout << Uf << " " << Ul << " " << Lf << " " << Ll << std::endl;
	  if(!notADSS && !easyUnion)
	    DSSres = NaiveDSS8<Integer>(a,b,S1.back(),S2.front(),Uf,Ul,Lf,Ll);
	}
    }
  
  return DSSres;
  
  

}

template <typename ConstPairIterator>
NaiveDSS8<Integer> computeMinDSSFromPreimage(typename StabbingLineComputer<ConstPairIterator>::PreimagePtr P, Point first, Point last)
{
  typedef typename StabbingLineComputer<ConstPairIterator>::Preimage::Container Hull; 
  Hull pH, qH;
  
  pH = P->pHull();
  qH = P->qHull();
  
  
  Point Uf = P->Uf(); // = pH.rbegin()
  Point Ul = P->Ul(); // = pH.begin()
  Point Lf = P->Lf(); // = qH.rbegin()
  Point Ll = P->Ll(); // = qH.begin()
  
  Point rUf,rUl,rLf,rLl;
  
  // std::cout << Uf << " " << Ul << " " << Lf <<  " " << Ll << std::endl;
  
  Vector Vlow = Lf - Ul;
  Vector Vup = Uf - Ll;
  
  typename Hull::iterator it;
  typename Hull::reverse_iterator rit;
  
  /***********************/
  /// Ul
  
  Point cur, next;
  bool ok = true;
  it = pH.begin();
  cur = *it;
  //std::cout << "first = " << cur;
  while(it != pH.end() && ok)
    {
      if(++it !=pH.end())
	{
	  next = *(it);
	  //std::cout << "next = " << next;
	  if(determinant(next-cur,Vlow)!=0)
	    ok = false;
	  else
	    cur = next;
	}
      else
	ok = false;
    }
  rUl = cur;
  //std::cout << rUl;
  /*************************/
  /// Uf

  ok = true;
  rit = pH.rbegin();
  //  std::cout << "first = " << *rit;
  cur = *rit;
  while(rit != pH.rend() && ok)
    {
      if(++rit !=pH.rend())
	{
	  next = *(rit);
	  //std::cout << "next =" << next;
	  if(determinant(cur-next,Vup)!=0)
	    ok = false;
	  else
	    cur = next;
	}
      else
	ok = false;
    }
  rUf = cur;
  //std::cout << rUf;
  /***********************/
  /// Ll
  
  it = qH.begin();
  cur = *it;
  ok = true;
  while(it != qH.end() && ok)
    {
      if(++it != qH.end())
	{
	  next = *(it);
	  if(determinant(next-cur,Vup)!=0)
	    ok = false;
	  else
	    cur = next;
	}
      else
	ok = false;
    }
  rLl = cur;

  /*************************/
  /// Lf

  ok = true;
  rit = qH.rbegin();
  cur = *rit;
  while(rit != pH.rend() && ok)
    {
      if(++rit != pH.rend())
	{
	  next = *rit;
	  if(determinant(cur-next,Vlow)!=0)
	    ok = false;
	  else
	    cur = next;
	}
      else
	ok = false;
    }
  rLf = cur;

  /****************************/
  
  Integer a, b;
  if(rUl != rUf)
    {
      a = rUf[1]-rUl[1];
      b = rUf[0]-rUl[0];
    }
  else
    {
      a = rLf[1]-rLl[1];
      b = rLf[0]-rLl[0];
    }
  
  IntegerComputer<Integer> ic;
  Integer g = ic.gcd(a,b);
  a = a/g;
  b = b/g;
		 
  NaiveDSS8<Integer> DSS(a,b,first, last, rLl, rLf, rUl+Point(0,-1), rUf+Point(0,-1));
  return DSS;

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
		  		 
		  
#ifdef DEBUG	  
		  std::cout << "DSS1: " << A << " " << B << " a =" << DSS1.a() << " b =" << DSS1.b() << " mu =" << DSS1.mu() << std::endl << std::endl;
		  std::cout << "DSS2: " << C << " " << D << " a =" << DSS2.a() << " b =" << DSS2.b() << " mu =" << DSS2.mu() << std::endl << std::endl;
#endif		  
                  
		  // Computation of the union of the two segments
		  
		  timeBegin = clock();
		  NaiveDSS8<Integer> DSS = fastUnionOfDSSs(DSS1,DSS2);
		  timeEnd = clock();
		  CPUTimeUnion += ((double)timeEnd-(double)timeBegin)/((double)CLOCKS_PER_SEC);
		  
#ifdef DEBUG
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
		  
		  //assert(DSS == DSSGroundTruth);
		  if(DSS != DSSGroundTruth)
		    {
		      std::cout << DSS1 << "\n" << DSS2 << std::endl; 
		      std::cout << DSS << std::endl;
		      std::cout << DSSGroundTruth << std::endl;
		      std::cout << "error\n-------------" << std::endl;
		  
		    }
		  //std::cout << "------------\n";
		  //std::cout << "ArithDSS " << DSSGroundTruth << std::endl << std::endl;
		  CPUTimeReco += ((double)timeEnd-(double)timeBegin)/((double)CLOCKS_PER_SEC);

		  /****************************************************************/
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
		  bool ok=true; bool Clock =false;
		  while ( s.end() != contour.end() && ok)
		    {
		      if(((s.end())->first)[0] > x2 && !Clock )
			{
			  Clock = true;
			  timeBegin=clock();
			}
		      ok = s.extendFront();
		    }

		  
		  NaiveDSS8<Integer> DSSStab(DSS1);
		  
		  typedef StabbingLineComputer<ConstPairIterator>::PreimagePtr myPreimagePtr; 
		  
		  DSSStab = computeMinDSSFromPreimage<ConstPairIterator>(s.getPreimage(),DSS1.back(),DSS2.front());
		  timeEnd = clock();
		  
		  if(DSSStab != DSSGroundTruth)
		    {
		      std::cout << "error DSSStab" << std::endl;
		      std::cout << DSS1 << std::endl;
		      std::cout << DSS2 << std::endl;
		      std::cout << DSSStab << std::endl;
		      std::cout << s << std::endl;
		      std::cout << DSSGroundTruth << std::endl;
		      std::cout << "----------------\n";
		    }
		  // std::cout << s.Uf() << " " << s.Ul() << " " << s.Lf() << " " << s.Ll() << std::endl;
		  //std::cout << "----------------\n";
		  CPUTimeStabbing += ((double)timeEnd-(double)timeBegin)/((double)CLOCKS_PER_SEC);
		  
		  
		  /*****************************************************/
		  /**** Computation with stabbing line computer with all the points ***/
		  
		  // contour.clear();
		  // myIt = contour1.begin();

		  // while(myIt!=contour1.end())
		  //   {
		  //     Pair p = make_pair(*myIt+Point(0,1),*myIt);
		  //     contour.push_back(p);
		  //     ++myIt;
		  //   }
		  // myIt = contour2.begin();
		  // while(myIt!= contour2.end())
		  //   {
		  //     Pair p = make_pair(*myIt+Point(0,1),*myIt);
		  //     contour.push_back(p);
		  //     ++myIt;
		  //   }
		 
		  // StabbingLineComputer<ConstPairIterator> sTot; 
		  // // //extension
		  // sTot.init( contour.begin() );
		  // ok=true; 
		  // while ( sTot.end() != contour.end() && ok)
		  //   {
		  //     ok = sTot.extendFront();
		  //   }
		  
		  // std::cout << sTot << std::endl;
		  
		 

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
  
  // int x;
  // for(x = 10;x<10000;x=x*2)
  //   testUnionOfConnectedDSSs(200,x,1000);
  
  testUnionOfConnectedDSSs(1000,100,100);

  //fastUnionOfDSSs();

  trace.info() << "Args:";
  for ( int i = 0; i < argc; ++i )
    trace.info() << " " << argv[ i ];
  trace.info() << endl;
  
  trace.endBlock();
  return 0;
  
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
