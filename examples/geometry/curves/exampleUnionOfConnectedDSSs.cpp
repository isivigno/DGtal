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
#include "DGtal/arithmetic/SternBrocot.h"
#include "DGtal/geometry/curves/ArithmeticalDSSComputer.h"
#include "DGtal/geometry/curves/ArithDSSIterator.h"
#include "DGtal/geometry/curves/StabbingLineComputer.h"
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
using namespace Z2i;

//#define DEBUG

#define NBINS 1000

Integer determinant(Point u, Point v)
{
  return (u[0]*v[1]-u[1]*v[0]);
}

Vector highestFractionInBetween(Vector v, Vector w)
{
  
  typedef SternBrocot<Integer, Integer > SB; // type of Stern-Brocot tree
  typedef SB::Fraction Fraction; // type of fraction
  typedef Fraction::ConstIterator ConstIterator; // iterator type for visiting quotients
  typedef Fraction::Value Value; // value of the iterator, a pair (quotient, depth)
    
  Fraction f(v[1],v[0]);
  
  Fraction g(w[1],w[0]); 
 
  // if(f>g)
  //   {
  //     Fraction tmp = f;
  //     f = g;
  //     g = tmp;
  //   }
  
   
  //std::cout << endl;

  int depthf = f.k();
  int depthg = g.k();
  
  //std::cout << "depths = " << depthf << " " << depthg << std::endl;
  
  //ConstIterator fbegin = f.begin(), fend = f.end();
  //ConstIterator gbegin = g.begin(), gend = g.end();
  
  Fraction res;
  
  bool found=false;
  int i = 0;
  bool minf;
  
  bool noRes = false;

  ConstIterator itf=f.begin(), itg=g.begin();
  for(itf = f.begin(), itg = g.begin(); itf != f.end() && itg != g.end() && !found && !noRes; ++itf, ++itg)
    {
      Value uf = *itf;
      Value ug = *itg;
      if(uf.first == ug.first)
	{
	  res.push_back(std::make_pair(uf.first,i));
	  i++;
	}
      else
	{
	  int coeff = min(uf.first,ug.first);
	  if(i==f.k() && i == g.k())
	    {
	      if(abs(uf.first-ug.first) <= 1)
		{
		  //std::cout << "No fraction of smaller denominator between f and g.\n";
		  noRes = true;
		  //exit(0);
		  //		  return 0;
		}
	      else
		{
		  res.push_back(std::make_pair(coeff+1,i));
		  //trace.info() << "found : fractions are ancestors but.." << std::endl;
		  found = true;
		}
	    }	
	  else
	    if((i!=f.k() && i !=g.k()) || (i==f.k() && uf.first > 2) || (i==g.k() && ug.first > 2))
	      {
		res.push_back(std::make_pair(coeff+1,i));
		//trace.info() << "found : fractions are not ancestors" << std::endl;
		found = true;
	      }	
	    else // we have itf++ = fend or itg++ = gend and exit the loop
	      res.push_back(std::make_pair(coeff,i));
	  
	}
    }
  
  if(!found && !noRes) // one fraction is the ancestor of the other
    {
      if(itf == f.end()) // either itf == fend, or itg == gend, but not both of them 
	minf = true;
      else
	minf = false;
      
      
      //std::cout << "ancestor" << std::endl;
      if(minf) // f is an ancestor of g
	{
	  //std::cout << "f ancestor of g" << std::endl;
	  if(g.k() == f.k() + 1)
	    res.push_back(std::make_pair(1+1,f.k()+1));
	  else
	    res.push_back(std::make_pair(1,f.k()+1));
	  ++itg;
	  if(itg!=g.end())
	    {
	      if(g.k() == f.k() +2) // no fraction of smaller denominator between f and g
		{
		  //std::cout << "No fraction of smaller denominator between f and g.\n";
		  noRes = true;
		  //res = f;
		  //exit(0);
		  //return 0;
		}
	      else
		res.push_back(std::make_pair((*itg).first+1,f.k()+2));
	    }
	  
	}
      else // g is an ancestor of f
	{
	  //std::cout << "g ancestor of f" << std::endl;
	  if(f.k() == g.k() + 1)
	    res.push_back(std::make_pair(1+1,g.k()+1));
	  else
	    res.push_back(std::make_pair(1,g.k()+1));
	  ++itf;
	  if(itf!=f.end())
	    {
	      //std::cout << (*itf).first;
	      if(f.k() == g.k() +2) // no fraction of smaller denominator between f and g
		{
		  //std::cout << "No fraction of smaller denominator between f and g.\n";
		  noRes = true;
		  //res = g;
		  //exit(0);
		  //return 0;
		}
	      else
		res.push_back(std::make_pair((*itf).first+1,g.k()+2));
	    }
	  
	}
    }
  
  if(noRes || res == f || res ==g)
    res = Fraction(f.p()+g.p(),f.q()+g.q());
  
  //std::cout << "res = " << res.p() << "/" << res.q() << std::endl;
  
  return Vector(res.q(),res.p());
  
}

// Does not work when the two DSS are not connected
// Very slightly faster than fastUnionOfDSSs when the DSSs are connected
NaiveDSS8<Integer> ArithmeticalUnionOfDSSs(NaiveDSS8<Integer> S1, NaiveDSS8<Integer> S2)
{
  Point Uf=S1.Uf(), Ul=S1.Ul(), Lf=S1.Lf(), Ll=S1.Ll();
    
  vector<Point> LPoints;
  
  // Test whether the two DSSs are connected or not, or if the last
  // point of S1 and the first point of S2 have the same ordinate. 
  Point A = S1.front();
  Point B = S2.back();
  Vector v = B-A;
  bool connected;
  if(v.normInfinity() > 1 && A[1] != B[1])
    connected = false;
  else
    connected = true;
      
  /// For all the cases, the first step is to compute the critical points

  bool easyUnion = false;
  Point Uff, Ull, Lff, Lll;

  NaiveDSS8<Integer> DSSres(S1);
  
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
  vector<Point>::const_iterator it;
  Integer r;
  
  bool inDSL = true;
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
  if(inDSL)
    return NaiveDSS8<Integer>(S2.a(),S2.b(),S1.back(),S2.front(),Uff,S2.Ul(),Lff,S2.Ll());
  
   
  
  LPoints.clear();
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
  
  
  bool notADSS = false;
  Integer a=S1.a(),b=S1.b(),mu=-(S1.dsl()).mu();
  it = LPoints.begin();
  while(it != LPoints.end() && !notADSS)
    {
      Point p = *it;
      r = a*p[0]-b*p[1] + mu;
           
      if(r>=b)
	{
	  Ll = p;
	  Uf = Ul;
	  b = Ll[0]-Lf[0];
	  a = Ll[1]-Lf[1];
	  mu = b-1-a*Lf[0]+b*Lf[1];
	  if(a*Uf[0]-b*Uf[1]+mu < 0 || a*Uf[0]-b*Uf[1]+mu>b-1)
	    {
	      std::cout << "union is not possible\n";
	      notADSS = true;
	    }
	}
      else
	if(r<0)
	  {
	    Ul = p;
	    Lf = Ll;
	    b = Ul[0]-Uf[0];
	    a = Ul[1]-Uf[1];
	    mu = -a*Uf[0]+b*Uf[1];
	    if(a*Lf[0]-b*Lf[1]+mu > b-1 || a*Lf[0]-b*Lf[1]+mu < 0)
	      {
		std::cout << "union is not possible\n";
		notADSS = true;
	      }
	  }
	else
	  {
	    if(r==0)
	      Ul = p;
	    
	    if(r==b-1)
	      Ll = p;
	  }
      // std::cout << a << " " << b << " " << mu << std::endl; 
      // std::cout << Uf << " " << Ul << " " << Lf <<  " " << Ll << std::endl;
      ++it;
    }
  
  if(!notADSS)
    DSSres = NaiveDSS8<Integer>(a,b,S1.back(),S2.front(),Uf,Ul,Lf,Ll);
  
  return DSSres;
  //std::cout << Uf << " " << Ul << " " << Lf <<  " " << Ll << "\n-----------" << std::endl;

}


NaiveDSS8<Integer> fastUnionOfDSSs(NaiveDSS8<Integer> S1, NaiveDSS8<Integer> S2, Integer *nEasy)
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
 
  // NaiveDSS8<Integer> S1( 1, 4,                   //slope
  // 			 Point(5,2), Point(9,3), //ending points 
  // 			 Point(6,3), Point(6,3), //upper points
  // 			 Point(5,2), Point(9,3)  //lower points
  // 			 );
  // NaiveDSS8<Integer> S2( 1, 5,                   //slope
  // 			 Point(15,4), Point(26,6), //ending points 
  // 			 Point(17,5), Point(22,6), //upper points
  // 			 Point(16,4), Point(21,5)  //lower points
  // 			 );
 

  /**************************************************/
  /**************************************************/
  
    
  // Test whether the two DSSs are connected or not, or if the last
  // point of S1 and the first point of S2 have the same ordinate. 
  Point A = S1.front();
  Point B = S2.back();
  Vector v = B-A;
  bool connected;
  if(v.normInfinity() > 1 && A[1] != B[1])
    connected = false;
  else
    connected = true;
      
  /// For all the cases, the first step is to compute the critical points

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

  // Check the leaning point of S2 with respect to S1
  vector<Point>::const_iterator it;
  Integer r;
  it = LPoints.begin();
  Ull = S1.Ul();
  Lll = S1.Ll();
  //bool below = false;
  //bool above = false;
  while(it != LPoints.end() && inDSL)
    {
      
      r = S1.remainder(*it)-S1.mu();
      
      // if(r>=S1.b())
      // 	below = true;
      // else
      // 	if(r< 0)
      // 	  above = true;
      //else
      
	if(r>= S1.b() || r< 0)
	  inDSL = false;
	else
	  {
	    if(r==0)
	      Ull = *it; // update leaning points in the case where S2
			 // in included in S1
	    
	    if(r==(S1.b()-1))
	      Lll = *it;
	    ++it;
	  }
      //++it;
    }

  // if(above && below)
  //   {
  //     std::cout << "union is not possible\n";
  //     // exit(0);
  //   }
  // else
  //   if(!above && !below)
  //     {
  // 	easyUnion = true;
  // 	//std::cout << "easy" << std::endl;
  // 	//DSSres = S1; // !!!! amend with the computation of the new leaning points: the characteristics do not change, but the leaning points may change. 
  // 	DSSres = NaiveDSS8<Integer>(S1.a(),S1.b(),S1.back(),S2.front(),S1.Uf(),Ull,S1.Lf(),Lll);
  //     }
  
  //   else
      //  {
      
      //std::cout << *it << std::endl;
  if(!inDSL)
    {
      
      //assert(r<0 || r>=S1.b());
      if(r>=S1.b()) // lower exterior point -> the slope decreases
	//	if(below)
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
      (*nEasy)++;
      //std::cout << "** " << S1.Uf() << " " << Ull << " " << S1.Lf() << " " << Lll << std::endl;
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
      //below = false;
      //above = false;
      while(inDSL && it != LPoints.end())
	{
	  r = S2.remainder(*it) - S2.mu();
	  
	  // if(r>=S2.b())
	  //   below = true;
	  // else
	  //   if(r<0)
	  //     above = true;
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
	  //++it;
	}
      
      // if(above && below)
      // 	std::cout << "union is not possible\n";
      // else
      // 	if(!above && !below)
      // 	  {
      // 	    easyUnion = true;
      // 	    DSSres = NaiveDSS8<Integer>(S2.a(),S2.b(),S1.back(),S2.front(),Uff,S2.Ul(),Lff,S2.Ll());
      // 	    //std::cout << "** " << Uff << " " << S2.Ul() << " " << Lff << " " << S2.Ll() << std::endl;
      // 	  }
      // 	else
      if(!inDSL)
	  {
	    // Integer r = S2.remainder(*it);
	    // std::cout << r << std::endl;
	    //assert(r<0 || r>=S2.b());
	    //if(below)
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
	  DSSres = NaiveDSS8<Integer>(S2.a(),S2.b(),S1.back(),S2.front(),Uff,S2.Ul(),Lff,S2.Ll());
	  //std::cout << "** " << Uff << " " << S2.Ul() << " " << Lff << " " << S2.Ll() << std::endl;
	  (*nEasy)++;
	}
    }
  
  bool notADSS = false;
  Integer a,b,mu;
  Integer aa,bb,muu;
  if(!easyUnion)
    {
      //#ifdef DEBUG
      //std::cout << "debug " << Uf << " " << Ul << " " << Lf << " " << Ll << std::endl;
      //#endif 
      
      Integer asol, bsol;
      
      b = Ll[0]-Lf[0];
      a = Ll[1]-Lf[1];
      mu = b-1-a*Lf[0]+b*Lf[1];
      Integer r1 = a*Uf[0]-b*Uf[1]+mu;
      Integer r2 = a*Ul[0]-b*Ul[1]+mu;
      
      bb = Ul[0]-Uf[0];
      aa = Ul[1]-Uf[1];
      muu = -a*Uf[0]+b*Uf[1];
      
      Integer r3 = a*Lf[0]-b*Lf[1]+mu;
      Integer r4 = a*Ll[0]-b*Ll[1]+mu;
      
      bool Lok = (r1 >= 0 && r1<b && r2>=0 && r2<b);
      bool Uok = (r3 >= 0 && r3<bb && r4>=0 && r4<bb);
      
      //test whether Uf, Ul, Lf and Ll are all critical support points
      //std::cout << Lok << " " << Uok << std::endl;
      // if(!Lok && !Uok)
      // 	{
      // 	  std::cout << "Union is not possible\n" << std::endl;
      // 	  notADSS = true;
      // 	}
      // else
	if(!(Lok && Uok))
	  {// if not
	    Integer d = determinant(Ul-Uf, Ll-Lf);
	    if(Uok)
	      {
		asol = aa;
		bsol = bb;
		// if(r3<r4)
		// 	Ll = Lf;
		// else
		// 	Lf = Ll;
		if(d>0)
		  Ll = Lf;
		else
		  Lf = Ll;
	      }
	    else
	      if(Lok)// Lok = true
		{
		  asol = a;
		  bsol = b;
		  if(d>0)
		    Ul = Uf;
		  else
		    Uf = Ul;
		}
	      else
	      	{
	      	std::cout << "Union is not possible\n" << std::endl;
	      	notADSS = true;
	      }
	  }
	else
	  {
	    //  std::cout << "là" << std::endl;
	    //std::cout << "debug " << Uf << " " << Ul << " " << Lf << " " << Ll << std::endl;
	    asol = a;
	    bsol = b;
	  }
      if(!notADSS && connected)
	DSSres = NaiveDSS8<Integer>(asol,bsol,S1.back(),S2.front(),Uf,Ul,Lf,Ll);
      
    }
  
    //   if(a == aa && b == bb) // then it's done, Uf, Ul, Lf, Ll are all leaning points
    //   	{
    //   	  DSSres = NaiveDSS8<Integer>(a,b,S1.back(),S2.front(),Uf,Ul,Lf,Ll);
    //   	}
    //   else
    // 	{
	  
    // 	  // Test whether Lf and Lf are both lower leaning points. Then
    // 	  // Uf and Ul must belong to the DSS defined by these lower
    // 	  // leaning points, and at least one of them is an upper
    // 	  // leaning point in the connected case 
	  
    // 	  mu = b-1-a*Lf[0]+b*Lf[1];
    // 	  Integer r1 = a*Uf[0]-b*Uf[1]+mu;
    // 	  Integer r2 = a*Ul[0]-b*Ul[1]+mu;
	  
    // 	  notADSS = false;
	  	  
    // 	  if(r1 >= 0 && r1<b && r2>=0 && r2<b) // Uf and Ul belong to the DSS defined by Lf and Ll
    // 	    {
    // 	      //assert(r1==0 || r2==0);
    // 	      // if(r1!=0 && r2 !=0)
    // 	      // std::cout << "pb with r1 r2 " << a << " " << b << " " << Uf << " " << r1 << " " << Ul << " " << r2 << std::endl;
    // 	      if(connected)
    // 		{
    // 		  if(r1!=0)
    // 		    Uf = Ul;
    // 		  if(r2!=0)
    // 		    Ul = Uf;
    // 		}
    // 	    }
    // 	  else
    // 	    {
	      
    // 	      std::cout << "ici " << r1 << " " << r2<< std::endl;
    // 	      std::cout << S1 << "\n" << S2 << std::endl;
    // 	      // Test whether Uf and Uf are both upper leaning points. Then
    // 	      // Lf and Ll must belong to the DSS defined by these upper
    // 	      // leaning points, and at least one of them is a lower
    // 	      // leaning point.  
    // 	      b = Ul[0]-Uf[0];
    // 	      a = Ul[1]-Uf[1];
	      
    // 	      mu = -a*Uf[0]+b*Uf[1];
	      
    // 	      r1 = a*Lf[0]-b*Lf[1]+mu;
    // 	      r2 = a*Ll[0]-b*Ll[1]+mu;
	      
    // 	      if(r1 >= 0 && r1<b && r2>=0 && r2<b)
    // 	      {
    // 		//assert(r1==b-1 || r2==b-1);
    // 		  if(connected)
    // 		    {
    // 		      if(r1!=b-1)
    // 			Lf = Ll;
    // 		      if(r2!=b-1)
    // 			Ll = Lf;
    // 		    }
    // 	      }
    // 	      else
    // 	      	{
    // 	      	  notADSS = true;
		  
    // 	      	  std::cout << S1 << "\n" << S2 << std::endl;
    // 	      	  std::cout << "union is not possible\n-------------\n\n";
    // 	      	}
    // 	    }
	  
	  
    //   // std::cout << a << " " << b << " " << mu << " " << r1 << " " << r2 << std::endl;
	 
    // 	  if(!notADSS && !easyUnion)
    // 	    {
    // 	      DSSres = NaiveDSS8<Integer>(a,b,S1.back(),S2.front(),Uf,Ul,Lf,Ll);
    // 	      //std::cout << "** " << Uf << " " << Ul << " " << Lf << " " << Ll << std::endl;
    // 	    }
    // 	}
    // }
  //std::cout << S1 << "\n" << S2 << std::endl;
  
  // If the two DSSs are connected or if the last point of S1 and the
  // first of S2 have the same y-coordinate, we are done, DSSres is
  // indeed the DSS of minimal characteristics.   
  //////////////////////////////////////
  /// Otherwise, some extra work needs to be done.
  if(!connected && !easyUnion)
    {
      Vector Vup = Ll+Point(0,1) - Uf;
      Vector Vlow = Ul-(Lf+Point(0,1));
      
      Vector res = highestFractionInBetween(Vlow,Vup);
      
      if(Uf == Ul)
	{
	  DSSres = NaiveDSS8<Integer>(res[1],res[0],S1.back(),S2.front(),Uf,Point(0,0),Point(0,0),Point(0,0));
	}
      else
	if(determinant(res,Ul-Uf)>0)
	  DSSres = NaiveDSS8<Integer>(res[1],res[0],S1.back(),S2.front(),Ul,Point(0,0),Point(0,0),Point(0,0));
	else
	  DSSres = NaiveDSS8<Integer>(res[1],res[0],S1.back(),S2.front(),Uf,Point(0,0),Point(0,0),Point(0,0));
    }
  
  //std::cout << DSSres;
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



double timeDis[NBINS];
int nbDis[NBINS];


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
  
  long double CPUTimeUnion = 0, CPUTimeReco = 0, CPUTimeStabbing = 0, CPUTimeAUnion = 0;
  
  Integer nEasy = 0;
  clock_t timeBegin, timeEnd;
  int nb = 0;
  //timeBegin = clock();
  for ( unsigned int i = 0; i < nbtries; ++i )
    {
      // Pick up the DSL slope
      Integer b( random() % modb + 1 );
      Integer a( random() % b +1);
      while(ic.gcd(a,b) !=1)
	a =random()%b +1;

      
      if ( ic.gcd( a, b ) == 1 )
        {
	 
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
                  //Integer x2 = x1 + 1 + (random() % modx);
		  Integer x2 = x1 + modx + (random() % modx);
                  
		  /************************************************/
		  /* Case of connected subsegments */
		  // The beginning of the second subsegment is set to x2. 
		  //Integer x3 = x2+1;
		  
		  /* Case of not-connected subsegments */
		   
		  Integer x3 = x2+2+ ( random() % modx );;
		  
		  
		  // Pick up the end of the second subsegment
		  //Integer x4 = x3 + 1 + ( random() % modx );
		  //Integer x4 = x3 + modx + ( random() % modx ) ;
		  Integer x4 = x3 + modx;
		  
                  //std::cout << a << " " << b << " " << mu << " " << x1 << " " << x2 << std::endl;
                  
                  Integer y1 = ic.floorDiv(a*x1+mu,b);
                  Integer y2 = ic.floorDiv(a*x2+mu,b);
		  Integer y3 = ic.floorDiv(a*x3+mu,b);
		  Integer y4 = ic.floorDiv(a*x4+mu,b);
                  
		  // keep only difficult cases when the two DSS are disconnected
		  // if(y2 != y3)
		    {
		       nb++;
	  
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
		  
		  int n = nEasy;
		  timeBegin = clock();
		  NaiveDSS8<Integer> DSS = fastUnionOfDSSs(DSS1,DSS2, &nEasy);
		  timeEnd = clock();
		  if(n == nEasy)
		    CPUTimeUnion += ((double)timeEnd-(double)timeBegin)/((double)CLOCKS_PER_SEC);
		  if(n==nEasy)
		    {
		      timeDis[(x4-x1)/10] += ((double)timeEnd-(double)timeBegin)/((double)CLOCKS_PER_SEC);
		      nbDis[(x4-x1)/10]++;
			//std::cout << x4-x1 << " " << std::setprecision(15) << (1000*((long double)timeEnd-(long double)timeBegin))/((long double)CLOCKS_PER_SEC) << std::endl;
		      
		    }
		  
		  // timeBegin = clock();
		  // NaiveDSS8<Integer> DSSArith = ArithmeticalUnionOfDSSs(DSS1,DSS2);
		  // timeEnd = clock();
		  // CPUTimeAUnion += ((double)timeEnd-(double)timeBegin)/((double)CLOCKS_PER_SEC);
		  
		  
#ifdef DEBUG
		  std::cout << "DSS Union " << DSS << std::endl << std::endl;
		  //std::cout << "---------------------------" << std::endl << std::endl;
#endif
		  // Verification and comparison using the arithmetical recognition algorithm
		  NaiveDSS8<Integer> DSSGroundTruth = DSS1;
		  if(x3==x2+1)
		    {
		      ConstIterator itt = contour2.begin();
		      
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
		      // if(DSSArith != DSSGroundTruth)
		      // 	{
		      // 	  std::cout << DSS1 << "\n" << DSS2 << std::endl; 
		      // 	  std::cout << DSSArith << std::endl;
		      // 	  std::cout << DSSGroundTruth << std::endl;
		      // 	  std::cout << "error Arith\n-------------" << std::endl;
			  
		      // 	}
		  

		      CPUTimeReco += ((double)timeEnd-(double)timeBegin)/((double)CLOCKS_PER_SEC);
		    }
		  /****************************************************************/
		  // Comparison with stabbingLineComputer adding only 
		  // DSS2 leaning points
		  
		  if(x3==x2+1)
		    {
		  
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
		    }
		    }	  
		    }
	    }
	    
	}
	
    }
    
  //std::cout << modb << " " << modx << " " << nb << " " << CPUTimeAUnion/(double) nb << " " << CPUTimeUnion/(double) nb << " " << CPUTimeReco/(double) nb << " " << CPUTimeStabbing/(double) nb << " " << nEasy << std::endl;
  
  // std::cout << modb << " " << modx << " " << nb << " " << " " << CPUTimeUnion/(double) (nb-nEasy) << " " << CPUTimeReco/(double) nb << " " << CPUTimeStabbing/(double) nb << " " << nEasy << std::endl;
 
}

///////////////////////////////////////////////////////////////////////////////

int main( int argc, char** argv )
{
  trace.beginBlock ( "Example exampleUnionOfConnectedDSSs" );
  
  for(int i=0;i<NBINS;i++)
    {
      timeDis[i] = 0;
      nbDis[i] = 0;
    }
  
  int x,b;
  //for(b = 20; b<5000;b=b*2)
    
  b = 2000;
  for(x = 10;x<2*b;x=x*1.2)
    testUnionOfConnectedDSSs(b,x,200);
  
  for(int i=0;i<NBINS;i++)
    {
      if(nbDis[i] !=0)
	std::cout << i << " " << (double) timeDis[i]/(double) nbDis[i] << std::endl;
    }
  
  //testUnionOfConnectedDSSs(100,50,10);

  // fastUnionOfDSSs();

  trace.info() << "Args:";
  for ( int i = 0; i < argc; ++i )
    trace.info() << " " << argv[ i ];
  trace.info() << endl;
  
  trace.endBlock();
  return 0;
  
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
