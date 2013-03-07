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
 * @file testFareyFan.cpp
 * @ingroup Tests
 * @author Isabelle Sivignon (\c isabelle.sivignon@gipsa-lab.grenoble-inp.fr )
 * gipsa-lab Grenoble Images Parole Signal Automatique (CNRS, UMR 5216), CNRS, France
 *
 * @date 2012/12/11
 *
 * Functions for testing class FareyFan.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/geometry/curves/FareyFan.h"
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;

///////////////////////////////////////////////////////////////////////////////
// Functions for testing class FareyFan.
///////////////////////////////////////////////////////////////////////////////
/**
 * Example of a test. To be completed.
 *
 */
bool testFareyFan()
{
  typedef DGtal::FareyFan<DGtal::int32_t> Fan;
  typedef Fan::Ray Ray;
  typedef Fan::PointR PointR;
  typedef Fan::Polygon Polygon;
  typedef Fan::Polygon::PointList PointList;

  Ray r1(6,4);
  Ray r2(1,1);
  PointR p1(3,5,1), p2(3,5,2), p3(3,5,3);
  
  std::cout << "Ray 1: " << r1.myX << " " << r1.myY << std::endl;
  std::cout << "Ray 2: " << r2.myX << " " << r2.myY << std::endl;

  PointR res = r1.intersect(r2);

  std::cout << "Intersection = (" << res[0] << "/" << res[1] << "," << res[2] << "/" << res[1] << ")" << std::endl;

  std::cout << "Position:" << "p1 ->" << r1.positionWrtRay(p1) << std::endl;
  std::cout << "Position:" << "p2 ->" << r1.positionWrtRay(p2) << std::endl;
  std::cout << "Position:" << "p3 ->" << r1.positionWrtRay(p3) << std::endl;

  PointR A(1,4,3), B(1,3,1), C(1,2,0), D(1,3,2);

  Polygon P(A,B,C,D);
  
  std::cout << P.myPoints[0]  << " (" << P.myRays[0].myX << "," << P.myRays[0].myY << ") " << P.myPoints[1] << std::endl;


  // PointR AA(1,4,4), BB(3,8,0), CC(3,7,0), DD(2,7,7);
  // Ray r(AA,BB);

  // std::cout << r.x << " " << r.y << std::endl;

  // PointList l;

  // // P.intersect(AA,BB,r,&l);

  // PointList::iterator it;

  // // for(it = l.begin();it != l.end();it++)
  // //   {
  // //     std::cout << *it;

  // //   }
  
  // Polygon Q(AA,BB,CC,DD);
  
  // l.clear();
  //  P.intersect(Q,&l);
  
  // for(it = l.begin();it != l.end();it++)
  //   {
  //     std::cout << *it;

  //   }


  Fan::Point pt(4,3);
  
  P.transform(pt);
  
  std::cout << P.myPoints[0]  << " (" << P.myRays[0].myX << "," << P.myRays[0].myY << ") " << P.myPoints[1] << std::endl;


  unsigned int nbok = 0;
  unsigned int nb = 0;
  
  trace.beginBlock ( "Testing block ..." );
  nbok += true ? 1 : 0; 
  nb++;
  trace.info() << "(" << nbok << "/" << nb << ") "
	       << "true == true" << std::endl;
  trace.endBlock();
  
  return nbok == nb;
}

bool testConvexIntersect()
{
  typedef DGtal::FareyFan<DGtal::int32_t> Fan;
  typedef Fan::Ray Ray;
  typedef Fan::PointR PointR;
  typedef Fan::Polygon Polygon;
  typedef Fan::Polygon::PointList PointList;

  //PointR A(0,1,1), B(1,3,0), C(1,2,0),D(1,3,1);
  //PointR AA(0,1,2), BB(1,2,-1), CC(1,1,-3),DD(1,2,0);
  
  PointR A(0,1,3), B(1,2,4), C(1,1,2),D(1,2,5);
  PointR AA(0,1,4), BB(1,2,5), CC(1,1,1),DD(1,2,6);
  

  
  Polygon P1(A,B,C,D);
  Polygon P2(AA,BB,CC,DD);
  
  std::cout << P1 << std::endl;
  std::cout << P2 << std::endl;

  PointList res;
  
  PointList::iterator it;
  
  P1.convexIntersectForDSSUnion(P2, &res);
  
  for(it=res.begin();it !=res.end() ; it++)
    std::cout << (*it) << std::endl;

  return true;
  
}



///////////////////////////////////////////////////////////////////////////////
// Standard services - public :

int main( int argc, char** argv )
{
  trace.beginBlock ( "Testing class FareyFan" );
  trace.info() << "Args:";
  for ( int i = 0; i < argc; ++i )
    trace.info() << " " << argv[ i ];
  trace.info() << endl;

  bool res = testFareyFan() && testConvexIntersect(); // && ... other tests
  trace.emphase() << ( res ? "Passed." : "Error." ) << endl;
  trace.endBlock();
  return res ? 0 : 1;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
