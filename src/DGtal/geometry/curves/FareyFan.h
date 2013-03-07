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

#pragma once

/**
 * @file FareyFan.h
 * @author Isabelle Sivignon (\c isabelle.sivignon@gipsa-lab.grenoble-inp.fr )
 * gipsa-lab Grenoble Images Parole Signal Automatique (CNRS, UMR 5216), CNRS, France
 *
 * @date 2012/12/11
 *
 * Header file for module FareyFan.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(FareyFan_RECURSES)
#error Recursive header files inclusion detected in FareyFan.h
#else // defined(FareyFan_RECURSES)
/** Prevents recursive inclusion of headers. */
#define FareyFan_RECURSES

#if !defined FareyFan_h
/** Prevents repeated inclusion of headers. */
#define FareyFan_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/kernel/PointVector.h"
#include "DGtal/arithmetic/IntegerComputer.h"
//////////////////////////////////////////////////////////////////////////////

//#define DEBUG

namespace DGtal
{

/////////////////////////////////////////////////////////////////////////////
// class FareyFan
/**
 * Description of class 'FareyFan' <p>
 * \brief Aim:
 */
template<typename TInteger>
class FareyFan
{
    // ----------------------- Standard services ------------------------------
public:
  
  typedef TInteger Integer;
  typedef long double FloatType;
  typedef DGtal::PointVector<2,Integer> Point; // integer point (a,b)
  typedef DGtal::PointVector<3,Integer> PointR; // rational point (p,q,r) = (p/q,r/q)
  typedef DGtal::PointVector<2,Integer> Vector;
  typedef DGtal::PointVector<2,Integer> Rational; // rational number
						  // (p,q) = p/q
  
  /**
     * Destructor.
     */
    ~FareyFan();

    
    typedef enum Position
    {
      BELOW = -1,
      ONTO,
      ABOVE
    } Position;
    

    class Ray
    {
    public :
      Integer myX;
      Integer myY;
      
      /**
       * Default constructor.
       * not valid
       */
      Ray();
      
      /**
       * Constructor with initialisation
       * @param x0, y0 two integers
       */
      Ray(const Integer x0, const Integer y0);
      
      /**
       * Constructor from a rational and a slope
       * @param p a rational point, a slope
       */
      Ray(const PointR p, const Integer slope);
      
      /**
       * Constructor from two rational points
       * @param p,q two rational points
       */
      Ray(const PointR p, const PointR q);
      
      
      /**
       * Equality operator.
       * @param other the object to compare with.
       * @return 'true' either if the leaning points perfectly match
       * or if the first leaning points match to the last ones
       * (same DSS scanned in the reverse way) 
       * and 'false' otherwise
       */
      bool operator==( const Ray & other ) const;
      

      /**
       * Compute the intersection point of this with another ray r.
       * @param r a Ray
       * @return a PointR 
       */
      PointR intersect(Ray r) const;
      
      // Compute the position of point with respect to a ray
      // Return BELOW, ABOVE or ONTO
      Position positionWrtRay(PointR p) const;
      
      ~Ray();
      

    };
    
    class Polygon
    {
    public :
      
      /* List of points of the polygon 
       *
       */
      typedef vector<PointR> PointList;
      PointList myPoints;
      
      
      /* List of rays of the polygon 
       *
       */
      typedef vector<Ray> RayList;
      
      RayList myRays;

      typedef enum { Pin, Qin, Unknown } InFlag;
      
      /**
       * Default constructor.
       * not valid
       */
      Polygon();
      
      /**
       * Constructor with initialisation of the vertices
       * The vertices are supposed to be given from the vertex of smallest abscissa, in counter clockwise order.
       * @param a,b,c,d four rational points
       */
      Polygon(const PointR a,const PointR b,const PointR c,const PointR d);
      
      
      /**
       * Constructor with initialisation of the vertices and the rays
       * @param a,b,c,d four rational points
       * @param rab,rbc,rcd,rda four rays
       */
      Polygon(const PointR a,const PointR b,const PointR c,const PointR d, const Ray rab,const Ray rbc,const Ray rcd,const Ray rda);
      
      /* /\**  */
      /*  * Intersection of the polygon with the edge [A,B] */
      /*  * Assume that the abscissa of A is greater than the left most point of the polygon.  */
      /*  *\/ */
      
      /* void intersect(const PointR A, const PointR B, const Ray rab, PointList *res); */
      
      /* void intersect(const Polygon P, PointList *res); */
      
      /* Affine tranformation of the polygon P according to the matrix 
       * (1 0 / -P[0] 1)(x y) + (0 p[1])
       * Corresponds to a translation of vector (P[0],P[1]) in the
       * primal space if the polygon is seen as a dual polygon.  
       */ 
      void transform(const Point P);
      
      /* Intersection of two convex polygons in a Farey fan for the
	 union of two DSSs: the result is the set of vertices that have
	 neither the maximum nor the minimum abscissa.*/     

      bool convexIntersectForDSSUnion(const Polygon & P, PointList *res) const;
    
      void addPoint(PointList *res, int a, int b, const Polygon P, const Polygon Q) const;

      /* Test whether the segments [ab] and [cd] intersect or not */
      char SegSegTest(PointR a, PointR b, Ray ra, int ia,PointR c, PointR d, Ray rb, int ib) const;
      void  Advance( int *a, int *aa, int n, bool inside,int i, const Polygon P, PointList *res) const;
      InFlag InOut(InFlag inflag, int aHB, int bHA) const;

      PointR directionVector(int i) const;
      //int relativePosition(int i, PointR p) const;
      int relativePosition(Ray r, int i, PointR p) const;
      int AreaSign(PointR A, PointR B, PointR C) const;
      int DotSign(PointR A, PointR B) const;
      
          /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
      
      friend std::ostream&  operator<< ( std::ostream & out, const Polygon & object )
	{
	  object.selfDisplay(out);
	  return out;
	};
      
      void selfDisplay ( std::ostream & out ) const ;
    
      
    bool Between( PointR a, PointR b, PointR c ) const;
    bool lowerThan(Rational a, Rational b) const ; 
    bool greaterThan(Rational a, Rational b) const;
    bool equalTo(Rational a, Rational b) const;
    int forwardEdge(int i) const;
  
      ~Polygon();
            
      
    };
    

    // ----------------------- Interface --------------------------------------
 public:

    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    //void selfDisplay ( std::ostream & out );

    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const;

    // ------------------------- Protected Datas ------------------------------
private:
    // ------------------------- Private Datas --------------------------------
private:

    // ------------------------- Hidden services ------------------------------
protected:

    /**
     * Constructor.
     * Forbidden by default (protected to avoid g++ warnings).
     */
    FareyFan();

private:

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    FareyFan ( const FareyFan & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    FareyFan & operator= ( const FareyFan & other );

    // ------------------------- Internals ------------------------------------
private:
    

 

}; // end of class FareyFan


/**
 * Overloads 'operator<<' for displaying objects of class 'FareyFan'.
 * @param out the output stream where the object is written.
 * @param object the object of class 'FareyFan' to write.
 * @return the output stream after the writing.
 */
// std::ostream&
//  operator<< ( std::ostream & out, const FareyFan & object );
 /* template <typename TInteger> */
 /*   std::ostream& */
 /*   operator<< ( std::ostream & out, const typename FareyFan<TInteger>::Polygon & object ) */
 /*   { */
 /*     object.selfDisplay( out); */
 /*     return out; */
 /*   } */
 /* ; */
  

 
 int forwardEdge(int i);

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#if !defined(BUILD_INLINE)
#include "DGtal/geometry/curves/FareyFan.ih"
#endif


//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined FareyFan_h

#undef FareyFan_RECURSES
#endif // else defined(FareyFan_RECURSES)
