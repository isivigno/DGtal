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
 * @file Polygon.h
 * @author Isabelle Sivignon (\c isabelle.sivignon@gipsa-lab.grenoble-inp.fr )
 * gipsa-lab Grenoble Images Parole Signal Automatique (CNRS, UMR 5216), CNRS, France
 *
 * @date 2013/01/28
 *
 * Header file for module Polygon.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(Polygon_RECURSES)
#error Recursive header files inclusion detected in Polygon.h
#else // defined(Polygon_RECURSES)
/** Prevents recursive inclusion of headers. */
#define Polygon_RECURSES

#if !defined Polygon_h
/** Prevents repeated inclusion of headers. */
#define Polygon_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
#include <vector>
//////////////////////////////////////////////////////////////////////////////



namespace DGtal
{

/////////////////////////////////////////////////////////////////////////////
// class Polygon
/**
 * Description of class 'Polygon' <p>
 * \brief Aim:
 */
  template<typename TCoordinate> 
    class Polygon
    {
    // ----------------------- Standard services ------------------------------
 public:
      
      typedef Polygon<TCoordinate> Self;
      typedef TCoordinate Coordinate;
      typedef std::vector<Coordinate> Point;
      typedef std::vector<Point> PointList;
      typedef typename PointList::iterator Iterator;


      /**
       * Default Constructor
       */
      Polygon();
      

    /**
     * Copy constructor.
     * @param other the object to clone.
     * Forbidden by default.
     */
    Polygon ( const Polygon & other );

    /**
     * Assignment.
     * @param other the object to copy.
     * @return a reference on 'this'.
     * Forbidden by default.
     */
    Self & operator= ( const Polygon & other );
    
    /**
     * Add a point to the polygon. Modifies private data myPoints.
     */
    void addPoint(const Point p);
    

    /**
     * Destructor.
     */
    ~Polygon(){};

      private:
    
    PointList myPoints;
    
      public:
    
    /**
     * Accessor
     */
    PointList getPoints();

    
    /**
     * Accessor
     */
    Point getLast();


    void convexIntersect(const Polygon & other, Polygon *Res);


    typedef enum { Pin, Qin, Unknown } InFlag;
    
    
/*
This code is described in "Computational Geometry in C" (Second Edition),
Chapter 7.  It is not written to be comprehensible without the
explanation in that book.

Written by Joseph O'Rourke.
Last modified: December 1997
Questions to orourke@cs.smith.edu.
--------------------------------------------------------------------
This code is Copyright 1997 by Joseph O'Rourke.  It may be freely
redistributed in its entirety provided that this copyright notice is
not removed.
--------------------------------------------------------------------
*/
/* Translated in C++ with template coordinates by I.Sivignon (01/2013) */

    /*---------------------------------------------------------------------
      Function prototypes.
      ---------------------------------------------------------------------*/
    void	PrintSharedSeg( Point p, Point q,Polygon *Res );
    Coordinate  Dot( Point a, Point b );
    int	AreaSign( Point a, Point b, Point c );
    char    SegSegInt( Point a, Point b, Point c, Point d, Point *p, Point *q );
    char    ParallelInt( Point a, Point b, Point c, Point d, Point *p, Point *q );
    bool    Between( Point a, Point b, Point c );
    Point    SubVec( Point a, Point b);
    bool    LeftOn( Point a, Point b, Point c );
    bool    Left( Point a, Point b, Point c );
    bool    Collinear( Point a, Point b, Point c );
    InFlag InOut( Point p, InFlag inflag, int aHB, int bHA, Polygon *Res );
    int     Advance( int a, int *aa, int n, bool inside, Point v, Polygon *Res );





    // ----------------------- Interface --------------------------------------
public:

    /**
     * Writes/Displays the object on an output stream.
     * @param out the output stream where the object is written.
     */
    void selfDisplay ( std::ostream & out );
    
    /**
     * Checks the validity/consistency of the object.
     * @return 'true' if the object is valid, 'false' otherwise.
     */
    bool isValid() const;

    // ------------------------- Protected Datas ------------------------------
private:
    // ------------------------- Private Datas --------------------------------
private:
    
    // ------------------------- Internals ------------------------------------
private:

}; // end of class Polygon


/**
 * Overloads 'operator<<' for displaying objects of class 'Polygon'.
 * @param out the output stream where the object is written.
 * @param object the object of class 'Polygon' to write.
 * @return the output stream after the writing.
 */
template <typename TCoordinate>
  std::ostream&
  operator<< ( std::ostream & out, Polygon<TCoordinate> & object )
  {
    object.selfDisplay( out);
    return out;
  }    
 ;
  

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#if !defined(BUILD_INLINE)
#include "DGtal/geometry/curves/Polygon.ih"
#endif


//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined Polygon_h

#undef Polygon_RECURSES
#endif // else defined(Polygon_RECURSES)
