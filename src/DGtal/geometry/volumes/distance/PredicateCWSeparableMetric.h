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
 * @file PredicateCWSeparableMetric.h
 * @author Isabelle Sivignon (\c isabelle.sivignon@gipsa-lab.grenoble-inp.fr )
 * gipsa-lab Grenoble Images Parole Signal Automatique (CNRS, UMR 5216), CNRS, France
 *
 * @date 2015/04/17
 *
 * Header file for module PredicateCWSeparableMetric.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(PredicateCWSeparableMetric_RECURSES)
#error Recursive header files inclusion detected in PredicateCWSeparableMetric.h
#else // defined(PredicateCWSeparableMetric_RECURSES)
/** Prevents recursive inclusion of headers. */
#define PredicateCWSeparableMetric_RECURSES

#if !defined PredicateCWSeparableMetric_h
/** Prevents repeated inclusion of headers. */
#define PredicateCWSeparableMetric_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class PredicateCWSeparableMetric
  /**
   * Description of template class 'PredicateCWSeparableMetric' <p>
   * \brief Aim:
   */
  template <typename TSpace,
            typename TPromoted=DGtal::int64_t>
  class PredicateCWSeparableMetric
  {
    // ----------------------- Standard services ------------------------------
  public:
  
  
  ///Copy the space type
  typedef TSpace Space;
  BOOST_CONCEPT_ASSERT(( concepts::CSpace<TSpace> ));
  
  ///Type for points
  typedef typename Space::Point Point;
  ///Type for points
  typedef typename Point::Coordinate Abscissa;
  ///Type for vectors
  typedef typename Space::Vector Vector;
  
  ///Type for internal distance values
  typedef TPromoted Promoted;
  
  ///Type for internal distance values
  typedef TPromoted Weight;
  BOOST_CONCEPT_ASSERT(( concepts::CInteger<Promoted> ));
  
  ///Type for Value (alias)
  typedef TPromoted Value;
  
  /**
   * Constructor.
   */
  PredicateCWSeparableMetric();

  
  /**
   * Destructor.
   */
  ~PredicateCWSeparableMetric();

  /**
   * Copy constructor.
   */
  PredicateCWSeparableMetric ( const PredicateCWSeparableMetric & /*other*/ ) {}

  /**
   * Assignment.
   * @return a reference on 'this'.
   */
  PredicateCWSeparableMetric & operator= ( const PredicateCWSeparableMetric & /*other*/ ) { return *this;}
  
  /**
   * Destructor.
   */
  ~PredicateCWSeparableMetric();
  
    // ----------------------- Interface --------------------------------------
  public:

  // ----------------------- CPowerMetric --------------------------------------
  
  
  /**
   *  Return the power distance of a point @a aPoint and a weighted
   *  point (@a aQ,@a aWq)
   *
   * @param aPoint a point
   * @param aQ a second point
   * @param aW1q first weight of the second point
   * @param aW2q second weight of the second point (aW1q <= aW2q)
   *
   * @return the power distance between aPoint and (Q,WQ)
   */
  Weight CWDistance(const Point &aPoint,
		    const Point &aQ,
		    const Weight &aW1q, const Weight &aW2q) const;
  
  
  
    /**
     * Given an origin and two points, this method decides which one
     * is closest to the origin. This method should be faster than
     * comparing distance values.
     *
     * @param origin the origin
     * @param first  the first point
     * @param w1F the first point first weight
     * @param w2F the first point second weight
     * @param second the second point
     * @param w1S the second point first weight
     * @param w2S the second point second weight
     *
     * @return a Closest enum: FIRST, SECOND or BOTH.
     */
  DGtal::Closest closestPower(const Point &origin,
			      const Point &first,
			      const Weight &w1F, const Weight &w2F,
			      const Point &second,
			      const Weight &w1S, const Weight &w2S) const;
  
  
  // ----------------------- CPowerSeparableMetric --------------------------------------
  
  /**
     * Given three weighted sites (u,v,w) and a straight segment
     * [startingPoint,endPoint] along dimension dim, we detect if the
     * cells of @a u and @a w @e strictly hide the cell of @a v on the
     * straight line according to the CW metric.
     *
     * This method is in @f$ O(log(n))@f$ if @a n is the size of the
     * straight segment. For @f$ l_2@f$ metric (p=2), the method is in
     * @f$ O(1)@f$.
     *
     * @pre u,v and w must be such that u[dim] < v[dim] < w[dim]
     *
     * @param u a site
     * @param w1u a weight
     * @param w2u a weight
     * @param v a site
     * @param w1v a weight
     * @param w2v a weight
     * @param w a site
     * @param w1w a weight
     * @param w2w a weight
     * @param startingPoint starting point of the segment
     * @param endPoint end point of the segment
     * @param dim direction of the straight line
     *
     * @return true if (u,w) hides v.
     */
    bool hiddenByPower(const Point &u,
		       const Weight &w1u,
		       const Weight &w2u,
		       const Point &v,
		       const Weight &w1v,
		       const Weight &w2v,   
		       const Point &w,
		       const Weight &w1w,
		       const Weight &w2w,
		       const Point &startingPoint,
		       const Point &endPoint,
		       const typename Point::UnsignedComponent dim) const;
  
 



  /**
   * Writes/Displays the object on an output stream.
   * @param out the output stream where the object is written.
   */
  void selfDisplay ( std::ostream & out ) const;
  
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
  PredicateCWSeparableMetric();
  
  private:
  
  /**
   * Copy constructor.
   * @param other the object to clone.
   * Forbidden by default.
   */
  PredicateCWSeparableMetric ( const PredicateCWSeparableMetric & other );
  
  /**
   * Assignment.
   * @param other the object to copy.
   * @return a reference on 'this'.
   * Forbidden by default.
   */
  PredicateCWSeparableMetric & operator= ( const PredicateCWSeparableMetric & other );
  
    // ------------------------- Internals ------------------------------------
  private:

  }; // end of class PredicateCWSeparableMetric


  /**
   * Overloads 'operator<<' for displaying objects of class 'PredicateCWSeparableMetric'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'PredicateCWSeparableMetric' to write.
   * @return the output stream after the writing.
   */
  template <typename T>
  std::ostream&
  operator<< ( std::ostream & out, const PredicateCWSeparableMetric<T> & object );

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "DGtal/geometry/volumes/distance/PredicateCWSeparableMetric.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined PredicateCWSeparableMetric_h

#undef PredicateCWSeparableMetric_RECURSES
#endif // else defined(PredicateCWSeparableMetric_RECURSES)
