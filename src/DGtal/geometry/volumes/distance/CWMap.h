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
 * @file CWMap.h
 * @author Isabelle Sivignon (\c isabelle.sivignon@gipsa-lab.grenoble-inp.fr )
 * gipsa-lab Grenoble Images Parole Signal Automatique (CNRS, UMR 5216), CNRS, France
 *
 * @date 2015/04/17
 *
 * Header file for module CWMap.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(CWMap_RECURSES)
#error Recursive header files inclusion detected in CWMap.h
#else // defined(CWMap_RECURSES)
/** Prevents recursive inclusion of headers. */
#define CWMap_RECURSES

#if !defined CWMap_h
/** Prevents repeated inclusion of headers. */
#define CWMap_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
#include <utility>
#include <vector>
#include "DGtal/base/CountedPtr.h"
#include "DGtal/base/ConstAlias.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/kernel/domains/HyperRectDomain.h"
#include "DGtal/images/CConstImage.h"
#include "DGtal/kernel/domains/HyperRectDomain.h"

//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // template class CWMap
  /**
   * Description of template class 'CWMap' <p>
   * \brief Aim:
   */
  template < typename TWeightImage,
    typename TCWSeparableMetric,
    typename TImageContainer = 
    ImageContainerBySTLVector<HyperRectDomain<typename TWeightImage::Domain::Space>,
    typename TWeightImage::Domain::Space::Vector> >
    class CWMap
    {
    // ----------------------- Standard services ------------------------------
  public:
  
  BOOST_CONCEPT_ASSERT(( concepts::CConstImage< TWeightImage > ));
  /// TODO  BOOST_CONCEPT_ASSERT(( concepts::CPowerSeparableMetric<TPowerSeparableMetric> ));
  
  ///Copy of the weight image types
  typedef TWeightImage WeightImage;
  typedef typename TWeightImage::Value Weight;
  typedef typename WeightImage::Domain::Space Space;
  typedef typename Space::Vector Vector;
  typedef typename Space::Point Point;
  typedef typename Space::Dimension Dimension;
  typedef typename Space::Size Size;
  typedef typename Space::Point::Coordinate Abscissa;
  
  //ImageContainer::Domain::Space must match with TSpace
  BOOST_STATIC_ASSERT ((boost::is_same< typename TWeightImage::Domain::Space,
			typename TImageContainer::Domain::Space >::value )); 
  
  //ImageContainer value type must be  TSpace::Vector
  BOOST_STATIC_ASSERT ((boost::is_same< typename TWeightImage::Domain::Space::Vector,
			typename TImageContainer::Value >::value )); 
  
  //ImageContainer domain type must be  HyperRectangular
  BOOST_STATIC_ASSERT ((boost::is_same< HyperRectDomain<typename TWeightImage::Domain::Space>,
			typename TImageContainer::Domain >::value )); 
 
  ///Definition of the underlying domain type.
  typedef typename TImageContainer::Domain Domain;
  
  ///We construct the type associated to the separable metric 
  typedef TCWSeparableMetric CWSeparableMetric;
    
    ///Type of resulting image
  typedef TImageContainer OutputImage;
  
  ///Definition of the image model value type.
  typedef Vector Value;
  ///Definition of the image value type.
  typedef typename OutputImage::ConstRange  ConstRange;
  
  ///Self type
  typedef CWMap<TWeightImage, TCWSeparableMetric, TImageContainer> Self;

  /**
   * Constructor.
   * 
   * This constructor computes the Compoundly Weighted Map of a set of point
   * sites using a SeparableMetric metric.  The method associates to
   * each point satisfying the foreground predicate, the closest
   * site for which the predicate is false. 
   *
   * All parameters are aliased in this class.
   *
   * @param aDomain defines the (hyper-rectangular) domain on which
   * the computation is performed.  
   * @param aWeightImage an image
   * returning the weight for some points
   * @param aMetric a CW
   * separable metric instance.
   */
  CWMap(ConstAlias<Domain> aDomain, ConstAlias<WeightImage> aWeightImage, ConstAlias<CWSeparableMetric> aMetric); 
   
  
  

  /** Basic constructor 
   * does nothing
   */
  CWMap(ConstAlias<Domain> aDomain);
  
   
  /**
   * Destructor.
   */
  ~CWMap();

    // ----------------------- Interface --------------------------------------
  public:

    // ------------------- ConstImage model ------------------------

    /**
     * Assignment operator from another CW map.
     *
     *  @param aOtherCWMap another instance of Self
     *  @return a reference to Self
     */
  Self &  operator=(const Self &aOtherCWMap );
  
  /**
   * Returns a reference (const) to the CW map domain.
   *  @return a domain
   */
  const Domain &  domain() const
  {
    return *myDomainPtr;
  }

    
  /**
   * Returns a const range on the CW map values.
   *  @return a const range
   */
  ConstRange constRange() const
  {
    return myImagePtr->constRange();
  }
  
  /**
   * Access to a CW value (a.k.a. vector to the closest site) at a point.
   *
   * @param aPoint the point to probe.
   */
  Value operator()(const Point &aPoint) const
  {
    return myImagePtr->operator()(aPoint);
  }    
  
  /** 
   * @return  Returns the underlying metric.
   */
  const CWSeparableMetric* metricPtr() const
  {
    return myMetricPtr;
  }
  
  /** 
   * @return  Returns the underlying weight image.
   */
  const WeightImage* weightImagePtr() const
  {
    return myWeightImagePtr;
  }


  
  /**
   * Writes/Displays the object on an output stream.
   * @param out the output stream where the object is written.
   */
  void selfDisplay ( std::ostream & out ) const;
  
  void restrictedCWMapBruteForce(ConstAlias<Domain> aDomain, ConstAlias<WeightImage> aWeightImage, ConstAlias<CWSeparableMetric> aMetric);
  
  void CWMapBruteForce(ConstAlias<Domain> aDomain, ConstAlias<WeightImage> aWeightImage, ConstAlias<CWSeparableMetric> aMetric);
  
  
  // ------------------------- Private functions --------------------------------
  private:
  
  /**
   * Compute the CW Map of a set of point sites using a
   * SeparableMetric metric.  The method associates to each point
   * satisfying the foreground predicate, the closest site for which
   * the predicate is false.  
   */
    void compute ( ) ;
  
  /** 
   *  Compute the other steps of the separable CW map.
   * 
   * @param dim the dimension to process
   */    
  void computeOtherSteps(const Dimension dim) const;
 
  /** 
   * Given  a CW map valid at dimension @a dim-1, this method
   * updates the map to make it consistent at dimension @a dim along
   * the 1D span starting at @a row along the dimension @a
   * dim.
   * 
   * @param row starting point of the 1D process.
   * @param dim dimension of the update.
   */
  void computeOtherStep1D (const Point &row, 
			   const Size dim) const;
  


  // ------------------------- Protected methods ------------------------------
  protected:

    /**
     * Constructor.
     * Forbidden by default (protected to avoid g++ warnings).
     */
  //    CWMap();

   // ------------------- Private members ------------------------
  private:

  ///Pointer to the computation domain
    const Domain * myDomainPtr;
  
  ///Copy of the image lower bound
  Point myLowerBoundCopy;
  
  ///Copy of the image lower bound
  Point myUpperBoundCopy;
  
  ///Value to act as a +infinity value
  Point myInfinity;

  protected:
  ///Pointer to the separable metric instance
  const CWSeparableMetric * myMetricPtr;
  
  ///CW map image
  CountedPtr<OutputImage> myImagePtr;
 
  ///Pointer to the point predicate
  const WeightImage * myWeightImagePtr;


  
  }; // end of class CWMap
  

  /**
   * Overloads 'operator<<' for displaying objects of class 'CWMap'.
   * @param out the output stream where the object is written.
   * @param object the object of class 'CWMap' to write.
   * @return the output stream after the writing.
   */
  /* template <typename T> */
  /* std::ostream& */
  /* operator<< ( std::ostream & out, const CWMap<T> & object ); */

} // namespace DGtal


///////////////////////////////////////////////////////////////////////////////
// Includes inline functions.
#include "DGtal/geometry/volumes/distance/CWMap.ih"

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined CWMap_h

#undef CWMap_RECURSES
#endif // else defined(CWMap_RECURSES)
