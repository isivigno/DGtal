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
 * @file CSurfaceLocalGeometricEstimator.h
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2014/04/01
 *
 * Header file for concept CSurfaceLocalGeometricEstimator.cpp
 *
 * This file is part of the DGtal library.
 */

#if defined(CSurfaceLocalGeometricEstimator_RECURSES)
#error Recursive header files inclusion detected in CSurfaceLocalGeometricEstimator.h
#else // defined(CSurfaceLocalGeometricEstimator_RECURSES)
/** Prevents recursive inclusion of headers. */
#define CSurfaceLocalGeometricEstimator_RECURSES

#if !defined CSurfaceLocalGeometricEstimator_h
/** Prevents repeated inclusion of headers. */
#define CSurfaceLocalGeometricEstimator_h

//////////////////////////////////////////////////////////////////////////////
// Inclusions
#include <iostream>
#include "DGtal/base/Common.h"
#include "DGtal/base/CQuantity.h"
#include <boost/iterator/iterator_archetypes.hpp>
//////////////////////////////////////////////////////////////////////////////

namespace DGtal
{

  /////////////////////////////////////////////////////////////////////////////
  // class CSurfaceLocalGeometricEstimator
  /**
     Description of \b concept '\b CSurfaceLocalGeometricEstimator' <p>
     @ingroup Concepts
     @brief Aim: This concept describes an object that can process 
     a range over some (abstract) surface so as to return one estimated quantity for each element 
     of the range (or a given subrange). 

     ### Refinement of 
     - boost::DefaultConstructible,  boost::CopyConstructible, boost::Assignable

     ### Associated types :
     - Surfel : the type of elements of the surface
     - Quantity : the type of the (geometric) quantity that is estimated.

     ### Notation
     - \e X : A type that is a model of CSurfaceLocalGeometricEstimator
     - \e x : object of type \e X
     - \e h : double
     - \e itb, \e ite, \e it,\e itb2, \e ite2 : instances of a model of forward iterators having Surfel as value type.
     - \e ito : an instance of a model of output iterator having Quantity as value type
 
     ### Definitions

     ### Valid expressions and semantics

     | Name           | Expression                    |   | Return type                  | Precondition | Semantics                                 |   | Complexity      |
     |----------------|-------------------------------|---|------------------------------|--------------|-------------------------------------------|---|-----------------|
     | Initialization | x.init( \e h, \e itb, \e ite )         |   | void                         | \e h > 0        | Grid step and range initialization in [\e itb,\e ite)        |   | model dependent |
     | Evaluation     | x.eval( \e it )                  |   | Quantity                     |              | Estimation of the quantity at \e \e it, which must belong to [\e itb,\e ite) |   | model dependent |
     | Evaluation     | \e ito = x.eval( \e itb2, \e ite2, \e ito ) |   | a model of output iterator   |              | Estimation for each element of [\e itb2, \e ite2), which must be a subrange of or the whole range [\e itb,\e ite) |   | model dependent |

     ### Invariants

     ### Models

     - TrueDigitalSurfaceLocalEstimator, VCMDigitalSurfaceEstimator.

     ### Notes
     
     - For now, a CSurfaceLocalGeometricEstimator is not a refinement
       of boost::DefaultConstructible, in opposition to
       CCurveLocalGeometricEstimator.

     @tparam T the type that should be a model of CSurfaceLocalGeometricEstimator.
  */
  template <typename T>
  struct CSurfaceLocalGeometricEstimator
    : boost::DefaultConstructible<T>,  boost::CopyConstructible<T>, boost::Assignable<T>

  {

    // ----------------------- Concept checks ------------------------------
  public:

    typedef typename T::Quantity Quantity;
    BOOST_CONCEPT_ASSERT(( CQuantity< Quantity > ));
    typedef typename T::Surfel Surfel;

    BOOST_CONCEPT_USAGE( CSurfaceLocalGeometricEstimator )
    {
      //init method
      myX.init( myH, myItb, myIte ); 

      ConceptUtils::sameType( myQ, myX.eval( myItb ) );
      ConceptUtils::sameType( myIto, myX.eval( myItb, myIte, myIto ) );
    }

    // ------------------------- Private Datas --------------------------------
  private:
    T myX;
    
    double myH; 
    Quantity myQ;
    boost::iterator_archetype<Surfel,
			      boost::iterator_archetypes::readable_iterator_t,
			      boost::forward_traversal_tag > myItb, myIte; 
    boost::iterator_archetype<Quantity,
			      boost::iterator_archetypes::writable_iterator_t,
			      boost::incrementable_traversal_tag > myIto; 

    // ------------------------- Internals ------------------------------------
  private:

  }; // end of concept CSurfaceLocalGeometricEstimator

} // namespace DGtal

//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#endif // !defined CSurfaceLocalGeometricEstimator_h

#undef CSurfaceLocalGeometricEstimator_RECURSES
#endif // else defined(CSurfaceLocalGeometricEstimator_RECURSES)
