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
 * @file testCWMap.cpp
 * @ingroup Tests
 * @author Isabelle Sivignon (\c isabelle.sivignon@gipsa-lab.grenoble-inp.fr )
 * gipsa-lab Grenoble Images Parole Signal Automatique (CNRS, UMR 5216), CNRS, France
 *
 * @date 2015/04/20
 *
 * Functions for testing class CWMap.
 *
 * This file is part of the DGtal library.
 */

///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "DGtal/base/Common.h"
#include "ConfigTest.h"
#include "DGtal/helpers/StdDefs.h"

#include "DGtal/io/boards/Board2D.h"
#include <DGtal/images/ImageContainerBySTLMap.h>
#include "DGtal/shapes/Shapes.h"
#include "DGtal/io/colormaps/GrayscaleColorMap.h"
#include "DGtal/io/colormaps/HueShadeColorMap.h"

#include <DGtal/geometry/volumes/distance/CWMap.h>
#include <DGtal/geometry/volumes/distance/PredicateCWSeparableMetric.h>



///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
using namespace Z2i;

///////////////////////////////////////////////////////////////////////////////
// Functions for testing class CWMap.
///////////////////////////////////////////////////////////////////////////////


bool testCWMap()
{
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

///////////////////////////////////////////////////////////////////////////////
// Standard services - public :

int main( int argc, char** argv )
{

  std::cout << DGtal::functions::power<DGtal::int64_t>(2,5) << std::endl;

  trace.beginBlock ( "Testing class CWMap" );
  trace.info() << "Args:";
  for ( int i = 0; i < argc; ++i )
    trace.info() << " " << argv[ i ];
  trace.info() << endl;


  Z2i::Domain domain(Z2i::Point(0,0),Z2i::Point(20,20));
    
  typedef ImageContainerBySTLVector<Z2i::Domain , Z2i::Vector > Image;
  Image image(domain);
 
  // typedef ImageContainerBySTLMap<DigitalSetDomain<Z2i::DigitalSet> , Vector > Image;
  // Image image(setDomain,Vector(0,0));
  
  //Setting some values
  image.setValue(Z2i::Point(7,7), Vector(3,6)); 
  image.setValue(Z2i::Point(10,10), Vector(2,8));
  //image.setValue(Z2i::Point(19,7), Vector(3,5));
  
  typedef PredicateCWSeparableMetric<Space, Z2i::Vector> CWMetric;
  CWMetric CWmetric;
  
  CWMap<Image,CWMetric> map(&domain);

  //map.restrictedCWMapBruteForce(&domain,&image,&CWmetric);
  
  map.CWMapBruteForce(&domain,&image,&CWmetric);


  //CWMap<Image, CWMetric> CWmap(&domain, &image, &CWmetric); // use separable algorithm

  Board2D board;
  
  // Display brute force result
  
  board.clear();

  typedef ImageContainerBySTLVector<Z2i::Domain, DGtal::int64_t> ImageDT;
  ImageDT CWDT(domain);
  
  
  for(CWMap<Image, CWMetric>::Domain::ConstIterator it = map.domain().begin(),
  	itend = map.domain().end(); it != itend; ++it)
    {
      CWMap<Image,CWMetric>::Value site = map( *it );   //closest site to (*it)
      DGtal::int64_t d = Z2i::l2Metric(site,*it);
      CWDT.setValue(*it,d);
      if (site != (*it))
  	Display2DFactory::draw( board,   site - (*it), (*it)); //Draw an arrow
    }

  
  board << domain;
  board.saveEPS("CWMapBruteForce.eps");
  
  board.clear();
  
    //! [DTColormaps]
  //Colormap used for the SVG output
  typedef HueShadeColorMap<long int, 2> HueTwice;
  
  ImageDT::Value maxv2=0;
  //We compute the maximum DT value on the CW map
  for ( ImageDT::ConstRange::ConstIterator it = CWDT.constRange().begin(), itend = CWDT.constRange().end();it != itend; ++it)
    if ( (*it) > maxv2)  maxv2 = (*it);
  
  trace.warning() << CWDT << " maxValue= "<<maxv2<< endl;
  board.clear();
  Display2DFactory::drawImage<HueTwice>(board, CWDT, 0.0, maxv2 + 1);
  board.saveEPS( "CWdt.eps" );
  
  
  
  for(Image::Domain::ConstIterator it = image.domain().begin(); it != image.domain().end() ; it ++)
    {
      if(image(*it) != Vector(0,0))
  	{
  	  Point center = *it;
  	  //board.setFillColor(Color());
  	  board.drawCircle(center[0],center[1],image(*it)[0]);
  	  board.drawCircle(center[0],center[1],image(*it)[1]);
  	}
      
    }

  // board.clear();
  // board.setLineStyle( Board2D::Shape::SolidStyle);

  // for(CWMap<Image, CWMetric>::Domain::ConstIterator it = CWmap.domain().begin(),
  // 	itend = CWmap.domain().end(); it != itend; ++it)
  //   {
  //     CWMap<Image,CWMetric>::Value site = CWmap( *it );   //closest site to (*it)
  //     if (site != (*it))
  // 	Display2DFactory::draw( board,   site - (*it), (*it)); //Draw an arrow
  //   }
  
  
  // // for(Image::Domain::ConstIterator it = image.domain().begin(); it != image.domain().end() ; it ++)
  // //   {
  // //     if(image(*it) != Vector(0,0))
  // // 	{
  // // 	  Point center = *it;
  // // 	  //board.setFillColor(Color());
  // // 	  board.drawCircle(center[0],center[1],image(*it)[0]);
  // // 	  board.drawCircle(center[0],center[1],image(*it)[1]);
  // // 	}
      
  // //   }
  // board << domain;
  
  // board.saveEPS("CWMapSeparable.eps");
  

  
  
  // bool res = testCWMap(); // && ... other tests
  // trace.emphase() << ( res ? "Passed." : "Error." ) << endl;
  // trace.endBlock();
  // return res ? 0 : 1;
  return 0;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
