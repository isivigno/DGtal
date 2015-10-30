#include <iostream>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/boards/Board2D.h"
#include <DGtal/images/ImageContainerBySTLMap.h>
#include "DGtal/shapes/Shapes.h"

#include "DGtal/io/colormaps/GrayscaleColorMap.h"

#include <DGtal/geometry/volumes/distance/CWMap.h>
#include <DGtal/geometry/volumes/distance/PredicateCWSeparableMetric.h>

using namespace std;
using namespace DGtal;
using namespace Z2i;

int main(int argc, char ** argv)
{

trace.beginBlock( "Testing CWDistanceTransform ..." );

Z2i::Domain domain(Z2i::Point(0,0),Z2i::Point(25,25));

Z2i::DigitalSet set(domain);
set.insertNew(Z2i::Point(7,7)); 
//set.insertNew(Z2i::Point(3,7)); 
set.insertNew(Z2i::Point(12,12));
DigitalSetDomain<Z2i::DigitalSet> setDomain(set); 

  typedef ImageContainerBySTLVector<DigitalSetDomain<Z2i::DigitalSet> , Vector > Image;
  Image image(setDomain,Vector(0,0));
  
  //Setting some values
  image.setValue(Z2i::Point(7,7), Vector(3,6)); 
//image.setValue(Z2i::Point(3,7), 0); 
  image.setValue(Z2i::Point(12,12), Vector(2,8));
 
  
Board2D board;
 
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
  board << domain << set;

  board.saveEPS("tolerancedBalls.eps");

  
  typedef PredicateCWSeparableMetric<Space, Z2i::Vector> CWMetric;
  CWMetric CWmetric;
  
  CWMap<Image,CWMetric> map(&domain);
  
  map.restrictedCWMapBruteForce(&domain,&image,&CWmetric);
  
  CWMap<Image, CWMetric> CWmap(&domain, &image, &CWmetric);


board.clear();

  for(CWMap<Image, CWMetric>::Domain::ConstIterator it = map.domain().begin(),
  	itend = map.domain().end(); it != itend; ++it)
    {
CWMap<Image,CWMetric>::Value site = map( *it );   //closest site to (*it)
float d=Z2i::l2Metric(site,*it);

      // if (site != (*it))
      // 	Display2DFactory::draw( board,   site - (*it), (*it)); //Draw an arrow
    }

  
  
  // Z2i::L2PowerMetric l2power;
  // PowerMap<Image, Z2i::L2PowerMetric> power(&domainLarge, &image, &l2power);
  // for(unsigned int i=0; i<11; i++)
  //   {
  //     for(unsigned int j=0; j<11; j++)
  // 	if (image.domain().isInside(Z2i::Point(i,j)))
  // 	  trace.info()<< image(Z2i::Point(i,j))<<" ";
  // 	else
  // 	  trace.info()<< "0 ";
  //     trace.info()<<std::endl;
  //   }



}
