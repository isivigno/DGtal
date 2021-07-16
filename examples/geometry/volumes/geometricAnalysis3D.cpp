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
 * @file geometry/volumes/geometricAnalysis3D.cpp
 * @ingroup Examples
 * @author Jacques-Olivier Lachaud (\c jacques-olivier.lachaud@univ-savoie.fr )
 * Laboratory of Mathematics (CNRS, UMR 5127), University of Savoie, France
 *
 * @date 2021/06/20
 *
 * An example file named geometricAnalysis3D
 *
 * This file is part of the DGtal library.
 */


/**
 * This example shows how to analyze the local geometry of 3D digital
 * sets with full convexity over cubical neighborhoods.
 *
 * @see \ref dgtal_topology_sec3_5
 *
 *
 * \example geometry/volumes/geometricAnalysis3D.cpp
 */


///////////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <queue>
#include "DGtal/base/Common.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/io/DrawWithDisplay3DModifier.h"
#include "DGtal/io/Color.h"
#include "DGtal/shapes/Shapes.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/helpers/Shortcuts.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
#include "DGtal/geometry/volumes/NeighborhoodConvexityAnalyzer.h"

///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
using namespace Z3i;
typedef Shortcuts< KSpace > SH3;

///////////////////////////////////////////////////////////////////////////////


template < typename KSpace, int N >
struct Analyzer {
  typedef typename KSpace::Point Point;
  typedef NeighborhoodConvexityAnalyzer< KSpace, N > NCA;

  template < typename ImagePtr >
  static
  std::vector<int>
  run( const KSpace& aK, std::vector<Point> pts, ImagePtr bimage )
  {
    NCA nca( aK.lowerBound(), aK.upperBound(),
             KSpace::dimension <= 2 ? 0 : 10000*KSpace::dimension*N );
    auto& image = *bimage;
    std::vector<int> result;
    std::map< Point, int > computed;
    int geom;
    int i  = 0;
    int nb = pts.size();
    for ( auto p : pts )
      {
        if ( i % 100 == 0 ) trace.progressBar( i, nb );
        auto it = computed.find( p );
        if ( it == computed.end() )
          {
            nca.setCenter( p, image );
            bool  cvx = nca.isFullyConvex( true );
            bool ccvx = nca.isComplementaryFullyConvex( false );
            geom = ( cvx ? 0x1 : 0x0 ) | ( ccvx ? 0x2 : 0x0 );
            computed[ p ] = geom;
          }
        else geom = it->second;
        result.push_back( geom );
        i++; 
      }
    return result;    
  }

  template < typename ImagePtr >
  static
  void
  run( std::vector<int> & to_update,
       const KSpace& aK, std::vector<Point> pts, ImagePtr bimage )
  {
    NCA nca( aK.lowerBound(), aK.upperBound(), 
             KSpace::dimension <= 2 ? 0 : 10000*KSpace::dimension*N );
    auto& image = *bimage;
    std::map< Point, int > computed;
    int geom;
    int i  = 0;
    int nb = pts.size();
    for ( auto p : pts )
      {
        if ( i % 100 == 0 ) trace.progressBar( i, nb );
        auto it = computed.find( p );
        if ( it == computed.end() )
          {
            nca.setCenter( p, image );
            bool  cvx = ( to_update[ i ] & 0x1 )
              ? nca.isFullyConvex( true )
              : false;
            bool ccvx = ( to_update[ i ] & 0x2 )
              ? nca.isComplementaryFullyConvex( false )
              : false;
            geom = ( cvx ? 0x1 : 0x0 ) | ( ccvx ? 0x2 : 0x0 );
            computed[ p ] = geom;
          }
        else geom = it->second;
        to_update[ i++ ] = geom;
      }
  }
};

template < typename KSpace, int N >
struct MultiScaleAnalyzer {
  typedef typename KSpace::Point Point;
  typedef NeighborhoodConvexityAnalyzer< KSpace, N > NCA;
  typedef std::pair< int, int > Geometry; // convex, concave

  template < typename ImagePtr >
  static
  std::vector< Geometry >
  multiscale_run( const KSpace& aK,
                  std::vector<Point> pts,
                  ImagePtr bimage )
  {
    auto prev_geometry
      = MultiScaleAnalyzer< KSpace, N-1>::multiscale_run( aK, pts, bimage );
    trace.info() << "------- Analyzing scale " << N << " --------" << std::endl;
    std::vector< int > geom( prev_geometry.size() );
    for ( int i = 0; i < geom.size(); i++ )
      geom[ i ] = ( prev_geometry[ i ].first == N-1 ? 0x1 : 0x0 )
        |  ( prev_geometry[ i ].second == N-1 ? 0x2 : 0x0 );
    Analyzer< KSpace, N>::run( geom, aK, pts, bimage );
    for ( int i = 0; i < geom.size(); i++ ) {
      prev_geometry[ i ].first  += ( geom[ i ] & 0x1 ) ? 1 : 0;
      prev_geometry[ i ].second += ( geom[ i ] & 0x2 ) ? 1 : 0;
    }
    return prev_geometry;
  }
};

/// Specialization
template < typename KSpace>
struct MultiScaleAnalyzer< KSpace, 0 > {
  typedef typename KSpace::Point Point;
  typedef std::pair< int, int > Geometry; // convex, concave

  template < typename ImagePtr >
  static
  std::vector< Geometry >
  multiscale_run( const KSpace& aK,
                  std::vector<Point> pts,
                  ImagePtr bimage )
  {
    return std::vector< Geometry >( pts.size(), std::make_pair( 0, 0 ) );
  }
};

int main( int argc, char** argv )
{
  if ( argc <= 2 )
    {
      trace.info() << "Usage: " << argv[ 0 ] << " <K> <input.vol> <m> <M>" << std::endl;
      trace.info() << "\tAnalyze the  shape with local full convexity" << std::endl;
      trace.info() << "\t- 1 <= K <= 5: analysis at scale K" << std::endl;
      trace.info() << "\t- K == 0: multiscale analysis (using scales 1-5)" << std::endl;
      trace.info() << "\t- input.vol: choose your favorite shape" << std::endl;
      trace.info() << "\t- m [==0], M [==255]: used to threshold input vol image" << std::endl;
      return 1;
    }
  int         N = argc > 1 ? atoi( argv[ 1 ] ) : 1;
  std::string fn= argc > 2 ? argv[ 2 ]         : "";
  int         m = argc > 3 ? atoi( argv[ 3 ] ) : 0;
  int         M = argc > 4 ? atoi( argv[ 4 ] ) : 255;

  QApplication application(argc,argv);
  Viewer3D<> viewer;
  viewer.setWindowTitle("simpleExample3DViewer");
  viewer.show();  

  auto   params  = SH3::defaultParameters();
  
  // Domain creation from two bounding points.
  trace.info() << "Building set or importing vol ... ";
  KSpace K;
  params( "thresholdMin", m );
  params( "thresholdMax", M );
  auto bimage = SH3::makeBinaryImage( fn, params );
  K = SH3::getKSpace( bimage );
  Point p1 = K.lowerBound();
  Point p2 = K.upperBound();
  Domain domain = Domain( p1, p2 );
  trace.info() << "  [Done]" << std::endl;
  // Compute surface
  params( "surfaceComponents" , "All" );
  auto surface = SH3::makeDigitalSurface( bimage, K, params );
  // Compute interior boundary points
  std::vector< Point > points;
  for ( auto s : (*surface) )
    {
      // get inside point on the border of the shape.
      Dimension k = K.sOrthDir( s );
      auto voxel  = K.sIncident( s, k, K.sDirect( s, k ) );
      Point p     = K.sCoords( voxel );
      points.push_back( p );
    } 
  trace.info() << "Shape has " << point.size() << " interior boundary points"
               << std::endl;
  if ( N != 0 )
    {
      std::vector< int > result;
      trace.beginBlock ( "Analyzing" );
      if ( N == 1 ) result = Analyzer< KSpace, 1 >::run( K, points, bimage );
      if ( N == 2 ) result = Analyzer< KSpace, 2 >::run( K, points, bimage );
      if ( N == 3 ) result = Analyzer< KSpace, 3 >::run( K, points, bimage );
      if ( N == 4 ) result = Analyzer< KSpace, 4 >::run( K, points, bimage );
      if ( N == 5 ) result = Analyzer< KSpace, 5 >::run( K, points, bimage );
      trace.endBlock();
      SCell dummy;
      Color colors[ 4 ] =
        { Color( 255, 0, 0, 255 ), Color( 0, 255, 0, 255 ),
          Color( 0, 0, 255, 255 ), Color( 255, 255, 255, 255 ) };
      auto surfels   = SH3::getSurfelRange( surface, params );
      SH3::Colors all_colors( surfels.size() );
      for ( int i = 0; i < surfels.size(); i++ )
        all_colors[ i ] = colors[ result[ i ] ];
      bool ok = SH3::saveOBJ( surface, SH3::RealVectors(), all_colors,
                              "geom-cvx.obj" );
      int i = 0;
      viewer << SetMode3D( dummy.className(), "Basic" );
      for ( auto s : (*surface) )
        {
          viewer << CustomColors3D( all_colors[ i ], all_colors[ i ] )
                 << s;
          i++;
        }
    }
  else
    {
      trace.beginBlock ( "Analyzing" );
      auto geometry =
        MultiScaleAnalyzer< KSpace, 5 >::multiscale_run( K, points, bimage );
      trace.endBlock();
      Color colors_planar[ 6 ] =
        { Color( 0, 255, 255, 255),
          Color( 50, 255, 255, 255), Color( 100, 255, 255, 255),
          Color( 150, 255, 255, 255), Color( 200, 255, 255, 255 ),
          Color( 255, 255, 255, 255 ) };
      Color color_atypical( 255, 0, 0, 255 ); 
      Color colors_cvx[ 5 ] =
        { Color( 0, 255, 0, 255 ), Color( 50, 255, 50, 255 ),
          Color( 100, 255, 100, 255 ), Color( 150, 255, 150, 255 ),
          Color( 200, 255, 200, 255 ) };
      Color colors_ccv[ 5 ] =
        { Color( 0, 0, 255, 255 ), Color( 50, 50, 255, 255 ),
          Color( 100, 100, 255, 255 ), Color( 150, 150, 255, 255 ),
          Color( 200, 200, 255, 255 ) };
      auto surfels   = SH3::getSurfelRange( surface, params );
      SH3::Colors all_colors( surfels.size() );
      for ( int i = 0; i < surfels.size(); i++ ) {
        int m0 = std::min( geometry[ i ].first, geometry[ i ].second );
        int m1 = std::max( geometry[ i ].first, geometry[ i ].second );
        if ( m1 == 0 ) all_colors[ i ] = color_atypical;
        else if ( m0 == m1 ) all_colors[ i ] = colors_planar[ 5 ];
        else if ( geometry[ i ].first > geometry[ i ].second )
          all_colors[ i ] = colors_cvx[ 5 - abs( m0 - m1 ) ];
        else
          all_colors[ i ] = colors_ccv[ 5 - abs( m0 - m1 ) ];
      }
      bool ok = SH3::saveOBJ( surface, SH3::RealVectors(), all_colors,
                              "geom-scale-cvx.obj" );
      SCell dummy;
      int i = 0;
      viewer << SetMode3D( dummy.className(), "Basic" );
      for ( auto s : (*surface) )
        {
          viewer << CustomColors3D( all_colors[ i ], all_colors[ i ] )
                 << s;
          i++;
        }
    }      
  viewer<< Viewer3D<>::updateDisplay;
  return application.exec();
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

