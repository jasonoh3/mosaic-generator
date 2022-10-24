/**
 * @file maptiles.cpp
 * Code for the maptiles function.
 */

#include <iostream>
#include <map>
#include "maptiles.h"
//#include "cs225/RGB_HSL.h"

using namespace std;


Point<3> convertToXYZ(LUVAPixel pixel) {
    return Point<3>( pixel.l, pixel.u, pixel.v );
}

MosaicCanvas* mapTiles(SourceImage const& theSource,
                       vector<TileImage>& theTiles)
{
    MosaicCanvas* mosaic = new MosaicCanvas(theSource.getRows(), theSource.getColumns());
    vector<Point<3>> pointTiles;
    map<Point<3>, TileImage*> tileMap;

    for (unsigned i = 0; i < theTiles.size(); ++i) {
        pointTiles.push_back(convertToXYZ(theTiles[i].getAverageColor()));
        tileMap[pointTiles[i]] = &theTiles[i];
    }

    KDTree<3> tree(pointTiles);

    for (int x = 0; x < theSource.getRows(); ++x) {
        for (int y = 0; y < theSource.getColumns(); ++y) {
            Point<3> colorPoint = convertToXYZ(theSource.getRegionColor(x, y));
            Point<3> tile = tree.findNearestNeighbor(colorPoint);
            mosaic->setTile(x, y, tileMap[tile]);
        }
    }

    return mosaic;
}

