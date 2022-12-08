//
// Created by Indra on 8/7/2020.
//

#ifndef FREEFLOW_MESHHANDLER_H
#define FREEFLOW_MESHHANDLER_H

#include "../include/CoordinateSystems.h"
#include "../include/EigenType.h"
#include <boost/bimap.hpp>
#include <boost/circular_buffer.hpp>
#include <chrono>
#include <fstream>
#include <gmsh.h>
#include <iostream>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

using boost::bimap;
using std::cout;
using std::map;
using std::string;

class MeshHandler {
public:
    struct nodeDefinition {
        Vector3d coordinate;
    };

    struct faceDefinition {
        std::vector<std::size_t> nodes;
        int boundaryTag = 0;
    };

    struct faceProperties {
        Vector2d d, norm, dn, midpoint;
        double ds, dsn;
    };

    struct cellDefinition {
        double area;
        Vector3d barycenter;
        int bcTag = 0, bcFaceRef;
        std::size_t bcCellRef, bcWallRef;

        Matrix<std::size_t, Dynamic, 1> stencil;
        std::vector<std::size_t> node, neighbor, faceNumber;

        std::vector<faceProperties> face;
        Vector2d h;
    };

    struct meshDefinition {
        std::size_t numNodes, nodeFirstTag, nodeLastTag;
        std::size_t faceFirstTag, faceLastTag;
        std::size_t numFaces;
        std::size_t internalCells, cellFirstTag, cellLastTag, totalCells;
        std::vector<nodeDefinition> nodes;
        std::vector<faceDefinition> faces;
        std::vector<cellDefinition> cell;
        bimap<string, int> boundaryNames;
        int stencilWidth, numCellFace, dimension;
        double minGridSize;
    };

    MeshHandler();

    void readMesh(const string &filename, const int &stencilWidth, const int &numPeriodicBc, meshDefinition &mesh);

    ~MeshHandler();

private:
    static bool isFileExist(const string &filename);
    static std::vector<int> unsupported();
    void getBoundary(meshDefinition &mesh);
    static void collectStencil(const int &stencilWidth, const int &numCellFace, meshDefinition &mesh);
    void cellGeometry(const int &cellType, meshDefinition &mesh);

    int numOfPeriodicBC;
    boost::circular_buffer<int> circularNumber;
    std::unordered_map<std::vector<std::size_t>, int, hasher> boundaryFaceMap;
};


#endif //FREEFLOW_MESHHANDLER_H
