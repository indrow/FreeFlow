/*
 * Created by Indra on 8/7/2020.
 *
 * Handling user mesh and generate domain for computation.
 * Periodic mesh is supported.
 *
 */
#include "MeshHandler.h"

MeshHandler::MeshHandler() {
    cout << "Mesh handler class has been called.\n";

    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 1);

    numOfPeriodicBC = 0;

    cout << "Warning: Only structured mesh supported!\n";
}

void MeshHandler::readMesh(const std::string &filename, const int &stencilWidth, const int &numPeriodicBc,
                           meshDefinition &mesh) {
    auto start = std::chrono::high_resolution_clock::now();

    this->numOfPeriodicBC = numPeriodicBc;

    /* Checking file integrity */
    if (!isFileExist(filename)) {
        gmsh::logger::write("Your mesh file is not found!",
                            "error");
        exit(-1);
    } else {
        gmsh::open(filename);
    }

    /* Saving user input: stencil width */
    mesh.stencilWidth = stencilWidth;

    int cellType = unsupported()[0];
    int dimension, order, numNodes, numPrimaryNodes;

    std::size_t cells;
    string name;
    std::vector<double> parametricCoord;
    std::vector<std::size_t> cellTags, nodeTags, cellsNodes;

    /* Add information to user */
    gmsh::model::mesh::getElementProperties(cellType, name, dimension, order, numNodes,
                                            parametricCoord, numPrimaryNodes);
    gmsh::logger::write("" + std::to_string(dimension) + "D cells are of type '" +
                        name + "' (type = " + std::to_string(cellType) + ") ");

    gmsh::model::mesh::preallocateElementsByType(cellType, true, true, cellTags, nodeTags, -1);
    gmsh::model::mesh::getElementsByType(cellType, cellTags, nodeTags, -1);
    gmsh::model::mesh::getElementFaceNodes(cellType, numPrimaryNodes, cellsNodes, -1, true);

    mesh.dimension = dimension;

    /* Creating circular number for node/faceNumber */
    circularNumber.resize(numPrimaryNodes, 0);
    for (int i = 0; i < numPrimaryNodes; ++i) {
        circularNumber[i] = i;
    }

    std::cout << "Performing mesh operation...\t\t";

    mesh.cellFirstTag = cellTags.front();
    mesh.cellLastTag = cellTags.back();
    mesh.internalCells = mesh.cellLastTag - mesh.cellFirstTag;
    mesh.cell.resize(cellTags.back() + 1);
    cells = 0;

    /* Collecting cell node from gmsh */
    for (std::size_t i = 0; i < cellsNodes.size(); i += numPrimaryNodes) {
        mesh.cell[mesh.cellFirstTag + cells].node =
                std::vector<std::size_t>(cellsNodes.begin() + i, cellsNodes.begin() + i + numPrimaryNodes);
        cells += 1;
    }

    int numFaceNode, faceType = 1;
    gmsh::model::mesh::getElementProperties(faceType, name, dimension, order, numNodes, parametricCoord, numFaceNode);

    int numCellFace = numPrimaryNodes;
    mesh.numCellFace = numCellFace;
    std::vector<std::size_t> cellFaceNodes(numFaceNode), periodicFaceNodes(numFaceNode);
    map<std::vector<std::size_t>, std::vector<std::size_t>> cellsNearFace;
    map<std::vector<std::size_t>, std::size_t> faces;
    std::size_t totalFace = 1;

    /* Collect data for cells tag arround the same faceNumber */
    for (std::size_t j = mesh.cellFirstTag; j <= mesh.cellLastTag; ++j) {
        for (int f = 0; f < numCellFace; ++f) {

            for (int nn = 0; nn < numFaceNode; ++nn) {
                cellFaceNodes[nn] = mesh.cell[j].node[circularNumber[f + nn]];
            }

            std::sort(cellFaceNodes.begin(), cellFaceNodes.end());

            if (!(faces.count(cellFaceNodes))) {
                faces.insert(std::pair<std::vector<std::size_t>, std::size_t>(cellFaceNodes, totalFace));
                totalFace += 1;
            }

            if (!(cellsNearFace.count(cellFaceNodes))) {
                cellsNearFace.insert(
                        std::pair<std::vector<std::size_t>, std::vector<std::size_t>>(cellFaceNodes, {j}));
            } else {
                cellsNearFace[cellFaceNodes].push_back(j);
            }
        }
    }

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "elapsed time: " << duration.count() << "ms.\n";

    std::cout << "Mapping the boundary condition...\t";
    start = std::chrono::high_resolution_clock::now();

    mesh.faceFirstTag = 1;
    mesh.faceLastTag = totalFace;
    mesh.numFaces = mesh.faceLastTag - mesh.faceFirstTag;
    mesh.faces.resize(totalFace);

    getBoundary(mesh);
    std::vector<std::unordered_map<std::size_t, std::size_t>> periodicMap;
    periodicMap.resize(numPeriodicBc);

    /* If periodic boundary condition/s are defined,
     * create map of the related node */
    if (numOfPeriodicBC > 0) {
        int masterTags;
        std::vector<std::size_t> pNodeTags, pNodeMasterTags;
        std::vector<double> tfo;
        for (int i = 0; i < numPeriodicBc; ++i) {
            gmsh::model::mesh::getPeriodicNodes(1, i + 1, masterTags, pNodeTags, pNodeMasterTags, tfo, false);
            for (std::size_t j = 0; j < pNodeMasterTags.size(); ++j) {
                periodicMap[i].insert(std::pair<std::size_t, std::size_t>(pNodeMasterTags[j], pNodeTags[j]));
                periodicMap[i].insert(std::pair<std::size_t, std::size_t>(pNodeTags[j], pNodeMasterTags[j]));
            }
            std::vector<size_t>().swap(pNodeTags);
            std::vector<size_t>().swap(pNodeMasterTags);
        }
    }

    std::vector<std::size_t> boundaryFace(numFaceNode, 0), knownFace(numFaceNode, 0);
    map<std::vector<std::size_t>, std::size_t>::iterator faceIterator;

    /* Boundary tagging for faceNumber */
    for (faceIterator = faces.begin(); faceIterator != faces.end(); ++faceIterator) {
        mesh.faces[faceIterator->second].nodes = faceIterator->first;

        knownFace = mesh.faces[faceIterator->second].nodes;
        mesh.faces[faceIterator->second].boundaryTag = boundaryFaceMap[knownFace];
    }

    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "elapsed time: " << duration.count() << "ms.\n";

    std::vector<std::size_t> v(numPrimaryNodes, 0), w(numCellFace, 0);

    std::cout << "Collecting stencil...\t\t\t\t";
    start = std::chrono::high_resolution_clock::now();
    bool periodic_pair;

    /* Collect the neighboring cells */
    for (std::size_t c = mesh.cellFirstTag; c <= mesh.cellLastTag; ++c) {
        mesh.cell[c].neighbor = v;
        mesh.cell[c].faceNumber = w;

        for (int f = 0; f < numCellFace; ++f) {
            for (int nn = 0; nn < numFaceNode; ++nn) {
                cellFaceNodes[nn] = mesh.cell[c].node[circularNumber[f + nn]];
            }

            periodic_pair = false;
            for (int i = 0; i < numOfPeriodicBC; ++i) {
                periodic_pair = true;
                for (auto const &cellFNode : cellFaceNodes) {
                    periodic_pair = periodic_pair && periodicMap[i].count(cellFNode);
                }
                if (periodic_pair) {
                    for (int j = 0; j < numFaceNode; ++j) {
                        periodicFaceNodes[j] = periodicMap[i][cellFaceNodes[j]];
                    }
                    break;
                }
            }

            if (periodic_pair) {
                std::sort(periodicFaceNodes.begin(), periodicFaceNodes.end());
                for (std::size_t allElm : cellsNearFace[periodicFaceNodes]) {
                    if (allElm != c && allElm != 0)
                        mesh.cell[c].neighbor[f] = allElm;
                }
            }

            std::sort(cellFaceNodes.begin(), cellFaceNodes.end());
            mesh.cell[c].faceNumber[f] = faces[cellFaceNodes];

            for (std::size_t allElm : cellsNearFace[cellFaceNodes]) {
                if (allElm != c)
                    mesh.cell[c].neighbor[f] = allElm;
            }
        }
    }

    int facePair;
    std::string bcName;
    std::size_t recNeighbor;
    cellDefinition tmpCell;
    mesh.totalCells = mesh.cellLastTag;

    /* Create new ghost cell (if boundary is not periodic) */
    for (std::size_t j = mesh.cellFirstTag; j <= mesh.cellLastTag; ++j) {
        for (int f = 0; f < numCellFace; ++f) {
            if (mesh.cell[j].neighbor[f] == 0) {
                bcName = mesh.boundaryNames.right.at(mesh.faces[mesh.cell[j].faceNumber[f]].boundaryTag);
                recNeighbor = j;

                for (int c = 0; c < stencilWidth; ++c) {

                    mesh.cell.push_back(tmpCell);
                    mesh.totalCells += 1;

                    facePair = circularNumber[f + 2];

                    mesh.cell[mesh.totalCells].neighbor.resize(mesh.numCellFace);
                    mesh.cell[mesh.totalCells].bcTag = mesh.faces[mesh.cell[j].faceNumber[f]].boundaryTag;

                    if (c == 0) {
                        mesh.cell[mesh.totalCells].neighbor[facePair] = j;
                        mesh.cell[mesh.totalCells].neighbor[f] = mesh.totalCells + 1;
                        mesh.cell[j].neighbor[f] = mesh.totalCells;
                    } else if (c < stencilWidth - 1) {
                        mesh.cell[mesh.totalCells].neighbor[f] = mesh.totalCells + 1;
                        mesh.cell[mesh.totalCells].neighbor[facePair] = mesh.totalCells - 1;
                    } else {
                        mesh.cell[mesh.totalCells].neighbor[f] = 0;
                        mesh.cell[mesh.totalCells].neighbor[facePair] = mesh.totalCells - 1;
                    }
                    if (bcName == "Wall") {
                        mesh.cell[mesh.totalCells].bcCellRef = recNeighbor;
                        mesh.cell[mesh.totalCells].bcWallRef = j;
                    } else {
                        mesh.cell[mesh.totalCells].bcCellRef = j;
                    }
                    mesh.cell[mesh.totalCells].bcFaceRef = f;

                    recNeighbor = mesh.cell[recNeighbor].neighbor[circularNumber[f + 2]];
                }
            }
        }
    }

    collectStencil(stencilWidth, numCellFace, mesh);

    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "elapsed time: " << duration.count() << "ms.\n";

    std::cout << "Saving mesh geometry...\t\t\t\t";
    start = std::chrono::high_resolution_clock::now();

    cellGeometry(cellType, mesh);

    stop = std::chrono::high_resolution_clock::now();
    duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
    std::cout << "elapsed time: " << duration.count() << "ms.\n";
}

bool MeshHandler::isFileExist(const string &filename) {
    std::ifstream infile(filename);
    return infile.good();
}

std::vector<int> MeshHandler::unsupported() {
    int modelDimension;
    std::vector<int> cellTypes;

    modelDimension = gmsh::model::getDimension();
    gmsh::model::mesh::getElementTypes(cellTypes, modelDimension);

    if (cellTypes.size() != 1) {
        gmsh::logger::write("Hybrid meshes not handled in this code!",
                            "error");
        exit(1);
    } else if (cellTypes[0] != 3) {
        gmsh::logger::write("Support only for quadrangle cell type!",
                            "error");
        exit(1);
    } else if (modelDimension != 2) {
        gmsh::logger::write("Sorry, the current solver support 2D geometry only!",
                            "error");
        exit(1);
    } else {
        return cellTypes;
    }
}

void MeshHandler::getBoundary(meshDefinition &mesh) {
    int dimension, order, numNodes, numPrimaryNodes;
    gmsh::vectorpair dimTags;
    std::vector<int> cellTypes, entities;
    std::vector<std::size_t> cellTags, nodeTags, tmpNodes;
    std::vector<double> nodeCoord;
    string physicalGroupName, cellName;

    gmsh::model::getPhysicalGroups(dimTags, -1);

    /* Map the boundary names and the faceNumber related to it */
    for (auto &dimTag : dimTags) {
        if (dimTag.first == 1) {
            gmsh::model::getEntitiesForPhysicalGroup(dimTag.first, dimTag.second, entities);

            for (auto &entitie : entities) {
                gmsh::model::mesh::getElementTypes(cellTypes, dimTag.first, entitie);

                for (auto &elmType : cellTypes) {
                    gmsh::model::mesh::preallocateElementsByType(elmType, true, true, cellTags, nodeTags, entitie);
                    gmsh::model::mesh::getElementsByType(elmType, cellTags, nodeTags, entitie, 0, 1);
                    gmsh::model::mesh::getElementProperties(elmType, cellName, dimension, order, numNodes,
                                                            nodeCoord, numPrimaryNodes);
                    gmsh::model::getPhysicalName(dimTag.first, dimTag.second, physicalGroupName);

                    if (!mesh.boundaryNames.right.count(dimTag.second)) {
                        mesh.boundaryNames.insert(bimap<string, int>::value_type(physicalGroupName, dimTag.second));
                    }

                    for (std::size_t i = 0; i < nodeTags.size(); i += numPrimaryNodes) {
                        tmpNodes = std::vector<std::size_t>
                                (nodeTags.begin() + i, nodeTags.begin() + i + numPrimaryNodes);
                        std::sort(tmpNodes.begin(), tmpNodes.end());

                        boundaryFaceMap.insert(std::pair<std::vector<std::size_t>, int>(tmpNodes, dimTag.second));
                    }
                }
            }
        }
    }

    if (mesh.boundaryNames.empty()) {
        std::cout << "Please set the boundary conditions in the mesh generation.\n";
        exit(1);
    }
}

void MeshHandler::collectStencil(const int &stencilWidth, const int &numCellFace, meshDefinition &mesh) {
    int k = stencilWidth, subStencilWidth = (k - 1) / 2;

    /* Collect the stencil from known neighboring data. */
    for (std::size_t i = mesh.cellFirstTag; i <= mesh.totalCells; ++i) {
        mesh.cell[i].stencil.resize(k * 2);
        mesh.cell[i].stencil[subStencilWidth] = i;
        mesh.cell[i].stencil[subStencilWidth + k] = i;

        for (int f = 0; f < numCellFace; ++f) {
            for (int j = 1; j < subStencilWidth + 1; ++j) {
                mesh.cell[i].stencil[subStencilWidth - j] =
                        mesh.cell[i].stencil[subStencilWidth - j + 1] == 0 ?
                        0 : mesh.cell[mesh.cell[i].stencil[subStencilWidth - j + 1]].neighbor[west];
                mesh.cell[i].stencil[subStencilWidth + j] =
                        mesh.cell[i].stencil[subStencilWidth + j - 1] == 0 ?
                        0 : mesh.cell[mesh.cell[i].stencil[subStencilWidth + j - 1]].neighbor[east];
                mesh.cell[i].stencil[subStencilWidth - j + k] =
                        mesh.cell[i].stencil[subStencilWidth - j + 1 + k] == 0 ?
                        0 : mesh.cell[mesh.cell[i].stencil[subStencilWidth - j + 1 + k]].neighbor[south];
                mesh.cell[i].stencil[subStencilWidth + j + k] =
                        mesh.cell[i].stencil[subStencilWidth + j - 1 + k] == 0 ?
                        0 : mesh.cell[mesh.cell[i].stencil[subStencilWidth + j - 1 + k]].neighbor[north];
            }
        }
    }
}

void MeshHandler::cellGeometry(const int &cellType, meshDefinition &mesh) {
    std::vector<std::size_t> nodeTags;
    std::vector<double> parametricCoord, nodeCoord, coord;
    gmsh::model::mesh::getNodes(nodeTags, nodeCoord, parametricCoord, -1, -1);

    mesh.nodeFirstTag = nodeTags.front();
    mesh.nodeLastTag = nodeTags.back();
    mesh.numNodes = mesh.nodeLastTag - mesh.nodeFirstTag;

    mesh.nodes.resize(mesh.nodeLastTag + 1);

    /* Save the node coordinate data */
    for (std::size_t i = 0; i <= mesh.numNodes; ++i) {
        coord = std::vector<double>(nodeCoord.begin() + i * 3, nodeCoord.begin() + (i + 1) * 3);
        mesh.nodes[i + mesh.nodeFirstTag].coordinate = Eigen::Map<Vector3d>(coord.data());
    }

    std::vector<double> barycenter, barycenterMap(3, 0.0);
    gmsh::model::mesh::preallocateBarycenters(cellType, barycenter, -1);
    gmsh::model::mesh::getBarycenters(cellType, -1, false, true, barycenter, 0, 1);

    /* Save the barycenter, transfMetrics of transformation and its determinant to the
     * mesh data. */
    std::size_t ii = 0;
    for (std::size_t i = mesh.cellFirstTag; i <= mesh.cellLastTag; ++i) {
        barycenterMap = std::vector<double>(barycenter.begin() + ii * 3, barycenter.begin() + (ii + 1) * 3);
        mesh.cell[i].barycenter = Eigen::Map<Vector3d>(barycenterMap.data());
        ii += 1;
    }

    std::size_t neighbor;
    /* Calculate the difference of each cell faceNumber. */
    mesh.minGridSize = 1.0E300;
    for (std::size_t i = mesh.cellFirstTag; i <= mesh.cellLastTag; ++i) {
        mesh.cell[i].face.resize(mesh.numCellFace);
        mesh.cell[i].area = 0.0;

        for (int f = 0; f < mesh.numCellFace; ++f) {
            mesh.cell[i].face[f].d = mesh.nodes[mesh.cell[i].node[circularNumber[f + 1]]].coordinate.segment<2>(0) -
                                     mesh.nodes[mesh.cell[i].node[circularNumber[f]]].coordinate.segment<2>(0);
            mesh.cell[i].face[f].ds = mesh.cell[i].face[f].d.norm();
            mesh.cell[i].face[f].norm(x) = mesh.cell[i].face[f].d(y) / mesh.cell[i].face[f].ds;
            mesh.cell[i].face[f].norm(y) = -mesh.cell[i].face[f].d(x) / mesh.cell[i].face[f].ds;
            mesh.minGridSize = fmin(mesh.minGridSize, mesh.cell[i].face[f].ds);
            mesh.cell[i].area += mesh.nodes[mesh.cell[i].node[circularNumber[f]]].coordinate(x) *
                                 mesh.nodes[mesh.cell[i].node[circularNumber[f + 1]]].coordinate(y) -
                                 mesh.nodes[mesh.cell[i].node[circularNumber[f + 1]]].coordinate(x) *
                                 mesh.nodes[mesh.cell[i].node[circularNumber[f]]].coordinate(y);

            neighbor = mesh.cell[i].neighbor[f] <= mesh.cellLastTag ? mesh.cell[i].neighbor[f] :
                       mesh.cell[i].neighbor[circularNumber[f + 2]];

            mesh.cell[i].face[f].dn = mesh.cell[neighbor].barycenter.segment<2>(0) -
                                      mesh.cell[neighbor].barycenter.segment<2>(0);
            mesh.cell[i].face[f].dsn = mesh.cell[i].face[f].dn.norm();
//            std::cout << mesh.nodes[mesh.cell[i].node[circularNumber[f + 1]]].coordinate.transpose() << "\n======\n";

        }
//        std::cout << i << " : ";
//        for (auto const &ng : mesh.cell[i].node) {
//            std::cout << mesh.nodes[ng].coordinate.transpose() << "\t";
//        }
//        std::cout << "\n";
        mesh.cell[i].area = 0.5 * fabs(mesh.cell[i].area);
    }

    std::vector<double> localCoord(12, 0.0), jacobians, determinants, points, midpoint(3, 0.0);
    localCoord[1] = -1.0;
    localCoord[3] = 1.0;
    localCoord[7] = 1.0;
    localCoord[9] = -1.0;

    gmsh::model::mesh::preallocateJacobians(cellType, 4, true, true, true, jacobians, determinants, points, -1);
    gmsh::model::mesh::getJacobians(cellType, localCoord, jacobians, determinants, points, -1, 0, 1);

    ii = 0;
    for (std::size_t i = mesh.cellFirstTag; i <= mesh.cellLastTag; ++i) {
        for (int f = 0; f < mesh.numCellFace; ++f) {
            midpoint = std::vector<double>(points.begin() + ii * 3, points.begin() + (ii + 1) * 3);
            mesh.cell[i].face[f].midpoint = Eigen::Map<Vector3d>(midpoint.data()).segment(0, 2);
            ii += 1;
        }
    }

    for (std::size_t i = mesh.cellFirstTag; i <= mesh.cellLastTag; ++i) {
        mesh.cell[i].h(x) = sqrt(pow(mesh.cell[i].face[1].midpoint(x) - mesh.cell[i].face[3].midpoint(x), 2) +
                                 pow(mesh.cell[i].face[1].midpoint(y) - mesh.cell[i].face[3].midpoint(y), 2));
        mesh.cell[i].h(y) = sqrt(pow(mesh.cell[i].face[2].midpoint(x) - mesh.cell[i].face[0].midpoint(x), 2) +
                                 pow(mesh.cell[i].face[2].midpoint(y) - mesh.cell[i].face[0].midpoint(y), 2));
    }
}

MeshHandler::~MeshHandler() {
    std::cout << "Destructor of mesh handler has been called.\n";
    gmsh::finalize();
}
