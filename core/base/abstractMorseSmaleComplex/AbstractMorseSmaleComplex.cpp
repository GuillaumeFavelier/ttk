#include<AbstractMorseSmaleComplex.h>

AbstractMorseSmaleComplex::AbstractMorseSmaleComplex():
  ReverveSaddleMaximumConnection{true},
  ReverveSaddleSaddleConnection{true},
  ComputeAscendingSeparatrices1{true},
  ComputeDescendingSeparatrices1{true},
  ComputeSaddleConnectors{true},
  ComputeAscendingSeparatrices2{false},
  ComputeDescendingSeparatrices2{false},
  ComputeAscendingSegmentation{true},
  ComputeDescendingSegmentation{true},
  ComputeFinalSegmentation{true},

  inputScalarField_{},
  inputTriangulation_{},
  inputOffsets_{},

  outputCriticalPoints_numberOfPoints_{},
  outputCriticalPoints_points_{},
  outputCriticalPoints_points_cellDimensions_{},
  outputCriticalPoints_points_cellIds_{},
  outputCriticalPoints_points_cellScalars_{},
  outputCriticalPoints_points_isOnBoundary_{},
  outputCriticalPoints_points_PLVertexIdentifiers_{},
  outputCriticalPoints_points_manifoldSize_{},

  outputSeparatrices1_numberOfPoints_{},
  outputSeparatrices1_points_{},
  outputSeparatrices1_points_smoothingMask_{},
  outputSeparatrices1_points_cellDimensions_{},
  outputSeparatrices1_points_cellIds_{},
  outputSeparatrices1_numberOfCells_{},
  outputSeparatrices1_cells_{},
  outputSeparatrices1_cells_sourceIds_{},
  outputSeparatrices1_cells_destinationIds_{},
  outputSeparatrices1_cells_separatrixIds_{},
  outputSeparatrices1_cells_separatrixTypes_{},
  outputSeparatrices1_cells_separatrixFunctionMaxima_{},
  outputSeparatrices1_cells_separatrixFunctionMinima_{},
  outputSeparatrices1_cells_separatrixFunctionDiffs_{},
  outputSeparatrices1_cells_isOnBoundary_{},

  outputSeparatrices2_numberOfPoints_{},
  outputSeparatrices2_points_{},
  outputSeparatrices2_numberOfCells_{},
  outputSeparatrices2_cells_{},
  outputSeparatrices2_cells_sourceIds_{},
  outputSeparatrices2_cells_separatrixIds_{},
  outputSeparatrices2_cells_separatrixTypes_{},
  outputSeparatrices2_cells_separatrixFunctionMaxima_{},
  outputSeparatrices2_cells_separatrixFunctionMinima_{},
  outputSeparatrices2_cells_separatrixFunctionDiffs_{},
  outputSeparatrices2_cells_isOnBoundary_{},

  outputAscendingManifold_{},
  outputDescendingManifold_{},
  outputMorseSmaleManifold_{}
{
  discreteGradient_.setReverseSaddleMaximumConnection(true);
  discreteGradient_.setReverseSaddleSaddleConnection(true);
}

AbstractMorseSmaleComplex::~AbstractMorseSmaleComplex(){
}

int AbstractMorseSmaleComplex::getDescendingSeparatrices1(const vector<Cell>& criticalPoints,
    vector<Separatrix>& separatrices,
    vector<vector<Cell>>& separatricesGeometry) const{

  vector<int> saddleIndexes;
  const int numberOfCriticalPoints=criticalPoints.size();
  for(int i=0; i<numberOfCriticalPoints; ++i){
    const Cell& criticalPoint=criticalPoints[i];

    if(criticalPoint.dim_==1)
      saddleIndexes.push_back(i);
  }
  const int numberOfSaddles=saddleIndexes.size();

  // estimation of the number of separatrices, apriori : numberOfAscendingPaths=2, numberOfDescendingPaths=2
  const int numberOfSeparatrices=4*numberOfSaddles;
  separatrices.resize(numberOfSeparatrices);
  separatricesGeometry.resize(numberOfSeparatrices);

  // apriori: by default construction, the separatrices are not valid
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i=0; i<numberOfSaddles; ++i){
    const int saddleIndex=saddleIndexes[i];
    const Cell& saddle=criticalPoints[saddleIndex];

    // add descending vpaths
    {
      const Cell& saddle1=saddle;

      for(int j=0; j<2; ++j){
        const int shift=j+2;

        int vertexId;
        inputTriangulation_->getEdgeVertex(saddle1.id_, j, vertexId);

        vector<Cell> vpath;
        vpath.push_back(saddle1);
        discreteGradient_.getDescendingPath(Cell(0,vertexId), vpath);

        const Cell& lastCell=vpath.back();
        if(lastCell.dim_==0 and discreteGradient_.isCellCritical(lastCell)){
          const int separatrixIndex=4*i+shift;

          separatricesGeometry[separatrixIndex]=std::move(vpath);
          separatrices[separatrixIndex]=std::move(Separatrix(true,saddle,lastCell,false,separatrixIndex));
        }
      }
    }
  }

  return 0;
}

int AbstractMorseSmaleComplex::setAscendingSegmentation(const vector<Cell>& criticalPoints,
    vector<int>& maxSeeds,
    int* const morseSmaleManifold,
    int& numberOfMaxima) const{
  const int numberOfVertices=inputTriangulation_->getNumberOfVertices();
  std::fill(morseSmaleManifold,morseSmaleManifold+numberOfVertices, -1);

  const int numberOfCells=inputTriangulation_->getNumberOfCells();
  vector<int> morseSmaleManifoldOnCells(numberOfCells, -1);
  const int cellDim=inputTriangulation_->getDimensionality();

  // get the seeds : maxima
  const int numberOfCriticalPoints=criticalPoints.size();
  for(int i=0; i<numberOfCriticalPoints; ++i){
    const Cell& criticalPoint=criticalPoints[i];

    if(criticalPoint.dim_==cellDim)
      maxSeeds.push_back(criticalPoint.id_);
  }
  const int numberOfSeeds=maxSeeds.size();
  numberOfMaxima=numberOfSeeds;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i=0; i<numberOfSeeds; ++i){
    queue<int> bfs;

    // push the seed
    {
      const int seedId=maxSeeds[i];
      bfs.push(seedId);
    }

    // BFS traversal
    while(!bfs.empty()){
      const int cofacetId=bfs.front();
      bfs.pop();

      if(morseSmaleManifoldOnCells[cofacetId]==-1){
        morseSmaleManifoldOnCells[cofacetId]=i;

        for(int j=0; j<(cellDim+1); ++j){
          int facetId=-1;
          if(cellDim==2)
            inputTriangulation_->getCellEdge(cofacetId, j, facetId);
          else if(cellDim==3)
            inputTriangulation_->getCellTriangle(cofacetId, j, facetId);

          int starNumber=0;
          if(cellDim==2)
            starNumber=inputTriangulation_->getEdgeStarNumber(facetId);
          else if(cellDim==3)
            starNumber=inputTriangulation_->getTriangleStarNumber(facetId);
          for(int k=0; k<starNumber; ++k){
            int neighborId=-1;
            if(cellDim==2)
              inputTriangulation_->getEdgeStar(facetId, k, neighborId);
            else if(cellDim==3)
              inputTriangulation_->getTriangleStar(facetId, k, neighborId);

            const int pairedCellId=discreteGradient_.getPairedCell(Cell(cellDim, neighborId), true);

            if(pairedCellId==facetId)
              bfs.push(neighborId);
          }
        }
      }
    }
  }

  // put segmentation infos from cells to points
#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i=0; i<numberOfVertices; ++i){
    int starId;
    inputTriangulation_->getVertexStar(i, 0, starId);
    morseSmaleManifold[i]=morseSmaleManifoldOnCells[starId];
  }

  return 0;
}

int AbstractMorseSmaleComplex::setDescendingSegmentation(const vector<Cell>& criticalPoints,
    int* const morseSmaleManifold,
    int& numberOfMinima) const{
  const int numberOfVertices=inputTriangulation_->getNumberOfVertices();
  std::fill(morseSmaleManifold,morseSmaleManifold+numberOfVertices, -1);

  // get the seeds : minima
  vector<int> seeds;
  const int numberOfCriticalPoints=criticalPoints.size();
  for(int i=0; i<numberOfCriticalPoints; ++i){
    const Cell& criticalPoint=criticalPoints[i];

    if(criticalPoint.dim_==0)
      seeds.push_back(criticalPoint.id_);
  }
  const int numberOfSeeds=seeds.size();
  numberOfMinima=numberOfSeeds;

#ifdef TTK_ENABLE_OPENMP
#pragma omp parallel for num_threads(threadNumber_)
#endif
  for(int i=0; i<numberOfSeeds; ++i){
    queue<int> bfs;

    // push the seed
    {
      const int seedId=seeds[i];
      bfs.push(seedId);
    }

    // BFS traversal
    while(!bfs.empty()){
      const int vertexId=bfs.front();
      bfs.pop();

      if(morseSmaleManifold[vertexId]==-1){
        morseSmaleManifold[vertexId]=i;

        const int edgeNumber=inputTriangulation_->getVertexEdgeNumber(vertexId);
        for(int j=0; j<edgeNumber; ++j){
          int edgeId;
          inputTriangulation_->getVertexEdge(vertexId, j, edgeId);

          for(int k=0; k<2; ++k){
            int neighborId;
            inputTriangulation_->getEdgeVertex(edgeId, k, neighborId);

            const int pairedCellId=discreteGradient_.getPairedCell(Cell(0, neighborId));

            if(pairedCellId==edgeId)
              bfs.push(neighborId);
          }
        }
      }
    }
  }

  return 0;
}

int AbstractMorseSmaleComplex::setFinalSegmentation(const int numberOfMaxima,
    const int numberOfMinima,
    const int* const ascendingManifold,
    const int* const descendingManifold,
    int* const morseSmaleManifold) const{
  vector<vector<pair<int,int>>> minTable(numberOfMinima);

  int id{};
  const int numberOfVertices=inputTriangulation_->getNumberOfVertices();
  for(int i=0; i<numberOfVertices; ++i){
    const int d=ascendingManifold[i];
    const int a=descendingManifold[i];

    if(a==-1 or d==-1){
      morseSmaleManifold[i]=-1;
      continue;
    }

    vector<pair<int,int>>& table=minTable[a];
    int foundId=-1;
    for(const pair<int,int>& p : table){
      if(p.first == d)
        foundId=p.second;
    }

    // add new association (a,d)
    if(foundId==-1){
      table.push_back(make_pair(d,id));
      morseSmaleManifold[i]=id;
      ++id;
    }
    // update to saved associationId
    else
      morseSmaleManifold[i]=foundId;
  }

  return 0;
}
