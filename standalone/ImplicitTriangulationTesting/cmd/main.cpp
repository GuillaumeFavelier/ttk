// include the local headers
#include                  <vtkProgramBase.h>

vtkUnstructuredGrid* build3dTriangulatedUnstructuredGrid(double origin[3],
    double spacing[3],
    int dimensions[3]){
  vtkUnstructuredGrid* ug=vtkUnstructuredGrid::New();
  vtkSmartPointer<vtkPoints> points=vtkSmartPointer<vtkPoints>::New();

  double p[3];
  for(int k=0; k<dimensions[2]; ++k){
    for(int j=0; j<dimensions[1]; ++j){
      for(int i=0; i<dimensions[0]; ++i){
        p[0]=origin[0]+spacing[0]*i;
        p[1]=origin[1]+spacing[1]*j;
        p[2]=origin[2]+spacing[2]*k;
        points->InsertNextPoint(p);
      }
    }
  }
  ug->SetPoints(points);

  ug->Allocate();

  int n=0;
  vtkIdType ids[4];
  int shiftJ=dimensions[0];
  int shiftK=dimensions[0]*dimensions[1];
  for(int k=0; k<dimensions[2]-1; ++k){
    for(int j=0; j<dimensions[1]-1; ++j){
      for(int i=0; i<dimensions[0]-1; ++i){

        n=i+j*shiftJ+k*shiftK;

        //T1
        ids[0]=n;
        ids[1]=n+1;
        ids[2]=n+shiftJ;
        ids[3]=n+shiftK+shiftJ;
        ug->InsertNextCell(VTK_TETRA,4,ids);

        //T2
        ids[0]=n+1;
        ids[1]=n+shiftJ;
        ids[2]=n+shiftJ+1;
        ids[3]=n+shiftK+shiftJ;
        ug->InsertNextCell(VTK_TETRA,4,ids);

        //T3
        ids[0]=n;
        ids[1]=n+1;
        ids[2]=n+shiftK;
        ids[3]=n+shiftK+shiftJ;
        ug->InsertNextCell(VTK_TETRA,4,ids);

        //T4
        ids[0]=n+1;
        ids[1]=n+shiftK;
        ids[2]=n+shiftK+1;
        ids[3]=n+shiftK+shiftJ;
        ug->InsertNextCell(VTK_TETRA,4,ids);

        //T5
        ids[0]=n+1;
        ids[1]=n+shiftK+1;
        ids[2]=n+shiftK+shiftJ;
        ids[3]=n+shiftK+shiftJ+1;
        ug->InsertNextCell(VTK_TETRA,4,ids);

        //T6
        ids[0]=n+1;
        ids[1]=n+shiftJ+1;
        ids[2]=n+shiftK+shiftJ;
        ids[3]=n+shiftK+shiftJ+1;
        ug->InsertNextCell(VTK_TETRA,4,ids);
      }
    }
  }

  return ug;
}

vtkUnstructuredGrid* build2dTriangulatedUnstructuredGrid(double origin[3],
    double spacing[3],
    int dimensions[3]){
  vtkUnstructuredGrid* ug=vtkUnstructuredGrid::New();
  vtkSmartPointer<vtkPoints> points=vtkSmartPointer<vtkPoints>::New();

  double p[3];
  for(int k=0; k<dimensions[2]; ++k){
    for(int j=0; j<dimensions[1]; ++j){
      for(int i=0; i<dimensions[0]; ++i){
        p[0]=origin[0]+spacing[0]*i;
        p[1]=origin[1]+spacing[1]*j;
        p[2]=origin[2]+spacing[2]*k;
        points->InsertNextPoint(p);
      }
    }
  }

  int Di{};
  int Dj{};
  if(dimensions[0]==1){
    Di=1; Dj=2;
  }
  else if(dimensions[1]==1){
    Di=0; Dj=2;
  }
  else{
    Di=0; Dj=1;
  }

  int n=0;
  vtkIdType ids[3];
  int shiftJ=dimensions[Di];
  for(int j=0; j<dimensions[Dj]-1; ++j){
    for(int i=0; i<dimensions[Di]-1; ++i){
      n=i+j*shiftJ;

      //T1
      ids[0]=n;
      ids[1]=n+1;
      ids[2]=n+shiftJ;
      ug->InsertNextCell(VTK_TRIANGLE,3,ids);

      //T2
      ids[0]=n+1;
      ids[1]=n+shiftJ+1;
      ids[2]=n+shiftJ;
      ug->InsertNextCell(VTK_TRIANGLE,3,ids);
    }
  }

  ug->SetPoints(points);

  return ug;
}

vtkImageData* buildImageData(double origin[3],
    double spacing[3],
    int dimensions[3]){
  vtkImageData* imageData=vtkImageData::New();
  vtkSmartPointer<vtkPoints> points=vtkSmartPointer<vtkPoints>::New();

  imageData->SetOrigin(origin);
  imageData->SetSpacing(spacing);
  imageData->SetDimensions(dimensions);

  return imageData;
}

int testEdgeLink3D(Triangulation& vtu_triangulation, Triangulation& vti_triangulation){
  vti_triangulation.preprocessEdges();
  vti_triangulation.preprocessEdgeLinks();
  vtu_triangulation.preprocessEdges();
  vtu_triangulation.preprocessEdgeLinks();

  cout << "3D: test getEdgeLink()...";

  const int vtu_numberOfEdges=vtu_triangulation.getNumberOfEdges();
  const int vti_numberOfEdges=vti_triangulation.getNumberOfEdges();

  if(vtu_numberOfEdges!=vti_numberOfEdges) return -1;

  vector<vector<vector<int>>> vertices0(vtu_numberOfEdges);
  vector<vector<vector<int>>> vertices1(vtu_numberOfEdges);

  auto cmp0=[](vector<int>& a, vector<int>& b){
    if(a[0]!=b[0]) return a[0]<b[0];
    else return a[1]<b[1];
  };

  auto cmp1=[](vector<vector<int>>& a, vector<vector<int>>& b){
    if(a[0][0]!=b[0][0]) return a[0][0]<b[0][0];
    else return a[0][1]<b[0][1];
  };

  for(int i=0; i<vtu_numberOfEdges; ++i){
    int v0;
    int v1;

    vector<int> tmp_v0;
    vtu_triangulation.getEdgeVertex(i, 0, v0);
    vtu_triangulation.getEdgeVertex(i, 1, v1);
    tmp_v0.push_back(v0);
    tmp_v0.push_back(v1);
    sort(tmp_v0.begin(), tmp_v0.end());
    vertices0[i].push_back(tmp_v0);
    const int vtu_linkNumber=vtu_triangulation.getEdgeLinkNumber(i);
    for(int j=0; j<vtu_linkNumber; ++j){
      vector<int> tmp_vertices0;

      int vtu_linkId;
      vtu_triangulation.getEdgeLink(i, j, vtu_linkId);
      vtu_triangulation.getEdgeVertex(vtu_linkId, 0, v0);
      vtu_triangulation.getEdgeVertex(vtu_linkId, 1, v1);
      tmp_vertices0.push_back(v0);
      tmp_vertices0.push_back(v1);

      sort(tmp_vertices0.begin(),tmp_vertices0.end());
      vertices0[i].push_back(tmp_vertices0);
    }
    sort(vertices0[i].begin()+1, vertices0[i].end(), cmp0);

    vector<int> tmp_v1;
    vti_triangulation.getEdgeVertex(i, 0, v0);
    vti_triangulation.getEdgeVertex(i, 1, v1);
    tmp_v1.push_back(v0);
    tmp_v1.push_back(v1);
    sort(tmp_v1.begin(), tmp_v1.end());
    vertices1[i].push_back(tmp_v1);
    const int vti_linkNumber=vti_triangulation.getEdgeLinkNumber(i);
    for(int j=0; j<vti_linkNumber; ++j){
      vector<int> tmp_vertices1;

      int vti_linkId;
      vti_triangulation.getEdgeLink(i, j, vti_linkId);
      vti_triangulation.getEdgeVertex(vti_linkId, 0, v0);
      vti_triangulation.getEdgeVertex(vti_linkId, 1, v1);
      tmp_vertices1.push_back(v0);
      tmp_vertices1.push_back(v1);

      sort(tmp_vertices1.begin(),tmp_vertices1.end());
      vertices1[i].push_back(tmp_vertices1);
    }
    sort(vertices1[i].begin()+1, vertices1[i].end(), cmp0);
  }
  sort(vertices0.begin(), vertices0.end(), cmp1);
  sort(vertices1.begin(), vertices1.end(), cmp1);

  for(int i=0; i<vtu_numberOfEdges; ++i){
    if(vertices0[i].size()!=vertices1[i].size()) return -1;
    const int vtu_linkNumber=vertices0[i].size();
    for(int j=0; j<vtu_linkNumber; ++j){
      if(vertices0[i][j][0]!=vertices1[i][j][0]) return -1;
      if(vertices0[i][j][1]!=vertices1[i][j][1]) return -1;
    }
  }

  cout << "OK" << endl;

  return 0;
}

int testTriangleLink3D(Triangulation& vtu_triangulation, Triangulation& vti_triangulation){
  vtu_triangulation.preprocessTriangles();
  vtu_triangulation.preprocessTriangleLinks();
  vti_triangulation.preprocessTriangles();
  vti_triangulation.preprocessTriangleLinks();

  cout << "3D: test getTriangleLink()...";
  const int vtu_numberOfTriangles=vtu_triangulation.getNumberOfTriangles();
  const int vti_numberOfTriangles=vti_triangulation.getNumberOfTriangles();

  if(vtu_numberOfTriangles!=vti_numberOfTriangles) return -1;

  vector<vector<vector<int>>> vertices0(vtu_numberOfTriangles);
  vector<vector<vector<int>>> vertices1(vtu_numberOfTriangles);

  auto cmp1=[](vector<vector<int>>& a, vector<vector<int>>& b){
    if(a[0][0]!=b[0][0]) return a[0][0]<b[0][0];
    else if(a[0][1]!=b[0][1]) return a[0][1]<b[0][1];
    else return a[0][2]<b[0][2];
  };

  for(int i=0; i<vtu_numberOfTriangles; ++i){
    int v0;
    int v1;
    int v2;

    vector<int> tmp_v0;
    vtu_triangulation.getTriangleVertex(i, 0, v0);
    vtu_triangulation.getTriangleVertex(i, 1, v1);
    vtu_triangulation.getTriangleVertex(i, 2, v2);
    tmp_v0.push_back(v0);
    tmp_v0.push_back(v1);
    tmp_v0.push_back(v2);
    sort(tmp_v0.begin(), tmp_v0.end());
    vertices0[i].push_back(tmp_v0);
    const int vtu_linkNumber=vtu_triangulation.getTriangleLinkNumber(i);
    vector<int> tmp_vertices0;
    for(int j=0; j<vtu_linkNumber; ++j){
      int vtu_linkId;
      vtu_triangulation.getTriangleLink(i, j, vtu_linkId);
      tmp_vertices0.push_back(vtu_linkId);
    }
    sort(tmp_vertices0.begin(), tmp_vertices0.end());
    vertices0[i].push_back(tmp_vertices0);

    vector<int> tmp_v1;
    vti_triangulation.getTriangleVertex(i, 0, v0);
    vti_triangulation.getTriangleVertex(i, 1, v1);
    vti_triangulation.getTriangleVertex(i, 2, v2);
    tmp_v1.push_back(v0);
    tmp_v1.push_back(v1);
    tmp_v1.push_back(v2);
    sort(tmp_v1.begin(), tmp_v1.end());
    vertices1[i].push_back(tmp_v1);
    const int vti_linkNumber=vti_triangulation.getTriangleLinkNumber(i);
    vector<int> tmp_vertices1;
    for(int j=0; j<vti_linkNumber; ++j){
      int vti_linkId;
      vti_triangulation.getTriangleLink(i, j, vti_linkId);
      tmp_vertices1.push_back(vti_linkId);
    }
    sort(tmp_vertices1.begin(), tmp_vertices1.end());
    vertices1[i].push_back(tmp_vertices1);
  }
  sort(vertices0.begin(), vertices0.end(), cmp1);
  sort(vertices1.begin(), vertices1.end(), cmp1);

  for(int i=0; i<vtu_numberOfTriangles; ++i){
    if(vertices0[i][0][0]!=vertices1[i][0][0]) return -1;
    if(vertices0[i][0][1]!=vertices1[i][0][1]) return -1;
    if(vertices0[i][0][2]!=vertices1[i][0][2]) return -1;

    const int vtu_linkNumber=vertices0[i][1].size();
    const int vti_linkNumber=vertices1[i][1].size();
    if(vtu_linkNumber!=vti_linkNumber) return -1;

    for(int j=0; j<vtu_linkNumber; ++j)
      if(vertices0[i][1][j]!=vertices1[i][1][j]) return -1;
  }

  cout << "OK" << endl;

  return 0;
}

int testEdgeLink2D(Triangulation& vtu_triangulation, Triangulation& vti_triangulation){
  vtu_triangulation.preprocessEdges();
  vtu_triangulation.preprocessEdgeLinks();
  vti_triangulation.preprocessEdges();
  vti_triangulation.preprocessEdgeLinks();

  cout << "2D: test getEdgeLink()...";
  const int vtu_numberOfEdges=vtu_triangulation.getNumberOfEdges();
  const int vti_numberOfEdges=vti_triangulation.getNumberOfEdges();

  if(vtu_numberOfEdges!=vti_numberOfEdges) return -1;

  vector<vector<vector<int>>> vertices0(vtu_numberOfEdges);
  vector<vector<vector<int>>> vertices1(vtu_numberOfEdges);

  auto cmp1=[](vector<vector<int>>& a, vector<vector<int>>& b){
    if(a[0][0]!=b[0][0]) return a[0][0]<b[0][0];
    else return a[0][1]<b[0][1];
  };

  for(int i=0; i<vtu_numberOfEdges; ++i){
    int v0;
    int v1;

    vector<int> tmp_v0;
    vtu_triangulation.getEdgeVertex(i, 0, v0);
    vtu_triangulation.getEdgeVertex(i, 1, v1);
    tmp_v0.push_back(v0);
    tmp_v0.push_back(v1);
    sort(tmp_v0.begin(), tmp_v0.end());
    vertices0[i].push_back(tmp_v0);
    const int vtu_linkNumber=vtu_triangulation.getEdgeLinkNumber(i);
    vector<int> tmp_vertices0;
    for(int j=0; j<vtu_linkNumber; ++j){
      int vtu_linkId;
      vtu_triangulation.getEdgeLink(i, j, vtu_linkId);
      tmp_vertices0.push_back(vtu_linkId);
    }
    sort(tmp_vertices0.begin(), tmp_vertices0.end());
    vertices0[i].push_back(tmp_vertices0);

    vector<int> tmp_v1;
    vti_triangulation.getEdgeVertex(i, 0, v0);
    vti_triangulation.getEdgeVertex(i, 1, v1);
    tmp_v1.push_back(v0);
    tmp_v1.push_back(v1);
    sort(tmp_v1.begin(), tmp_v1.end());
    vertices1[i].push_back(tmp_v1);
    const int vti_linkNumber=vti_triangulation.getEdgeLinkNumber(i);
    vector<int> tmp_vertices1;
    for(int j=0; j<vti_linkNumber; ++j){
      int vti_linkId;
      vti_triangulation.getEdgeLink(i, j, vti_linkId);
      tmp_vertices1.push_back(vti_linkId);
    }
    sort(tmp_vertices1.begin(), tmp_vertices1.end());
    vertices1[i].push_back(tmp_vertices1);
  }
  sort(vertices0.begin(), vertices0.end(), cmp1);
  sort(vertices1.begin(), vertices1.end(), cmp1);

  for(int i=0; i<vtu_numberOfEdges; ++i){
    if(vertices0[i][0][0]!=vertices1[i][0][0]) return -1;
    if(vertices0[i][0][1]!=vertices1[i][0][1]) return -1;

    const int vtu_linkNumber=vertices0[i][1].size();
    const int vti_linkNumber=vertices1[i][1].size();
    if(vtu_linkNumber!=vti_linkNumber) return -1;

    for(int j=0; j<vtu_linkNumber; ++j)
      if(vertices0[i][1][j]!=vertices1[i][1][j]) return -1;
  }

  cout << "OK" << endl;

  return 0;
}

int test3D(double origin[3], double spacing[3], int dimension[3]){
  vtkUnstructuredGrid* vtu=build3dTriangulatedUnstructuredGrid(origin, spacing, dimension);
  vtkImageData* vti=buildImageData(origin, spacing, dimension);

  Triangulation vtu_triangulation;
  vtu_triangulation.setInputPoints(vtu->GetNumberOfPoints(), (float*)(vtu->GetPoints()->GetVoidPointer(0)));
  vtu_triangulation.setInputCells(vtu->GetNumberOfCells(), vtu->GetCells()->GetPointer());

  Triangulation vti_triangulation;
  vti_triangulation.setInputGrid(origin[0], origin[1], origin[2],
      spacing[0], spacing[1], spacing[2],
      dimension[0], dimension[1], dimension[2]);

  int ret=0;

  ret=testEdgeLink3D(vtu_triangulation, vti_triangulation);
  if(ret){
    cout << "\nError for dimensions: " <<  dimension[0] << " " << dimension[1] << " " << dimension[2] << endl;
    return ret;
  }

  ret=testTriangleLink3D(vtu_triangulation, vti_triangulation);
  if(ret){
    cout << "\nError for dimensions : " <<  dimension[0] << " " << dimension[1] << " " << dimension[2] << endl;
    return ret;
  }

  vtu->Delete();
  vti->Delete();

  return 0;
}

int test2D(double origin[3], double spacing[3], int dimension[3]){
  vtkUnstructuredGrid* vtu=build2dTriangulatedUnstructuredGrid(origin, spacing, dimension);
  vtkImageData* vti=buildImageData(origin, spacing, dimension);

  Triangulation vtu_triangulation;
  vtu_triangulation.setInputPoints(vtu->GetNumberOfPoints(), (float*)(vtu->GetPoints()->GetVoidPointer(0)));
  vtu_triangulation.setInputCells(vtu->GetNumberOfCells(), vtu->GetCells()->GetPointer());

  Triangulation vti_triangulation;
  vti_triangulation.setInputGrid(origin[0], origin[1], origin[2],
      spacing[0], spacing[1], spacing[2],
      dimension[0], dimension[1], dimension[2]);

  int ret=0;

  ret=testEdgeLink2D(vtu_triangulation, vti_triangulation);
  if(ret){
    cout << "\nError for dimensions: " <<  dimension[0] << " " << dimension[1] << " " << dimension[2] << endl;
    return ret;
  }

  vtu->Delete();
  vti->Delete();

  return 0;
}

int main(int argc, char **argv) {
  double origin[3]{0,0,0};
  double spacing[3]{1,1,1};
  int dimension[3];

  const int limit=256;

  if(argc!=2) return 0;
  srand(time(0));

  const int dimensionality=atoi(argv[1]);

  if(dimensionality==2){
    const int i=rand()%3;
    const int j=(i+1)%3;
    const int k=(i+2)%3;
    dimension[i]=1;
    dimension[j]=(rand()%(limit-1))+2;
    dimension[k]=(rand()%(limit-1))+2;
    return test2D(origin, spacing, dimension);
  }
  else if(dimensionality==3){
    dimension[0]=(rand()%(limit-1))+2;
    dimension[1]=(rand()%(limit-1))+2;
    dimension[2]=(rand()%(limit-1))+2;
    return test3D(origin, spacing, dimension);
  }

  return 0;
}
