#include <iostream>
#include <fstream>
#include "vtk.h"
#include <iomanip> // pro std::setprecision

void VTKfile(const std::string& filename, const string& meshFile, Eigen::VectorXd& solution_x) {
    MatrixVectorResult result = GetTriangles(meshFile); //naètení sítì
    vector<vector<double>> points = result.matrix3;
    vector<vector<int>> triangles = result.matrix1;
    int nTri = result.integerResult;
    Eigen::VectorXd x = solution_x;

    // Otevøení souboru pro zápis
    ofstream vtkFile(filename.c_str());

   // Zápis hlavièky VTK souboru
    vtkFile << "# vtk DataFile Version 3.0" << endl;
    vtkFile << "    " << endl;
    vtkFile << "ASCII" << endl;
    vtkFile << "DATASET UNSTRUCTURED_GRID" << endl;


    int numPoints = result.matrix3.size();    // Poèet bodù

    vtkFile << "POINTS " << numPoints << " float" << endl;  // Zápis informací o bodech
    for (int i = 0; i < numPoints; ++i) {
        vtkFile << points[i][0] << " " << points[i][1] << " " << points[i][2] << endl;
    }

    // Zápis informací o buòkách 
    vtkFile << "CELLS " << nTri << " " << nTri * 4 << endl;
    //vtkFile << numPoints;
    for (int i = 0; i < nTri; ++i) {
        vtkFile << 3 << " " << triangles[i][0] - 1 << " " << triangles[i][1] - 1 << " " << triangles[i][2] - 1 << endl;
    }


    // Zápis typu buòky
    vtkFile << "CELL_TYPES " << nTri << endl;
    for (int i = 0; i < nTri; ++i) {
        vtkFile << 5 << endl;
    }
    // Zápis hodnot øešení 
    vtkFile << "POINT_DATA " << numPoints << endl;

    //vtkFile << "SCALARS scalars float 1" << endl;
    //vtkFile << "LOOKUP_TABLE default" << endl;
    //for (int i = 0; i < numPoints; ++i) {
    //    vtkFile << x[i] << endl;
    //}
    
    // Nastavit počet desetinných míst
    vtkFile << std::fixed << std::setprecision(20);
    

    vtkFile << "VECTORS displacement double" << endl;
    for (int i = 0; i < numPoints; ++i) {
        vtkFile << x[i] << " " << x[i + numPoints] << " " << 0.0 << endl;
    }


    // Závìr souboru
    vtkFile.close();

    cout << "Soubor " << filename << " byl vytvoøen." << endl;
}


