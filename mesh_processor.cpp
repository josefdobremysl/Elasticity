#include "mesh_processor.h"
#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <tuple>


MatrixVectorResult GetTriangles(const string& meshFile) {
    MatrixVectorResult result;
    ifstream file(meshFile);

    if (!file.is_open()) {
        cerr << "Nelze otevřít soubor " << meshFile << endl;
        return {};
    }

    string line;
    bool isReadingElements = false;

    vector < vector<double>> nodes_coordinates;
    int i = 0;

    while (getline(file, line)) {
        if (line.find("$Nodes") != string::npos) {
            isReadingElements = true;
            continue;
        }

        if (line.find("$EndNodes") != string::npos) {
            break;
        }

        if (isReadingElements) {
            istringstream iss(line);
            int nodeId;
            double cooX, cooY, cooZ;
            vector<double> coordinates(3);
            iss >> nodeId >> cooX >> cooY >> cooZ;


            if (i > 0) {
                nodes_coordinates.push_back({ cooX, cooY, cooZ });
            }
            i++;
        }
    }
    //cout << "  " << endl;
    //cout << " nodes coordinates: " << endl;

    //for (const auto& row : nodes_coordinates) {    // vytiskne matici
    //    for (const double value : row) {
    //        std::cout << value << ' ';
    //    }
    //    cout << endl;
    //}

    isReadingElements = false; // Resetujeme načítání prvků
    file.clear(); // Resetujeme stav souboru
    file.seekg(0, std::ios::beg); // Nastavíme pozici na začátek souboru




    vector<int> element_type_line_V;      // příprav vektorů pro úsečky
    vector<int> element_id_line_V;
    vector<int> num_tags_line_V;
    vector<int> physical_tag_line_V;
    vector<int> elementary_tag_line_V;

    vector<int> element_type_tri_V;       // příprav vektorů pro trojúhelníky
    vector<int> element_id_tri_V;
    vector<int> num_tags_tri_V;
    vector<int> physical_tag_tri_V;
    vector<int> elementary_tag_tri_V;

    vector<vector<int>> triangle_vertex_M;
    vector<vector<int>> line_vertex_M;



    int skip_first = 0;

    while (getline(file, line)) {
        if (line.find("$Elements") != string::npos) {
            isReadingElements = true;
            continue;
        }

        if (line.find("$EndElements") != string::npos) {
            break;
        }

        if (isReadingElements) {
            istringstream iss(line);
            int elementType, elementId, numTags, physicalTag, elementaryTag, nodeA, nodeB, nodeC;
            iss >> elementId >> elementType >> numTags >> physicalTag >> elementaryTag >> nodeA >> nodeB >> nodeC;

            //std::cout << elementId << " " << elementType << " " << numTags << " " << physicalTag << " " << elementaryTag << " " << nodeA << " " << nodeB << " " << nodeC << " " << std::endl;


            if (skip_first > 0) {

                if (elementType == 1) {                              // 1 je označení úseček
                    elementary_tag_line_V.push_back(elementaryTag);
                    element_type_line_V.push_back(elementType);
                    element_id_line_V.push_back(elementId);
                    num_tags_line_V.push_back(numTags);
                    physical_tag_line_V.push_back(physicalTag);
                    elementary_tag_line_V.push_back(elementaryTag);
                    line_vertex_M.push_back({ nodeA, nodeB, nodeC });
                }
                if (elementType == 2) {                              // 2 je označení trojúhelníků
                    elementary_tag_tri_V.push_back(elementaryTag);
                    element_type_tri_V.push_back(elementType);
                    element_id_tri_V.push_back(elementId);
                    num_tags_tri_V.push_back(numTags);
                    physical_tag_tri_V.push_back(physicalTag);
                    elementary_tag_tri_V.push_back(elementaryTag);
                    triangle_vertex_M.push_back({ nodeA, nodeB, nodeC });
                }
            }
            skip_first++;



        }

    }
    for (int i = 0; i < physical_tag_line_V.size(); i++) {           // Tiskne vektor
        //std::cout << physical_tag_line_V[i] << ' ';
    }



    result.matrix1 = triangle_vertex_M;          // ukládání výstupu funkce
    result.matrix2 = line_vertex_M;
    result.matrix3 = nodes_coordinates;
    result.vector1 = physical_tag_line_V;
    result.vector2 = physical_tag_tri_V;

    int n_Tri = physical_tag_tri_V.size();      // počet trojúhelníků
    int n_Line = physical_tag_line_V.size();    // počet úseček

    result.integerResult = n_Tri;
    result.integerResult2 = n_Line;



    file.close();
    return result;
}
