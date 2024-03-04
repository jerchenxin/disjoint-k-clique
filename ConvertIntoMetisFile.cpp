// g++ -O3 -std=c++17 ConvertIntoMetisFile.cpp -o ConvertIntoMetisFile

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

void convertToMetisGraph(const std::string& inputFile, const std::string& outputFile) {
    std::ifstream inFile(inputFile);
    std::ofstream outFile(outputFile);

    if (!inFile || !outFile) {
        std::cout << "Failed to open files!" << std::endl;
        return;
    }

    int numVertices, numEdges;
    inFile >> numVertices >> numEdges;

    // Write the header information to the output file
    outFile << numVertices << " " << numEdges << " 0 0" << std::endl;

    std::vector<std::vector<int>> adjacencyList(numVertices);

    int vertex1, vertex2;
    while (inFile >> vertex1 >> vertex2) {
        // Adjust vertex numbering to 0-based index
        // vertex1--;
        // vertex2--;
        
        // Add edges to the adjacency list
        adjacencyList[vertex1].push_back(vertex2);
        adjacencyList[vertex2].push_back(vertex1);
    }

    // Write the adjacency list to the output file
    for (const auto& neighbors : adjacencyList) {
        for (const auto& neighbor : neighbors) {
            outFile << neighbor << " ";
        }
        outFile << std::endl;
    }

    inFile.close();
    outFile.close();

    std::cout << "Conversion completed successfully!" << std::endl;
}

int main(int argc, char* argv[]) {
    std::string inputFile(argv[1]); 
    std::string outputFile(argv[2]);

    convertToMetisGraph(inputFile, outputFile);

    return 0;
}