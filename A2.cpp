// Primary file for the A2 code

#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

void main(int argc, char *argv[])
{
    //cout << "Hello, world, from Visual C++!" << endl;
    ifstream pointsFile;
    pointsFile.open(argv[1]);
    if(pointsFile.is_open()) {
        string line;      
        getline(pointsFile, line);
        int totalLines = stoi(line);      
        vector<vector<double>> inputPoints(totalLines);      
        for(int i = 0; i < totalLines; i++) {
            getline(pointsFile, line);
            vector<double> tempArray;
            for(int j = 0; j < 5; j++) {
                string token = line.substr(0, line.find(","));
                double tempVal = stod(token);
                tempArray.push_back(tempVal);
                line.erase(0, (line.find(",") + 1));
            }
            inputPoints.push_back(tempArray);
        }
        //cout << inputPoints.front().front() <<endl;
    }
    pointsFile.close();
}