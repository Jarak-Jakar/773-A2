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
		int totalLines = stoi(line.substr(0, line.find(",")));
		line.erase(0, (line.find(",") + 1));
		double nummerA = stod(line.substr(0, line.find(",")));
		line.erase(0, (line.find(",") + 1));
		double nummerB = stod(line.substr(0, line.find(",")));
		line.erase(0, (line.find(",") + 1));
		double pixelSize = stod(line.substr(0, line.find(",")));
		//line.erase(0, (line.find(",") + 1));
		cout << "pixelSize = " << pixelSize << endl;
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

		/*vector<vector<double>>::iterator outerIter;
		outerIter = inputPoints.begin();

		for (; outerIter != inputPoints.end(); outerIter++)
		{
			vector<double> tempArray;
			tempArray = *outerIter;
			vector<double>::iterator innerIter;
			innerIter = tempArray.begin();
			for (; innerIter != tempArray.end(); innerIter++)
			{
				cout << "Element: " << *innerIter << endl;
			}

			cout << "End of entry" << endl << endl;
		}*/

		vector<vector<double>>::iterator grabXnY;
		cout << "Declared grab" << endl;
		grabXnY = inputPoints.begin();
		cout << "Initialised grab" << endl;
		for (; grabXnY != inputPoints.end(); grabXnY++)
		{
			cout << "Entered for loop" << endl;
			vector<double> tempArray = *grabXnY;
			cout << "Dereferenced grab" << endl;
			double temp;
			double Xd; // = tempArray[3] * pixelSize;
			temp = tempArray[3];
			Xd = temp * pixelSize;
			cout << "Made Xd" << endl;
			double Yd; // = tempArray[4] * pixelSize;
			temp = tempArray[4];
			Yd = temp * pixelSize;
			cout << "Made Yd" << endl;
			cout << "Xd & Yd are: " << Xd << ", " << Yd << endl;
		}
    }
    pointsFile.close();
}