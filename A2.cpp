// Primary file for the A2 code

#include <iostream>
#include <armadillo>
#include <math.h>

using namespace std;
using namespace arma;

void main(int argc, char *argv[])
{
    ifstream pointsFile;
    pointsFile.open(argv[1]);

	//After first line, each line of pointsFile is laid out as Xw, Yw, Zw, Xc, Yc

    if(pointsFile.is_open()) {
        string line;      
        getline(pointsFile, line);
		int totalLines = stoi(line.substr(0, line.find(",")));
		line.erase(0, (line.find(",") + 1));
		double centreX = stod(line.substr(0, line.find(",")));
		line.erase(0, (line.find(",") + 1));
		double centreY = stod(line.substr(0, line.find(",")));
		line.erase(0, (line.find(",") + 1));
		double pixelSize = stod(line.substr(0, line.find(",")));
		line.erase(0, (line.find(",") + 1));
		double Zed = stod(line.substr(0, line.find(",")));     
		vector<vector<double>> inputPoints;
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
			tempArray.~vector();
        }

		Col<double> XdVec(totalLines, fill::zeros);
		Col<double> Xd2Yd2Vec(totalLines, fill::zeros);
		Mat<double> MMat(0, 7);
		Row<double> BuildM(7, fill::zeros);
		uword i = 0;

		vector<vector<double>>::iterator inputPointsIterator;
		inputPointsIterator = inputPoints.begin();
		for (; inputPointsIterator != inputPoints.end(); inputPointsIterator++)
		{
			vector<double> tempArray = *inputPointsIterator;
			double Xd = pixelSize * (tempArray[3] - centreX);
			tempArray.push_back(Xd); //Place Xd onto the back of the vector, at index 5
			double Yd = pixelSize * (tempArray[4] - centreY);
			tempArray.push_back(Yd);  //Place Yd onto the back of the vector, at index 6
			*inputPointsIterator = tempArray;

			BuildM(0) = Yd * tempArray[0];
			BuildM(1) = Yd * tempArray[1];
			BuildM(2) = Yd * tempArray[2];
			BuildM(3) = Yd;
			BuildM(4) = -1.0 * Xd * tempArray[0];
			BuildM(5) = -1.0 * Xd * tempArray[1];
			BuildM(6) = -1.0 * Xd * tempArray[2];

			MMat.insert_rows(i, BuildM);
			XdVec(i) = Xd;
			Xd2Yd2Vec(i) = pow(Xd, 2.0) + pow(Yd, 2.0);
			i++;
		}

		/*Mat<double> MMatPInv = (((MMat.t() * MMat).i()) * MMat.t());
		//Mat<double> MMatPInv = pinv(MMat);
		//Mat<double> LMat = MMatPInv * XdVec;
		//Col<double> LVec(LMat);
		Col<double> LVec = MMatPInv * XdVec;*/
		double detM = det(MMat.t() * MMat);
		Col<double> LVec = (((MMat.t() * MMat).i()) * MMat.t()) * XdVec;


		double Ty;
		double AbsTy = abs(1.0 / (sqrt(pow(LVec(4), 2.0) + pow(LVec(5), 2.0) + pow(LVec(6), 2.0))));
		double Sx = AbsTy * sqrt(pow(LVec(0), 2.0) + pow(LVec(1), 2.0) + pow(LVec(2), 2.0));

		double TempRs[7];
		TempRs[0] = AbsTy * LVec(0); // Corresponds to temporary r11
		TempRs[1] = AbsTy * LVec(1); // Corresponds to temporary r12
		TempRs[2] = AbsTy * LVec(2); // Corresponds to temporary r13
		TempRs[3] = AbsTy * LVec(3); // Corresponds to temporary tx
		TempRs[4] = AbsTy * LVec(4); // Corresponds to temporary r21
		TempRs[5] = AbsTy * LVec(5); // Corresponds to temporary r22
		TempRs[6] = AbsTy * LVec(6); // Corresponds to temporary r23

		uword LargestIndex;
		Xd2Yd2Vec.max(LargestIndex);
		vector<double> LargestPoint = inputPoints.at(LargestIndex);
		
		double TempXt = TempRs[0] * LargestPoint[0] + TempRs[1] * LargestPoint[1] + TempRs[2] * LargestPoint[2] + TempRs[3];
		double TempYt = TempRs[4] * LargestPoint[0] + TempRs[5] * LargestPoint[1] + TempRs[6] * LargestPoint[2] + AbsTy;

		if ((signbit(TempXt) == signbit(LargestPoint[5])) && (signbit(TempYt) == signbit(LargestPoint[6]))) {
			Ty = AbsTy;
		}
		else
		{
			Ty = AbsTy * -1;
		}

		Mat<double> Rs(2,3, fill::zeros); double Tx; //Only two rows for Rs for now, because the third row is inserted later after the cross-product stage
		Rs(0,0) = (Ty * LVec(0)) / Sx; // Corresponds to r11
		Rs(0,1) = (Ty * LVec(1)) / Sx; // Corresponds to r12
		Rs(0,2) = (Ty * LVec(2)) / Sx; // Corresponds to r13
		Rs(1,0) = Ty * LVec(4); // Corresponds to r23
		Rs(1,1) = Ty * LVec(5); // Corresponds to r21
		Rs(1,2) = Ty * LVec(6); // Corresponds to r22
		Tx = (Ty * LVec(3)) / Sx; // Corresponds to tx

		Row<double> R3 = cross(Rs.row(0), Rs.row(1));
		Rs.insert_rows(2, R3);

		double dotR1nR2 = dot(Rs.row(0), Rs.row(1));

		Mat<double> UyYd(totalLines, 2, fill::zeros);
		Col<double> Uz(totalLines, fill::zeros);

		i = 0;
		inputPointsIterator = inputPoints.begin();
		for (; inputPointsIterator != inputPoints.end(); inputPointsIterator++) {
			vector<double> tempArray = *inputPointsIterator;
			UyYd(i, 0) = Rs(1, 0) * tempArray[0] + Rs(1, 1) * tempArray[1] + Rs(1, 2) * tempArray[2] + Ty;
			UyYd(i, 1) = tempArray[6] * -1.0;
			Uz(i) = (Rs(2, 0) * tempArray[0] + Rs(2, 1) * tempArray[1] + Rs(2, 2) * tempArray[2]) * tempArray[6];
			i++;
		}

		Col<double> fandTz(2, fill::zeros);
		fandTz = pinv(UyYd) * Uz;

		ofstream OutputToFile("CalibrationResults.csv", ios::trunc);
		OutputToFile << "L1, L2, L3, L4, L5, L6, L7" << endl;
		OutputToFile << LVec(0) << "," << LVec(1) << "," << LVec(2) << "," << LVec(3) << "," << LVec(4) << "," << LVec(5) << "," << LVec(6) << endl << endl;
		OutputToFile << "r11, r12, r13, r21, r22, r23, r31, r32, r33" << endl;
		OutputToFile << Rs(0,0) << "," << Rs(0, 1) << "," << Rs(0, 2) << "," << Rs(1, 0) << "," << Rs(1, 1) << "," << Rs(1, 2) << "," 
			<< Rs(2, 0) << "," << Rs(2, 1) << "," << Rs(2, 2) << endl << endl;
		OutputToFile << "Tx, Ty, Tz, f, Sx, Det M, pixelSize, R1.R2" << endl;
		OutputToFile << Tx << "," << Ty << "," << fandTz(1) << "," << fandTz(0) << "," << Sx << "," << detM << "," << pixelSize << "," << dotR1nR2 << endl << endl;

		// Generate estimated image coordinates based on measured world coordinates and generated values
		// Produce pixel error estimates, and then use those to calculate first estimate for Kappa1
		// Sc =  RSw + T

		Col<double> Ts(3);
		Ts(0) = Tx;
		Ts(1) = Ty;
		Ts(2) = fandTz(1);
		i = 0;
		inputPointsIterator = inputPoints.begin();
		Col<double> Sw(3);
		Col<double> Sc(3);
		OutputToFile << "Xw, Yw, Zw, Xcm, Ycm, Xce, Yce, Xcm - Xce, Ycm - Yce, ErrMagnitude, Kappa1 estimate, Xd2Yd2" << endl;
		for (; inputPointsIterator != inputPoints.end(); inputPointsIterator++) 
		{
			vector<double> tempArray = *inputPointsIterator;
			Sw(0) = tempArray[0];
			Sw(1) = tempArray[1];
			Sw(2) = tempArray[2];
			Sc = (Rs * Sw) + Ts;
			double Xu = (fandTz(0) / pixelSize) * (Sc(0) / Sc(2));
			double Yu = (fandTz(0) / pixelSize) * (Sc(1) / Sc(2));
			double Cx = (Sx * Xu) + centreX;
			double Cy = (Yu) + centreY;
			double errXpix = tempArray[3] - Cx;
			double errYpix = tempArray[4] - Cy;
			double errMagnitude = sqrt(pow(errXpix, 2.0) + pow(errYpix, 2.0));
			double Kappa1 = (errMagnitude * pixelSize) / (pow(sqrt(Xd2Yd2Vec(i)), 3.0));
			tempArray.push_back(Cx); // Store the estimated camera coordinate in inputPoints index 7
			tempArray.push_back(Cy); // Store the estimated camera coordinate in inputPoints index 8
			tempArray.push_back(errXpix); // Store the pixel error in inputPoints index 9
			tempArray.push_back(errYpix); // Store the pixel error in inputPoints index 10
			tempArray.push_back(errMagnitude); // Store the pixel error magnitude in inputPoints index 11
			tempArray.push_back(Kappa1); // Store Kappa1 estimate in inputPoints index 12
			*inputPointsIterator = tempArray; // Redefine inputPoints(i) with the extra elements
			OutputToFile << Sw(0) << "," << Sw(1) << "," << Sw(2) << "," << tempArray[3] << "," << tempArray[4] << "," << Cx << "," << Cy << "," << errXpix << "," << errYpix << 
				"," << errMagnitude << "," << Kappa1 << "," << Xd2Yd2Vec(i) << endl;
			i++;
		}

		// Generate estimated world coordinates based on measured pixel values

		/*i = 0;
		inputPointsIterator = inputPoints.begin();
		OutputToFile << endl << "Xwe, Ywe, Zwe, Xi, Yi" << endl;

		for (; inputPointsIterator != inputPoints.end(); inputPointsIterator++)
		{
			vector<double> tempArray = *inputPointsIterator;
			
		}*/

		OutputToFile.close();
    }
    pointsFile.close();
	cout << "Success!";
}