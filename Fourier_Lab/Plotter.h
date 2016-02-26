//
//  Plotter.h
//  Fourier_Lab
//
//  Created by kayama on 2016/02/24.
//  Copyright © 2016年 kayama. All rights reserved.
//

#ifndef Plotter_h
#define Plotter_h
#include <fstream>
#include <math.h>
#include "MyTaylorFunc.h"
#include "MyTableFunc.h"
#include "NormalFunc.h"
void Plotter()
{
    size_t dataN = 1024;
    double dx = 2 * M_PI / dataN;
    initTable(1024);
    
    ofstream ofs("output.csv");
    
    vector<double> rData(dataN, 0);
    for (int i=dataN/2; i < dataN/2 + 15; i++)
    {
        //rData[i] = sin(dx * i *100.1234) + sin(dx * i *101.1234) + sin(dx * i *99.1234);
        rData[i] = 1;
    }

    
    vector<double> iData(dataN, 0);
    vector<double> rDest1, iDest1, spec1;
    vector<double> rDest2, iDest2, spec2;
    vector<double> rDest3, iDest3, spec3;
    
    
    NormalDFT(rData, iData, rDest1, iDest1, spec1);
    MyTaylorDFT(rData, iData, rDest2, iDest2, spec2);
    
    for (size_t i = 0; i< dataN ; i++)
    {
        double RAD = dx * i / M_PI;
        
        double theta = dx * i;
        int signSin, signCos;
        rotateTheta(&theta, &signSin, &signCos);
        double myCosData = myTableCos(theta);
        double Cos = cos(dx * i);
        double mySinData = myTaylorSin(theta, signSin);
        double Sin = sin(dx * i);
        double diff = fabs(spec1[i] - spec2[i]);
        
        ofs
        << "," << RAD
        << "," << rData[i]
        << "," << myCosData
        << "," << Cos
        << "," << mySinData
        << "," << Sin
        << "," << diff
        << "," << spec1[i]
        << "," << spec2[i]
        << endl;
        
    }
    
    ofs.close();
    system("python Plot.py output.csv");
}

#endif /* Plotter_h */
