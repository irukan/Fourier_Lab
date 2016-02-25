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
#include "MyFunc.h"
#include "NormalFunc.h"
void Plotter()
{
    size_t dataN = 1024;
    double dx = 2 * M_PI / dataN;
    
    ofstream ofs("output.csv");
    
    vector<double> rData(dataN);
    for (int i=0; i< rData.size(); i++)
        rData[i] = sin(dx * i * 100) + cos(dx * i * 200);
    vector<double> iData(dataN, 0);
    vector<double> rDest1, iDest1, spec1;
    vector<double> rDest2, iDest2, spec2;
    vector<double> rDest3, iDest3, spec3;
    
    
    NormalDFT(rData, iData, rDest1, iDest1, spec1);
    MyDFT2(rData, iData, rDest2, iDest2, spec2);
    
    for (size_t i = 0; i< dataN; i++)
    {
        double X = dx * i / M_PI;
        double myCos = myTaylorCos(dx * i);
        double Cos = cos(dx * i);
        double mySin = myTaylorSin(dx * i);
        double Sin = sin(dx * i);
        double diff = fabs(spec1[i] - spec2[i]);
        
        ofs
        << "," << X
        << "," << myCos
        << "," << Cos
        << "," << mySin
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
