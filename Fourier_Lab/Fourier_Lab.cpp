//
//  Fourier_Lab.cpp
//  Fourier_Lab
//
//  Created by kayama on 2016/02/14.
//  Copyright © 2016年 kayama. All rights reserved.
//

#include <iostream>
#include <vector>
#include <math.h>

#include "MyTableFunc.h"
#include "MyTaylorFunc.h"
#include "NormalFunc.h"
#include "MyTimer.h"
#include "Plotter.h"

using namespace std;

int main(int argc, const char * argv[]) {
    //Plotter();
    
    const size_t dataN = 1024;
    const size_t loopN = 1000;
    vector<double> rData(dataN);
    double dx = 2 * M_PI / dataN;
    
    for (int i=0; i< rData.size(); i++)
        rData[i] = sin(dx * i) + cos(dx * 1.2 * i);
    vector<double> iData(dataN, 0);
    
    initTable(dataN);
    vector<double> rDest1(dataN), iDest1(dataN), spec1(dataN);
    vector<double> rDest2(dataN), iDest2(dataN), spec2(dataN);
    vector<double> rDest3(dataN), iDest3(dataN), spec3(dataN);
    
    for (size_t lp = 0; lp < loopN; lp++)
    {
        TIMER.start("NormalDFT", MICRO);
        NormalDFT(rData, iData, rDest1, iDest1, spec1);
        TIMER.stop();
        
        TIMER.start("MyTaylorDFT", MICRO);
        MyTaylorDFT(rData, iData, rDest2, iDest2, spec2);
        TIMER.stop();
      
        TIMER.start("MyTableDFT", MICRO);
        MyTableDFT(rData, iData, rDest3, iDest3, spec3);
        TIMER.stop();
    }
    
    TIMER.output("execTime.csv");
    
    return 0;
}
