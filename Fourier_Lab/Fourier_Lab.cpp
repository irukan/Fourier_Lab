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

#include "MyFunc.h"
#include "NormalFunc.h"
#include "MyTimer.h"

using namespace std;

int main(int argc, const char * argv[]) {
    
    const size_t dataN = 1024;
    
    vector<double> rData(dataN);
    double dx = 2 * M_PI / dataN;
    
    for (int i=0; i< rData.size(); i++)
        rData[i] = sin(dx * i);
    
    vector<double> iData(dataN, 0);
    vector<double> rDest1, iDest1, spec1;
    vector<double> rDest2, iDest2, spec2;
    vector<double> rDest3, iDest3, spec3;
    
    TIMER.start("NormalDFT", MICRO);
    NormalDFT(rData, iData, rDest1, iDest1, spec1);
    TIMER.stop();
    
    initTable(1024);
    TIMER.start("MyDFT1", MICRO);
    MyDFT1(rData, iData, rDest2, iDest2, spec2);
    TIMER.stop();
   
  
    TIMER.start("MyDFT2", MICRO);
    MyDFT2(rData, iData, rDest3, iDest3, spec3);
    TIMER.stop();
    
    TIMER.output();
    
    return 0;
}
