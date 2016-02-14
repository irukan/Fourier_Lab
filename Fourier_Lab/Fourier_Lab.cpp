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

using namespace std;

int main(int argc, const char * argv[]) {
    
    const size_t dataN = 1024;
    
    vector<double> rData(dataN);
    double dx = 2 * M_PI / dataN;
    
    for (int i=0; i< rData.size(); i++)
        rData[i] = sin(dx * i);
    
    vector<double> iData(dataN, 0);
    vector<double> rDest, iDest, spec;
    
    NormalDFT(rData, iData, rDest, iDest, spec);
    
    
    
    return 0;
}
