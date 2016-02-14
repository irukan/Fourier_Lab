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

using namespace std;

int main(int argc, const char * argv[]) {
    
    vector<double> data(1024);
    double dx = 2 * M_PI / data.size();
    
    for (int i=0; i< data.size(); i++)
        data[i] = sin(dx * i);
    
    
    return 0;
}
