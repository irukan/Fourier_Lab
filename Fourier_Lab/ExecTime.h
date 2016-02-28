//
//  ExecTime.h
//  Fourier_Lab
//
//  Created by kayama on 2016/02/27.
//  Copyright © 2016年 kayama. All rights reserved.
//

#ifndef ExecTime_h
#define ExecTime_h

void ExecTime()
{
    const size_t dataN = 1024;
    initTable(dataN);
    const size_t loopN = 500;
    vector<double> rData(dataN);
    double dx = 2 * M_PI / dataN;
    
    for (int i=0; i< rData.size(); i++)
        rData[i] = sin(dx * i) + cos(dx * 1.2 * i);
    vector<double> iData(dataN, 0);
    
    
    vector<double> rDest1(dataN), iDest1(dataN), spec1(dataN);
    vector<double> rDest2(dataN), iDest2(dataN), spec2(dataN);
    vector<double> rDest3(dataN), iDest3(dataN), spec3(dataN);
    vector<double> rDest4(dataN), iDest4(dataN), spec4(dataN);
    vector<double> rDest5(dataN), iDest5(dataN), spec5(dataN);
    
    for (size_t lp = 0; lp < loopN; lp++)
    {
        TIMER.start("NormalDFT", MICRO);
        NormalDFT(rData, iData, rDest1, iDest1, spec1);
        TIMER.stop();

        TIMER.start("NormalDFT2", MICRO);
        NormalDFT(rData, iData, rDest5, iDest5, spec5);
        TIMER.stop();

        TIMER.start("MyTaylorDFT", MICRO);
        MyTaylorDFT(rData, iData, rDest2, iDest2, spec2);
        TIMER.stop();
        
        TIMER.start("MyTableDFT", MICRO);
        MyTableDFT(rData, iData, rDest3, iDest3, spec3);
        TIMER.stop();
        
        TIMER.start("MyTableDFT_SSE", MICRO);
        MyTableDFT_SSE(rData, iData, rDest4, iDest4, spec4);
        TIMER.stop();
    }
    
    TIMER.output("execTime.csv");
}

#endif /* ExecTime_h */
