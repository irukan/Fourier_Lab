//
//  ExecTime.h
//  Fourier_Lab
//
//  Created by kayama on 2016/02/27.
//  Copyright © 2016年 kayama. All rights reserved.
//

#ifndef ExecTime_h
#define ExecTime_h

void ExecTime(int dataN, int loopN)
{
    initTable(dataN);
    vector<double> rData(dataN);
    double dx = 2 * M_PI / dataN;
    
    for (int i=0; i< rData.size(); i++)
        rData[i] = sin(dx * i) + cos(dx * 1.2 * i);
    vector<double> iData(dataN, 0);
    
    
    vector<double> rDest(dataN), iDest(dataN), spec(dataN);
    
    for (size_t lp = 0; lp < loopN; lp++)
    {
//        TIMER.start("NormalDFT", MICRO);
//        NormalDFT(rData, iData, rDest, iDest, spec);
//        TIMER.stop();
//        
//        TIMER.start("NormalFFT", MICRO);
        NormalFFT(rData, iData, rDest, iDest, spec);
//        TIMER.stop();

//        TIMER.start("NormalDFT2", MICRO);
//        NormalDFT2(rData, iData, rDest, iDest, spec);
//        TIMER.stop();
//
//        TIMER.start("MyTaylorDFT", MICRO);
//        MyTaylorDFT(rData, iData, rDest, iDest, spec);
//        TIMER.stop();
//
//        TIMER.start("MyTaylorDFT2", MICRO);
//        MyTaylorDFT2(rData, iData, rDest, iDest, spec);
//        TIMER.stop();
//        
//        TIMER.start("MyTableDFT", MICRO);
//        MyTableDFT(rData, iData, rDest, iDest, spec);
//        TIMER.stop();

        TIMER.start("MyTableDFT_SSE", MICRO);
        MyTableDFT_SSE(rData, iData, rDest, iDest, spec);
        TIMER.stop();
//
        TIMER.start("MyTableDFT_SSE2", MICRO);
        MyTableDFT_SSE2(rData, iData, rDest, iDest, spec);
        TIMER.stop();
        
        TIMER.start("FFT", MICRO);
        FFT(rData, iData, rDest, iDest, spec);
        TIMER.stop();
    }
    
    TIMER.output("execTime.csv");
}

#endif /* ExecTime_h */
