//
//  PlotData.h
//  Fourier_Lab
//
//  Created by kayama on 2016/02/24.
//  Copyright © 2016年 kayama. All rights reserved.
//

#ifndef PlotData_h
#define PlotData_h




void PlotData()
{
    int dataN = 1024;
    double dx = 2 * M_PI / dataN;
    initTable(dataN);
    
    
    vector<double> Source(dataN, 0);
    for (int i=dataN/2; i < dataN/2 + 15; i++)
    {
        //Source[i] = sin(dx * i *100.1234) + sin(dx * i *101.1234) + sin(dx * i *99.1234);
        Source[i] = 1;
    }
//    for (int i =0; i < dataN; i++)
//    {
//        Source[i] = sin(dx * i * 100.1234);
//    }
    PLOTER.SetData("Source", Source);
    
    vector<double> Rad, Cos, Sin;
    for (size_t i = 0; i< dataN ; i++)
    {
        PLOTER.SetData("Rad", dx * i * M_PI);
        PLOTER.SetData("Cos", cos(dx *i));
        PLOTER.SetData("Sin", sin(dx *i));
    }
    
    vector<double> iData(dataN, 0);
    vector<double> rDest, iDest, spec;
    
    
    NormalDFT(Source, iData, rDest, iDest, spec);
    PLOTER.SetData("NormalDFT-Real", rDest);
    PLOTER.SetData("NormalDFT-Imag", iDest);
    PLOTER.SetData("NormalDFT-Spec", spec);
    
    NormalDFT2(Source, iData, rDest, iDest, spec);
    PLOTER.SetData("NormalDFT2-Real", rDest);
    PLOTER.SetData("NormalDFT2-Imag", iDest);
    PLOTER.SetData("NormalDFT2-Spec", spec);

    MyTableDFT(Source, iData, rDest, iDest, spec);
    PLOTER.SetData("MyTableDFT-Real", rDest);
    PLOTER.SetData("MyTableDFT-Imag", iDest);
    PLOTER.SetData("MyTableDFT-Spec", spec);
    
    MyTableDFT_SSE(Source, iData, rDest, iDest, spec);
    PLOTER.SetData("MyTableDFT_SSE-Real", rDest);
    PLOTER.SetData("MyTableDFT_SSE-Imag", iDest);
    PLOTER.SetData("MyTableDFT_SSE-Spec", spec);

    MyTaylorDFT(Source, iData, rDest, iDest, spec);
    PLOTER.SetData("MyTaylorDFT-Real", rDest);
    PLOTER.SetData("MyTaylorDFT-Imag", iDest);
    PLOTER.SetData("MyTaylorDFT-Spec", spec);
    
    MyTaylorDFT2(Source, iData, rDest, iDest, spec);
    PLOTER.SetData("MyTaylorDFT2-Real", rDest);
    PLOTER.SetData("MyTaylorDFT2-Imag", iDest);
    PLOTER.SetData("MyTaylorDFT2-Spec", spec);
    
    NormalFFT(Source, iData, rDest, iDest, spec);
    PLOTER.SetData("NormalFFT-Real", rDest);
    PLOTER.SetData("NormalFFT-Imag", iDest);
    PLOTER.SetData("NormalFFT-Spec", spec);
   
}



#endif /* PlotData_h */
