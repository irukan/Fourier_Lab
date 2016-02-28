//
//  PlotData.h
//  Fourier_Lab
//
//  Created by kayama on 2016/02/24.
//  Copyright © 2016年 kayama. All rights reserved.
//

#ifndef PlotData_h
#define PlotData_h
#include <fstream>
#include <math.h>
#include <map>
#include <vector>
#include <string>
#include <sstream>
#include "MyTaylorFunc.h"
#include "MyTableFunc.h"
#include "NormalFunc.h"

using namespace::std;
class CSVOutput
{
private:
    map<string, vector<double>&> m_data;
    size_t maxDataN;
public:
    CSVOutput()
    {
        maxDataN = 0;
    }
    
    void
    SetData(string name, vector<double>& data)
    {
        //m_data[name] = data;
        m_data.insert( map<string, vector<double>&>::value_type( name, data ) );
        maxDataN = maxDataN < data.size() ? data.size() : maxDataN;
    }
    
    void
    OutputData(string fileName)
    {
        ofstream ofs(fileName);
        
        map<string, vector<double>&>::iterator itr_data;
        
        // key名の行　出力
        for (itr_data = m_data.begin(); itr_data != m_data.end(); ++itr_data)
        {
            if (itr_data != m_data.begin())
                ofs << ",";
            
            ofs << itr_data->first;
        }
        ofs << endl;
        
        // データの出力
        for (size_t i = 0; i < maxDataN; i++)
        {
            ostringstream oss;
            for (itr_data = m_data.begin(); itr_data != m_data.end(); ++itr_data)
            {
                if (itr_data != m_data.begin())
                    oss << ",";
                
                if ((itr_data->second).size() > i)
                    oss << (itr_data->second)[i];
            }
            oss << endl;
            
            ofs << oss.str();
        }
        // ファイル出力
        ofs.close();
    }
};


void Plotter()
{
    size_t dataN = 1024;
    double dx = 2 * M_PI / dataN;
    initTable(dataN);
    
    CSVOutput csv;
    
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
    csv.SetData("Source", Source);
    
    vector<double> Rad, Cos, Sin;
    for (size_t i = 0; i< dataN/2 ; i++)
    {
        Rad.push_back(dx * i / M_PI);
        Cos.push_back(cos(dx * i));
        Sin.push_back(sin(dx * i));
    }
    csv.SetData("Rad", Rad);
    csv.SetData("Cos", Cos);
    csv.SetData("Sin", Sin);
    
    vector<double> iData(dataN, 0);
    vector<double> rDest, iDest, spec;
    
    
    NormalDFT(Source, iData, rDest, iDest, spec);
    csv.SetData("NormalDFT-Real", rDest);
    csv.SetData("NormalDFT-Imag", iDest);
    csv.SetData("NormalDFT-Spec", spec);
    
    NormalDFT2(Source, iData, rDest, iDest, spec);
    csv.SetData("NormalDFT2-Real", rDest);
    csv.SetData("NormalDFT2-Imag", iDest);
    csv.SetData("NormalDFT2-Spec", spec);

    MyTableDFT(Source, iData, rDest, iDest, spec);
    csv.SetData("MyTableDFT-Real", rDest);
    csv.SetData("MyTableDFT-Imag", iDest);
    csv.SetData("MyTableDFT-Spec", spec);
    
    MyTableDFT_SSE(Source, iData, rDest, iDest, spec);
    csv.SetData("MyTableDFT_SSE-Real", rDest);
    csv.SetData("MyTableDFT_SSE-Imag", iDest);
    csv.SetData("MyTableDFT_SSE-Spec", spec);

    MyTaylorDFT(Source, iData, rDest, iDest, spec);
    csv.SetData("MyTaylorDFT-Real", rDest);
    csv.SetData("MyTaylorDFT-Imag", iDest);
    csv.SetData("MyTaylorDFT-Spec", spec);
    
    
    csv.OutputData("output.csv");
    system("python Plot.py output.csv");

}



#endif /* PlotData_h */
