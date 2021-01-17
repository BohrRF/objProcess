#include <iostream>
#include <list>
#include "objpoints.h"
#include "pointdata.h"

const string EXTRACT = "output_120_35_0.3lr_th0.1279_extract.txt";
//const string AFTEREF = "MC.txt";
const string AFTEREF = "origin.txt";
const string BTRACE = "FFHR-d1A_R3.50_20181105_extract.txt";



int main()
{
    CPointDataBase data;
    CMethod m;
    CSlicePoints Ex;
    Cslice Af;

    data.addPointFromFile(BTRACE);
    //Af.importFromObj(AFTEREF);
    int Seg[] = {45};
    //int Seg[] = {0, 1, 2, 5, 10, 20, 45, 90, 180};
    for(int Deg: Seg)
    {
        //Af.Angle = Deg;
        //Af.findPointOnSlice();
        //Af.outSlice("Objslice");

        Ex.Angle = Deg;
        Ex.findPointsOnSlice(data);
        Ex.outSliceFile();

        Ex.addLar();
        m.mapCreate(Ex);
        m.convexHullFind();
        m.method();

        m.compareOut("AFCompare", Af, Ex);

        Af.m_empty();
        Ex.m_empty();
        m.m_empty();

    }
}
/*
    vector<int> a{1, 2, 3, 4, 5, 6};
    list<int> b{10, 11, 12, 13, 14};
    auto it = b.begin();
    it = b.insert(b.begin(), cbegin(a), cend(a));
    cout << *it;

*/
