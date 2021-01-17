#ifndef POINTDATA_H_INCLUDED
#define POINTDATA_H_INCLUDED
#include "objpoints.h"
#include <vector>
#include <string>
#include <math.h>

using std::vector;
using std::string;

class CPoint
{
public:
    int Bindex;
    int Tindex;
    double m_R;
    double m_w;
    double m_x;
    double m_y;
    double m_z;
    double m_length;
    double m_br;
    double m_bz;
    double m_bt;
    double m_b;
    CPoint()
    {
        Bindex = 0;
        Tindex = 0;
        m_R = 0;
        m_w = 0;
        m_x = 0;
        m_y = 0;
        m_z = 0;
        m_length = 0;
        m_br = 0;
        m_bz = 0;
        m_bt = 0;
        m_b = 0;
    }
    void m_show() const;
    void operator =(const CPoint &rval)
    {
        this->Bindex = rval.Bindex;
        this->Tindex = rval.Tindex;
        this->m_R = rval.m_R;
        this->m_w = rval.m_w;
        this->m_x = rval.m_x;
        this->m_y = rval.m_y;
        this->m_z = rval.m_z;
        this->m_length = rval.m_length;
        this->m_br = rval.m_br;
        this->m_bz = rval.m_bz;
        this->m_bt = rval.m_bt;
        this->m_b = rval.m_b;
    }
};

class CPointDataBase
{
public:
    vector<CPoint> PointData;
    vector<CPoint>::size_type addPointFromFile(const string &FileName);
};


class CSlicePoints
{
public:
    double Angle;
    std::vector<CPoint> m_points;
    vector<CPoint>::size_type findPointsOnSlice(const CPointDataBase &Database);
    void addLar();
    void outSliceFile() const;
    void showAll();
    void compareOut(const string &FNAME, const Cslice &Aft);
    void m_empty();
    CSlicePoints(double deg = 0)
    {
        Angle = deg;
    }
};

class XYPoint
{
public:
    int index;
    double m_r;
    double m_z;
    double m_w;
    XYPoint(){}
    XYPoint(const CPoint &point)
    {
        m_r = point.m_R;
        this->m_z = point.m_z;
        m_w = atan2(this->m_z, m_r);
    }
    bool operator < (XYPoint point)
    {
        return this->m_w < point.m_w;
    }
    bool operator == (XYPoint point)
    {
        return this->index == point.index;
    }
};

class CMethod
{
    double Angle;
    int pnum;
    double mju;
    double delt;
    double alpha;
    double areaCal() const;
    double lenCal() const;
    void paraCal();
    double cross(const XYPoint &P1, const XYPoint &P2, const XYPoint &P3);
    double angCal(const XYPoint &P1, const XYPoint &P2, const XYPoint &P3);
    double disCal(const XYPoint &P1, const XYPoint &P2, const XYPoint &P3);
    vector<XYPoint> uP;
    vector<XYPoint> P;

public:
    vector<XYPoint>::size_type mapCreate(const CSlicePoints &slicePoints);
    vector<XYPoint>::size_type convexHullFind();

    void method();
    void showHull() const;
    double getAlpha() const;
    void outUp() const;
    void outHullFile() const;
    void compareOut(const string &FNAME, const Cslice &Aft, const CSlicePoints &La);
    void m_empty();
};



#endif // POINTDATA_H_INCLUDED










