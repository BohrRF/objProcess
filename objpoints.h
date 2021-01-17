#ifndef OBJPOINTS_H_INCLUDED
#define OBJPOINTS_H_INCLUDED

#include <string>
#include <vector>
using std::string;
using std::vector;


class objPoint
{
public:
    double m_x;
    double m_y;
    double m_z;
    double m_ang;
    objPoint()
    {
        m_x = 0;
        m_y = 0;
        m_z = 0;
        m_ang = 0;
    }
    objPoint(const objPoint &rval)
    {
        m_x = rval.m_x;
        m_y = rval.m_y;
        m_z = rval.m_z;
        m_ang = rval.m_ang;
    }
};

class polyPoint
{
public:
    objPoint m_P1;
    objPoint m_P2;
    objPoint m_P3;
};

class Cplain
{
public:
    double m_x;
    double m_y;
    double m_z;
    Cplain(const double &x, const double &y, const double &z)
    {
        m_x = x;
        m_y = y;
        m_z = z;
    }
};

class Cline
{
public:
    objPoint P1;
    objPoint P2;
    Cline(const objPoint &P1, const objPoint &P2)
    {
        this->P1 = P1;
        this->P2 = P2;
    }
    Cline(){}
};

class Cslice
{
    vector<objPoint> points;
    vector<polyPoint> polypoints;
public:
    double Angle;
    vector<Cline> slice;
    Cslice(double ang = 0){Angle = ang;}
    void importFromObj(const string &FILENAME);
    void findPointOnSlice();
    void outSlice(const string &FNA);
    void sortByAng(polyPoint &poly);
    objPoint linePlainInter(const Cplain &plain, const Cline &line);
    void CompareOut(const string &FNAME, const Cslice &Aft);
    void m_empty();
};


#endif // OBJPOINTS_H_INCLUDED
