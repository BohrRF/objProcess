#include "objpoints.h"
#include <fstream>
#include <iostream>
#include <math.h>

const double PI = 3.1415926536;
const double center_x = 120;
const double center_y = 120;
const double center_z = 35;
const double mul_x = 0.2;//0.2
const double mul_y = 0.2;//0.2
const double mul_z = 7.55/35;//7.55/35

using std::ios;
using std::ifstream;
using std::ofstream;
using std::cout;
using std::endl;
using std::cin;
using std::cerr;
using std::clog;

void Cslice::sortByAng(polyPoint &poly)
{
    objPoint temp;
    if(poly.m_P1.m_ang > poly.m_P2.m_ang)
    {
        temp = poly.m_P1;
        poly.m_P1 = poly.m_P2;
        poly.m_P2 = temp;
    }
    if(poly.m_P2.m_ang > poly.m_P3.m_ang)
    {
        temp = poly.m_P2;
        poly.m_P2 = poly.m_P3;
        poly.m_P3 = temp;
    }
}

void Cslice::importFromObj(const string &FILENAME)
{
    objPoint temp;
    polyPoint ptemp;
    vector<objPoint>::size_type p1, p2, p3;
    char c;
    ifstream fin(FILENAME, ios::in);
    if(fin.fail())
    {
        cerr << "failed to open.";
        exit(0);
    }
    clog << "Reading Point Information From " << '<' << FILENAME << '>' << endl;
    while(fin >> c)
    {
        if(c == 'v')
        {
            fin >> temp.m_x >> temp.m_y >> temp.m_z;
            temp.m_x -= center_x;
            temp.m_y = -temp.m_y + center_y;
            temp.m_z -= center_z;
            temp.m_x *= mul_x;
            temp.m_y *= mul_y;
            temp.m_z *= mul_z;
            //cout << "y=" << temp.m_y << " x=" << temp.m_x;
            temp.m_ang = atan2(temp.m_y, temp.m_x);

            if(temp.m_ang < 0)
            {
                temp.m_ang += 2 * PI;
                temp.m_ang = (temp.m_ang >= 2 * PI)? 0: temp.m_ang;
            }
            /*
            if(temp.m_ang > PI)
            {
                cout << " tan(y/x)=" << temp.m_ang << endl;
                cout << temp.m_x << ',' << temp.m_y << ',' << temp.m_z << ' ' << temp.m_ang << endl << endl;
            }
            */
            points.push_back(temp);
        }
        else if(c == 'f')
        {
            fin >> p1 >> p2 >> p3;
            ptemp.m_P1 = points[p1 - 1];
            ptemp.m_P2 = points[p2 - 1];
            ptemp.m_P3 = points[p3 - 1];
            polypoints.push_back(ptemp);
        }
    }
    clog << points.size() << " Points and " << polypoints.size() << " Polygons were Loaded" << endl << endl;
}

objPoint Cslice::linePlainInter(const Cplain &plain, const Cline &line)
{
    objPoint temp;
    double vp1, vp2, vp3, v1, v2, v3, m1, m2, m3, t, vpt;
    vp1 = plain.m_x;
    vp2 = plain.m_y;
    vp3 = plain.m_z;
    v1 = line.P1.m_x - line.P2.m_x;
    v2 = line.P1.m_y - line.P2.m_y;
    v3 = line.P1.m_z - line.P2.m_z;
    m1 = line.P2.m_x;
    m2 = line.P2.m_y;
    m3 = line.P2.m_z;
    vpt = v1 * vp1 + v2 * vp2 + v3 * vp3;

    if (vpt == 0)
    {
        clog << "parallel" << endl;
        temp.m_x = m1 + v1 / 2;
        temp.m_y = m2 + v2 / 2;
        temp.m_z = m3 + v3 / 2;
        temp.m_ang = atan2(temp.m_y, temp.m_x);
        if(temp.m_ang < 0)
        {
            temp.m_ang += 2 * PI;
            temp.m_ang = (temp.m_ang >= 2 * PI)? 0: temp.m_ang;
        }
    }
    else
    {
        t = -(m1 * vp1 + m2 * vp2 + m3 * vp3) / vpt;
        if(t > 1 || t < 0)   //avoid parallel number bug
        {
            temp.m_x = m1 + v1 / 2;
            temp.m_y = m2 + v2 / 2;
            temp.m_z = m3 + v3 / 2;
            temp.m_ang = atan2(temp.m_y, temp.m_x);
            if(temp.m_ang < 0)
            {
                temp.m_ang += 2 * PI;
                temp.m_ang = (temp.m_ang >= 2 * PI)? 0: temp.m_ang;
            }
            return temp;
        }

        temp.m_x = m1 + v1 * t;
        temp.m_y = m2 + v2 * t;
        temp.m_z = m3 + v3 * t;
        temp.m_ang = atan2(temp.m_y, temp.m_x);
        if(temp.m_ang < 0)
        {
            temp.m_ang += 2 * PI;
            temp.m_ang = (temp.m_ang >= 2 * PI)? 0: temp.m_ang;
        }
    }
    return temp;
}

void Cslice::findPointOnSlice()
{
    Cline temp, line;
    double ArcAngle = Angle * PI / 180;
    Cplain plain(cos(ArcAngle + PI / 2), sin(ArcAngle + PI / 2), 0);

    clog << "Calculating interact between " << Angle << " Degree Surface with Polygons" << endl;
    for(vector<polyPoint>::iterator it = polypoints.begin(); it != polypoints.end(); it++)
    {
        sortByAng(*it);
        //clog << " p1:" << it->m_P1.m_ang << " p2:" << it->m_P2.m_ang << " p3:" << it->m_P3.m_ang << endl;
        if(it->m_P1.m_ang < PI / 4)
        {
            if(it->m_P2.m_ang > PI * 7 / 4)
            {
                //clog << " p1:" << it->m_P1.m_ang << " p2:" << it->m_P2.m_ang << " p3:" << it->m_P3.m_ang << endl;
                it->m_P2.m_ang -= 2 * PI;
                it->m_P3.m_ang -= 2 * PI;
                //clog << " p1:" << it->m_P1.m_ang << " p2:" << it->m_P2.m_ang << " p3:" << it->m_P3.m_ang << endl << endl;

            }
            else if(it->m_P3.m_ang > PI * 7 / 4)
            {
                //clog << " p1:" << it->m_P1.m_ang << " p2:" << it->m_P2.m_ang << " p3:" << it->m_P3.m_ang << endl;
                it->m_P3.m_ang -= 2 * PI;
                //clog << " p1:" << it->m_P1.m_ang << " p2:" << it->m_P2.m_ang << " p3:" << it->m_P3.m_ang << endl << endl;
            }
            sortByAng(*it);
        }

        if((ArcAngle >= it->m_P1.m_ang) && (ArcAngle <= it->m_P3.m_ang))
        {
            if(ArcAngle >= it->m_P2.m_ang)
            {
                //clog << "ArcAngle >= it->m_P2.m_ang" << "Arc:" << ArcAngle << " p1:" << it->m_P1.m_ang << " p2:" << it->m_P2.m_ang << " p3:" << it->m_P3.m_ang << endl;
                line.P1 = it->m_P1;
                line.P2 = it->m_P3;
                temp.P1 = linePlainInter(plain, line);
                line.P1 = it->m_P2;
                line.P2 = it->m_P3;
                temp.P2 = linePlainInter(plain, line);
                if(temp.P1.m_z > 100 || temp.P2.m_z > 100 || abs(temp.P1.m_x) < 1 || abs(temp.P2.m_x) < 1)
                {
                    cout << it->m_P1.m_x << ',' << it->m_P1.m_y << ',' << it->m_P1.m_z << endl;
                    cout << it->m_P3.m_x << ',' << it->m_P3.m_y << ',' << it->m_P3.m_z << endl;
                    cout << temp.P1.m_x << ',' << temp.P1.m_y << ',' << temp.P1.m_z << ' ' << temp.P1.m_ang << endl << endl;
                    cout << it->m_P2.m_x << ',' << it->m_P2.m_y << ',' << it->m_P2.m_z << endl;
                    cout << it->m_P3.m_x << ',' << it->m_P3.m_y << ',' << it->m_P3.m_z << endl;
                    cout << temp.P2.m_x << ',' << temp.P2.m_y << ',' << temp.P2.m_z << ' ' << temp.P2.m_ang << endl << endl;
                }
            }
            else
            {
                //clog << "ArcAngle >= it->m_P2.m_ang" << "Arc:" << ArcAngle << " p1:" << it->m_P1.m_ang << " p2:" << it->m_P2.m_ang << " p3:" << it->m_P3.m_ang << endl;
                line.P1 = it->m_P1;
                line.P2 = it->m_P3;
                temp.P1 = linePlainInter(plain, line);
                line.P1 = it->m_P1;
                line.P2 = it->m_P2;
                temp.P2 = linePlainInter(plain, line);
                if(temp.P1.m_z > 100 || temp.P2.m_z > 100 || abs(temp.P1.m_x) < 1 || abs(temp.P2.m_x) < 1)
                {
                    cout << it->m_P1.m_x << ',' << it->m_P1.m_y << ',' << it->m_P1.m_z << endl;
                    cout << it->m_P3.m_x << ',' << it->m_P3.m_y << ',' << it->m_P3.m_z << endl;
                    cout << temp.P1.m_x << ',' << temp.P1.m_y << ',' << temp.P1.m_z << ' ' << temp.P1.m_ang << endl << endl;
                    cout << it->m_P1.m_x << ',' << it->m_P1.m_y << ',' << it->m_P1.m_z << endl;
                    cout << it->m_P2.m_x << ',' << it->m_P2.m_y << ',' << it->m_P2.m_z << endl;
                    cout << temp.P2.m_x << ',' << temp.P2.m_y << ',' << temp.P2.m_z << ' ' << temp.P2.m_ang << endl << endl;
                }
            }
            slice.push_back(temp);
        }
    }
    clog << slice.size() << " Points on " << Angle << " Degree Slice were Found" << endl << endl;
}

void Cslice::outSlice(const string &FNAME)
{
    ofstream fout(FNAME + "-" + std::to_string(Angle) + "degree.txt", ios::out);
    if(fout.fail())
    {
        cerr << "fail to open";
    }
    double p1x, p1y, p2x, p2y;
    for(auto a: slice)
    {
        p1x = sqrt(a.P1.m_x * a.P1.m_x + a.P1.m_y * a.P1.m_y);
        p2x = sqrt(a.P2.m_x * a.P2.m_x + a.P2.m_y * a.P2.m_y);
        p1y = a.P1.m_z;
        p2y = a.P2.m_z;
        fout << p1x << '\t' << p1y << endl;
        fout << p2x << '\t' << p2y << endl;
    }
    fout.close();
}

void Cslice::CompareOut(const string &FNAME, const Cslice &Aft)
{
    vector<Cline>::const_iterator a = slice.begin();
    vector<Cline>::const_iterator b = Aft.slice.begin();
    ofstream fout(FNAME + "-" + std::to_string(Angle) + "degree.txt", ios::out);
    if(fout.fail())
    {
        cerr << "fail to open";
    }
    while(1)
    {
        if(a != slice.end())
        {
            fout <<  sqrt(a->P1.m_x * a->P1.m_x + a->P1.m_y * a->P1.m_y) << '\t' << a->P1.m_z << '\t';
            if(b != Aft.slice.end())
            {
                // Ori(a) and Aft(b) neither reached end
                fout << sqrt(b->P1.m_x * b->P1.m_x + b->P1.m_y * b->P1.m_y) << '\t' << b->P1.m_z << endl;
                fout << sqrt(a->P2.m_x * a->P2.m_x + a->P2.m_y * a->P2.m_y) << '\t' << a->P2.m_z << '\t'<< sqrt(b->P2.m_x * b->P2.m_x + b->P2.m_y * b->P2.m_y) << '\t' << b->P2.m_z << endl;
                a++;b++;
            }
            else
            {
                // Aft(b) reached end but Ori(a) doesn't
                fout << '\t' << '\t' << endl;
                fout << sqrt(a->P2.m_x * a->P2.m_x + a->P2.m_y * a->P2.m_y) << '\t' << a->P2.m_z << '\t' << '\t' << endl;
                a++;
            }
        }
        else
        {
            //Ori(a) reached end but Aft(b) doesn't
            fout << '\t' << '\t' << endl;
            if(b != Aft.slice.end())
            {
                fout << sqrt(b->P1.m_x * b->P1.m_x + b->P1.m_y * b->P1.m_y) << '\t' << b->P1.m_z << endl;
                fout << '\t' << '\t' << sqrt(b->P2.m_x * b->P2.m_x + b->P2.m_y * b->P2.m_y) << '\t' << b->P2.m_z << endl;
                b++;
            }
            else
                //both reached end, finish
                break;
        }
    }
}


void Cslice::m_empty()
{
    slice.clear();
}





