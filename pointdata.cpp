#include "pointdata.h"

#include <iostream>
#include <vector>
#include <list>
#include <string>
#include <fstream>
#include <conio.h>
#include <math.h>
#include <algorithm>

const double PI = 3.14159265358979323846;

const double PARA = 1;
const double MUL = 1;

using std::vector;
using std::list;
using std::string;
using std::ifstream;
using std::ofstream;
using std::ios;
using std::cerr;
using std::cout;
using std::clog;
using std::endl;

vector<CPoint>::size_type CPointDataBase::addPointFromFile(const string& FileName)
{
    CPoint ptemp;
    ifstream f(FileName, ios::in);
    if(!f)
    {
        cerr << "fail";
        exit(0);
    }
    clog << "Reading Magnetic Field Trace Data From " << '<' << FileName << '>' << endl;
    while(f >> ptemp.Bindex >> ptemp.Tindex >> ptemp.m_R >> ptemp.m_w >> ptemp.m_x >> ptemp.m_y >> ptemp.m_z >> ptemp.m_length >> ptemp.m_br >> ptemp.m_bz >> ptemp.m_bt >> ptemp.m_b && !f.eof())
    {
        ptemp.m_R *= MUL;
        ptemp.m_z *= MUL;
        PointData.push_back(ptemp);
    }

    f.close();
    //m_points.begin()->m_show();
    clog << "Loaded " << PointData.size() << " Points From File." << endl << endl;
    return PointData.size();
}

void CPoint::m_show() const
{
    printf("%d (w,r,z) = (%.15f, %.15f, %.15f)\n", Bindex, m_w, m_R, m_z);
}

vector<CPoint>::size_type CSlicePoints::findPointsOnSlice(const CPointDataBase &Database)
{
    CPoint ptemp;

    clog << "Finding Points On " << Angle << " Degree Surface" << endl;
    for(auto a: Database.PointData)
    {
        //ptemp.m_show();
        ptemp = a;
        while(ptemp.m_w>=360)
            ptemp.m_w -= 360;
        while(ptemp.m_w < 0)
            ptemp.m_w += 360;
        if(ptemp.m_w - Angle > 1E-7 || ptemp.m_w - Angle < (-1E-7))
            continue;
        else
            m_points.push_back(ptemp);
    }
    //m_points.begin()->m_show();
    clog << m_points.size() << " Points Found" << endl << endl;
    return m_points.size();
}




void CSlicePoints::showAll()
{
    cout << m_points.size() << "points in this slice. Show all of them? (return for continue)" << endl;
    if(getch() == '\r')
    {
        for(auto a : m_points)
        {
            a.m_show();
        }
    }
}

void CSlicePoints::addLar()
{
    vector<CPoint> temp;
    CPoint Ptemp;
    for(auto a: m_points)
    {
        for(int i = 0; i < 16; i++)
        {
            Ptemp.m_R = a.m_R + MUL * 0.3 * cos(2 * PI * i / 16);
            Ptemp.m_z = a.m_z + MUL * 0.3 * sin(2 * PI * i / 16);
            temp.emplace_back(Ptemp);
        }
    }
    m_points.clear();
    for(auto a: temp)
    {
        m_points.push_back(a);
    }
}



vector<XYPoint>::size_type CMethod::mapCreate(const CSlicePoints &slicePoints)
{
    int sIndex = 0;
    Angle = slicePoints.Angle;
    for(auto a: slicePoints.m_points)
    {
        uP.push_back(a);
        uP.back().index = sIndex++;
    }
    return pnum = uP.size();
}

double CMethod::cross(const XYPoint &P1, const XYPoint &P2, const XYPoint &P3)
{
    return (P2.m_r-P1.m_r)*(P3.m_z-P1.m_z)-(P3.m_r-P1.m_r)*(P2.m_z-P1.m_z);
}

double CMethod::angCal(const XYPoint &P1, const XYPoint &P2, const XYPoint &P3)
{
    double x1, x2, y1, y2, c, s;
    x1 = P1.m_r - P2.m_r;
    y1 = P1.m_z - P2.m_z;
    x2 = P3.m_r - P2.m_r;
    y2 = P3.m_z - P2.m_z;
    c = acos(((x1 * x2 + y1 * y2) / (sqrt(x1 * x1 + y1 * y1) * sqrt(x2 * x2 + y2 * y2))>=1)?1:(x1 * x2 + y1 * y2) / (sqrt(x1 * x1 + y1 * y1) * sqrt(x2 * x2 + y2 * y2)));
    s = asin(((x1 * y2 - y1 * x2) / (sqrt(x1 * x1 + y1 * y1) * sqrt(x2 * x2 + y2 * y2))>=1)?1:(x1 * y2 - y1 * x2) / (sqrt(x1 * x1 + y1 * y1) * sqrt(x2 * x2 + y2 * y2)));
    /*
    printf("p1(%f, %f)\n", P1.m_r, P1.m_z);
    printf("p2(%f, %f)\n", P2.m_r, P2.m_z);
    printf("p3(%f, %f)\n", P3.m_r, P3.m_z);
    printf("x1: %f, y1: %f, x2: %f, y2: %f\n", x1, y1, x2, y2);
    printf("x*y: %f\n", (x1 * x2 + y1 * y2) / (sqrt(x1 * x1 + y1 * y1) * sqrt(x2 * x2 + y2 * y2)));

    printf("c: %f, s: %f\n\n", c, s);
    */
    if(s >= 0)
        return c;
    else
        return 2 * PI - c;
}

double CMethod::disCal(const XYPoint &P1, const XYPoint &P2, const XYPoint &P3)
{
    double x1, x2, y1, y2;
    x1 = P3.m_r - P1.m_r;
    y1 = P3.m_z - P1.m_z;
    x2 = P2.m_r - P1.m_r;
    y2 = P2.m_z - P1.m_z;
    /*
    printf("p1(%f, %f)\n", P1.m_r, P1.m_z);
    printf("p2(%f, %f)\n", P2.m_r, P2.m_z);
    printf("p3(%f, %f)\n", P3.m_r, P3.m_z);
    printf("distance: %f\n", sqrt((x1 - cx) * (x1 - cx) + (y1 - cy) * (y1 - cy)));
    */
    return sqrt(x1 * x1 + y1 * y1 - (x1 * x2 + y1 * y2) * (x1 * x2 + y1 * y2) / (x2 * x2 + y2 * y2));
}


vector<XYPoint>::size_type CMethod::convexHullFind()
{
    vector<XYPoint>::iterator it = uP.begin();
    vector<XYPoint>::iterator pit;
    vector<XYPoint> temp;
    vector<XYPoint>::iterator miny;

    clog << "Finding Convex Hull." << endl;

    for(auto a: uP)
    {
        temp.push_back(a);
    }

    miny = temp.begin();
    for(auto it = temp.begin(); it != temp.end(); it++) //find min y
    {
        if((*it).m_z < (*miny).m_z)
            miny = it;
    }
    //cout << "miny->m_r:" << miny->m_r << " miny->m_z:"<< miny->m_z << endl;

    for(auto it = temp.begin(); it != temp.end(); it++) //plot recalculate for sorting
    {
        double tempr = it->m_r - (*miny).m_r;
        double tempz = it->m_z - (*miny).m_z;
        it->m_w = atan2(tempz, tempr);
        //cout << tempr << "," << tempz << "," << a.m_w << endl;
    }
    (*miny).m_w = -5;

    sort(temp.begin(), temp.end(), [=](const XYPoint &para1, const XYPoint &para2){return para1.m_w <= para2.m_w;});

    P.clear();
    P.push_back(*(temp.begin()));
    P.push_back(*(temp.begin() + 1));
    for(auto itt = temp.begin() + 2; itt != temp.end(); itt++)
    {
        pit = P.end();
        while(cross(*(pit - 2), *(pit - 1), *itt) < 0)
        {
            pit = P.erase(pit - 1);
        }
        P.push_back(*itt);
    }

    temp.clear();
    it = P.begin();
    while(it != P.end())
    {
        temp.push_back(*it);
        ++it;
    }

    sort(temp.begin(), temp.end(), [=](const XYPoint para1, const XYPoint para2){return para1.index < para2.index;});
    pit = temp.begin();
    auto puP = uP.begin();
    while(puP != uP.end())
    {
        //cout << "uP:" << (*puP).index << "P:" << (*pit).index << endl;
        if(*puP == *pit)
        {
            puP = uP.erase(puP) - 1;
            if(++pit == temp.end())break;
        }
        puP++;
    }

    clog << "Found " << P.size() << " Points on Convex Hull." << endl << endl;
    return P.size();
}


void CMethod::showHull() const
{
    for(auto a: P)
    {
        printf("%d (w,r,z) = (%.15f, %.15f, %.15f)\n", a.index, a.m_w, a.m_r, a.m_z);
    }
}

void CSlicePoints::outSliceFile() const
{
    ofstream fout("sliceData" + std::to_string(Angle)+ ".txt", ios::out);
    if(!fout)
    {
        cerr << "fail.";
    }

    clog << "Outputting " << '<' << "sliceData" << std::to_string(Angle) << ".txt" << '>' << endl;
    for(auto a: m_points)
    {
        fout << a.m_R << '\t' << a.m_z << endl;
    }
    fout.close();

    clog << "Output finished."<< endl << endl;
}


void CSlicePoints::compareOut(const string &FNAME, const Cslice &Aft)
{
    vector<CPoint>::const_iterator a = m_points.begin();
    vector<Cline>::const_iterator b = Aft.slice.begin();
    ofstream fout(FNAME + "-" + std::to_string(Angle) + "degree.txt", ios::out);
    if(fout.fail())
    {
        cerr << "fail to open";
    }

    clog << "Outputting " << FNAME << "-" << std::to_string(Angle) << "degree.txt" << endl;
    while(1)
    {
        if(a != m_points.end())
        {
            fout << a->m_R << '\t' << a->m_z << '\t';
            a++;
            if(b != Aft.slice.end())
            {
                // Ori(a) and Aft(b) neither reached end
                fout << sqrt(b->P1.m_x * b->P1.m_x + b->P1.m_y * b->P1.m_y) << '\t' << b->P1.m_z << endl;
                if(a != m_points.end())
                {
                    // Ori(a) and Aft(b) neither reached end
                    fout << a->m_R << '\t' << a->m_z << '\t'<< sqrt(b->P2.m_x * b->P2.m_x + b->P2.m_y * b->P2.m_y) << '\t' << b->P2.m_z << endl;
                    a++;b++;
                }
                else
                {
                    // Ori(a) reached end but Aft(b) doesn't
                    fout << '\t' << '\t' << sqrt(b->P2.m_x * b->P2.m_x + b->P2.m_y * b->P2.m_y) << '\t' << b->P2.m_z << endl;
                    b++;
                }
            }
            else
            {
                // Aft(b) reached end but Ori(a) doesn't
                fout << '\t' << '\t' << endl;
                if(a != m_points.end())
                {
                    fout << a->m_R << '\t' << a->m_z << '\t' << '\t' << endl;
                    a++;
                }
                else
                    //both reached end, finish
                    break;
            }
        }
        else
        {
            fout << '\t' << '\t';
            if(b != Aft.slice.end())
            {
                //Ori(a) reached end but Aft(b) doesn't
                fout << sqrt(b->P1.m_x * b->P1.m_x + b->P1.m_y * b->P1.m_y) << '\t' << b->P1.m_z << endl;
                fout << '\t' << '\t' << sqrt(b->P2.m_x * b->P2.m_x + b->P2.m_y * b->P2.m_y) << '\t' << b->P2.m_z << endl;
                b++;
            }
            else
                //both reached end, finish
                break;
        }
    }
    clog << "Output finished." << endl << endl;
}

void CMethod::outUp() const
{
    static int s = 0;
    s++;
    ofstream fout("uP" + std::to_string(s)+ ".txt", ios::out);
    if(!fout)
    {
        cerr << "fail.";
    }
    for(auto a: uP)
    {
        fout << a.m_r << '\t' << a.m_z << endl;
    }
    fout.close();
}


void CMethod::outHullFile() const
{
    ofstream fout("hullData" + std::to_string(Angle) + ".txt", ios::out);
    if(!fout)
    {
        cerr << "fail.";
    }
    for(auto a: P)
    {
        fout << a.m_r << '\t' << a.m_z << endl;
    }
    fout.close();
}

double CMethod::areaCal() const
{
    double area = 0;
    for(auto it = P.begin(); (it + 1)!= P.end(); it++)
    {
        area += it->m_r * (it + 1)->m_z - (it + 1)->m_r * it -> m_z;
    }
    return PARA * (area + P.back().m_r * P.front().m_z - P.front().m_r * P.back().m_z);
}

double CMethod::lenCal() const
{
    double len = 0;
    for(auto it = P.begin(); it != P.end(); it++)
    {
        len += sqrt(((it + 1)->m_z - it->m_z) * ((it + 1)->m_z - it->m_z) + ((it + 1)->m_r - it->m_r) * ((it + 1)->m_r - it->m_r));
    }
    return len + sqrt((P.back().m_z - P.front().m_z) * (P.back().m_z - P.front().m_z) + (P.back().m_r - P.front().m_r) * (P.back().m_r - P.front().m_r));
}

void CMethod::paraCal()
{
    double aveL = 0, amin = 0, amax = 0;
    aveL = sqrt(areaCal() / pnum);
    mju = lenCal() / (aveL * pnum);

    for(auto it = P.begin(); it != P.end(); it++)
    {
        delt += pow(sqrt(((it + 1)->m_z - it->m_z) * ((it + 1)->m_z - it->m_z) + ((it + 1)->m_r - it->m_r) * ((it + 1)->m_r - it->m_r)) / aveL - 1, 2);
    }
    delt += pow(sqrt((P.back().m_z - P.front().m_z) * (P.back().m_z - P.front().m_z) + (P.back().m_r - P.front().m_r) * (P.back().m_r - P.front().m_r)) / aveL - 1, 2);
    delt /= P.size();
    amin = ((PI - 2 * asin(mju) - 0.5 * PI / delt) > 0)?(PI - 2 * asin(mju) - 0.5 * PI / delt): 0;
    amax = ((PI - 2 * asin(mju) - 0.5 * PI / delt) < PI)?(PI - 2 * asin(mju) + 0.5 * PI / delt): PI;
    alpha = 0.5 * (amax - amin) + amin;
}

double CMethod::getAlpha() const
{
    return alpha;
}

void CMethod::method()
{
    vector<vector<XYPoint>::iterator> candi;
    vector<double> dis;
    list<XYPoint> temp;
    double ang;
    clog << "Calculating Parameter" << endl;
    paraCal();
    //alpha = PI / 2;
    clog << "Parameter Calculate Finished" << alpha << endl << endl;

    clog << "Calculating Outmost Edge of Lamar Points" << endl;
    for(auto a: P)
    {
        temp.push_back(a);
    }
    temp.push_back(P.front());

    auto it = temp.begin();
    while(it != temp.end())
    {
        for(auto ita = uP.begin(); ita != uP.end(); ita++)
        {
            auto next = it;
            next++;
            ang = angCal(*it, *ita, *next);
            if(ang <= (PI - alpha / 2) || ang > PI)
                continue;
            else
            {
                //cout << ang << endl;
                candi.push_back(ita);
                dis.push_back(disCal(*it, *next, *ita));
            }
        }
        if(dis.size() == 0)
        {
            //cout << "go to next point on hull" << endl;
            it++;
        }
        else
        {
            auto pmin = dis.begin();
            for(auto a = dis.begin(); a != dis.end(); a++) //find min
            {
                //cout << *a << endl;
                if(*a < *pmin)
                    pmin = a;
            }
            //cout << endl;
            auto next = it;
            next++;
            temp.insert(next, *(candi[pmin - dis.begin()]) ); // add point
            uP.erase(candi[pmin - dis.begin()]);          //delete point
            candi.clear();
            dis.clear();
            //it = P.begin();
        }
    }
    P.clear();
    for(auto a: temp)
    {
        P.push_back(a);
    }
    clog << "Calculation finished" << endl << endl;
}

void CMethod::compareOut(const string& FNAME, const Cslice &Aft, const CSlicePoints &La)
{
    vector<CPoint>::const_iterator a = La.m_points.begin(); //lamer points
    vector<XYPoint>::const_iterator b = P.begin();          //lamer edge
    vector<Cline>::const_iterator c = Aft.slice.begin();    //MC after edited

    ofstream fout(FNAME + "-" + std::to_string(Angle) + "degree.txt", ios::out);
    if(fout.fail())
    {
        cerr << "failed to open";
    }

    clog << "Outputting " << '<' << FNAME + "-" + std::to_string(Angle) + "degree.txt" << '>' << endl;

    fout << "Lamer Points" << '\t' << '\t' << "Lamer Edge" << '\t' << '\t' << "MC After Edited" << endl;
    while(a != La.m_points.end())
    {
        fout << a->m_R << '\t' << a->m_z << '\t';
        a++;

        if(b != P.end())
        {
            fout << b->m_r << '\t' << b->m_z << '\t';
            b++;

            if(c != Aft.slice.end())
            {
                fout << sqrt(c->P1.m_x * c->P1.m_x + c->P1.m_y * c->P1.m_y) << '\t' << c->P1.m_z << endl;
                if(b != P.end())
                {
                    fout << a->m_R << '\t' << a->m_z << '\t' << b->m_r << '\t' << b->m_z << '\t' << sqrt(c->P2.m_x * c->P2.m_x + c->P2.m_y * c->P2.m_y) << '\t' << c->P2.m_z << endl;
                    b++;
                }
                else
                {
                    fout << a->m_R << '\t' << a->m_z << '\t' << '\t' << '\t' << sqrt(c->P2.m_x * c->P2.m_x + c->P2.m_y * c->P2.m_y) << '\t' << c->P2.m_z << endl;
                }
                a++;
                c++;
            }
            else
                fout << endl;
        }
        else
        {
            if(c != Aft.slice.end())
            {
                fout << '\t' << '\t' << sqrt(c->P1.m_x * c->P1.m_x + c->P1.m_y * c->P1.m_y) << '\t' << c->P1.m_z << endl;
                fout << a->m_R << '\t' << a->m_z << '\t' << '\t' << '\t' << sqrt(c->P2.m_x * c->P2.m_x + c->P2.m_y * c->P2.m_y) << '\t' << c->P2.m_z << endl;
                a++;
                c++;
            }
            else
                fout << endl;
        }
    }
    fout.close();
    clog << "Output Finished" << endl << endl;
}

void CSlicePoints::m_empty()
{
    m_points.clear();
}

void CMethod::m_empty()
{
    P.clear();
    uP.clear();
}




