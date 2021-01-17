# objProcess
extract out most layer of given 2D point cloud.

In the processing of extracting out most points, member variable in class CPoint are used partly, which are m_r and m_z as x and y for a single point.

class CSlicePoints and Cmethod are the main container and processor for the extraction.
member variable m_points in CSlicePoints is filled by points as CPoint objects. Becasue vector m_points is a public member of class CSlicePoints, a sigle loop could write points in easily.

Class Cmethod could be used as below

 m.mapCreate(const CSlicePoints& points);
 m.convexHullFind();
 m.method();
 m.outHullFile();
 
last function output points ou the out most layer as file named "hullData(angle).txt" .(angle) states for the value of member variable angle.
