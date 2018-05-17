#ifndef VISIBILITY_H
#define VISIBILITY_H

#include <bits/stdc++.h>
#include "fundamentals.h"

using namespace std;

//! Class which stores a polygon, and the visibility of points can be queried.
class visibility
{
	vector<point> poly;
	point inf;
	stack<int> s;
	vector<point> ret;
	vector< pair<long double, long double> > pts;
	vector<long double> alpha;
	pair<long double, long double> pt;
	int flag;

	string advance(int&, int&);
	string retard(int&, int&);
	string scan(int&, int&);

	public:
	//! Constructor to initialize the class with the vertices of a polygon
    /*!
        \param v a vector of points representing the vertices of polygon in ccw order.
    */
	visibility(const vector<point>&);
	//! Finds visibility polygon with respect to a point.
    /*!
        \param p the point wrt to which the visibility polygon is found.
    */
	vector<point> find_visibility_polygon(point);
};

#include "visibility.cpp"

#endif
