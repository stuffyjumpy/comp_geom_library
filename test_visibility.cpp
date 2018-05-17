#include<bits/stdc++.h>
#include "api/fundamentals.h"
#include "api/visibility.h"
#include <GL/glut.h>

using namespace std;

pair<long double, long double> pt;

vector<point> file_read()
{
    freopen("testcases/point_10000_0.txt","r",stdin);
    vector<point> vertices;

    int n;
    cin>>n;

    for(int i=0;i<n;i++)
    {
        long double val1;
        cin>>val1;

        long double val2;
        cin>>val2;

        point pt;
        pt.set_point(val1*78000.0-39000.0,val2*78000.0-39000.0);

        vertices.push_back(pt);
    }
    cin >> pt.first >> pt.second;
    return vertices;
}

vector<pair<point,point>> edges,diagonals;
vector<point> V;

void init2D(float r, float g, float b)
{
    glClearColor(r,g,b,0.0);
    glMatrixMode (GL_PROJECTION);
    gluOrtho2D (-40000.0, 40000.0, -40000.0, 40000.0);
}

vector<point> ret;
void display()
{
    glClear(GL_COLOR_BUFFER_BIT);
    glColor3f(1.0, 0.0, 0.5);
    GLfloat thickness = 3.0;

    glLineWidth(thickness);
    glBegin(GL_LINES);
    for(int i=0; i<V.size(); i++)
    {
        glVertex2f(V[i].get_point_cartesian().first,
                   V[i].get_point_cartesian().second);
        glVertex2f(V[(i+1)%V.size()].get_point_cartesian().first,
                   V[(i+1)%V.size()].get_point_cartesian().second);
    }
    glEnd();

    glColor3f(0.0, 1.0, 0.0);
    thickness = 1.0;
    glLineWidth(thickness);
    glBegin(GL_LINES);
    for(int i=0; i<ret.size(); i++)
    {
        glVertex2f(ret[i].get_point_cartesian().first,
                   ret[i].get_point_cartesian().second);
        glVertex2f(ret[(i+1)%ret.size()].get_point_cartesian().first,
                   ret[(i+1)%ret.size()].get_point_cartesian().second);
    }
    glEnd();

    glColor3f(1.0, 1.0, 1.0);
    thickness = 7;
    glPointSize(thickness);
    glBegin(GL_POINTS);
        glVertex2f(pt.first, pt.second);
    glEnd();
    glFlush();
}

int main(int argc, char** argv)
{
    //Compile as g++ filename.cpp -std=c++1y -lpthread -lGL -lGLU -lglut
    V = file_read();
    visibility v(V);
    ret = v.find_visibility_polygon(point(pt.first, pt.second));
    for(auto x: ret) cout << x << endl;
    glutInit(&argc, argv);

    glutInitDisplayMode (GLUT_SINGLE | GLUT_RGB);
    glutInitWindowSize (1500, 1500);
    glutInitWindowPosition (1000, 1000);
    glutCreateWindow ("Polygon Triangulation");
    init2D(0.0,0.0,0.0);
    glutDisplayFunc(display);
    glutMainLoop();
}
