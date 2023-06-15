#include "vect.h"

//vector methods

vect::vect()
{
    //default constructor
    x = 0;
    y = 0;
    z = 0;
}

vect::vect(double x_p, double y_p, double z_p) {
    x = x_p;
    y = y_p;
    z = z_p;
}

void vect::operator=(vect v) {
    //copying another vector
    x = v.x;
    y = v.y;
    z = v.z;
}


vect vect::operator+(vect a) {
    //addition of vectors
    vect ret_vec(x+a.x, y+a.y, z+a.z);
    return ret_vec;
}
void vect::operator+=(vect a) {
    //add to vector
    x += a.x;
    y += a.y;
    z += a.z;
}


vect vect::operator * (double l) {
    //multiplication on double  
    vect ret_vec(l*x,y*l,z*l);
    return ret_vec;
}

void vect::operator *= (double l) {
    //multyply vector on double
    x*=l;
    y*=l;
    z*=l;
}

double vect::sqr_len() {
    //returns squared length of vector
    return x*x + y*y + z*z;
}



vect vect::operator - (vect a) {
    vect ret(x - a.x, y - a.y, z - a.z);
    return ret;
}
void vect::operator -= (vect a) {
    x -= a.x;
    y -= a.y;
    z -= a.z;
}
vect vect::operator / (double l) {
    vect ret(x/l, y/l, z/l);
    return ret;
}
void vect::operator /= (double l) {
    x /= l;
    y /= l;
    z /= l;
}

