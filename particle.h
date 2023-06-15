#pragma once
#include "vect.h"
#include <cmath>

//class of a single particle

class particle
{
private:
    //r(t)
    vect position;
    //r(t-dt)
    vect pre_position;
    vect velocity;
    vect aceleration;
    //counters for average velocity
    vect average_velocity;
    int16_t count;
    //offset from initial position
    vect dr_from_zero;
    //periodic boundary conditions realization
    void period(vect & r, double cell_size);
    void periodPosition(double & axis,  double & pre_axis, double cell_size);
public:
    particle();
    /*basic part*/
    //getters
    vect getPosition();
    vect getVelocity();
    vect getAceleration();
    vect getAverageVelocity();
    //setters
    void setAceleration(vect & force);
    void setPosition(vect & pos);
    void setPosition(double x, double y, double z);
    void setVelocity(vect & vel);
    void setVelocity(double x, double y, double z);
    void setRandomVelocity(double max_projection);
    void updateAverageVelocity();
    double countDistanceSqr(particle * p);
    vect countForce(particle * p, double cell_size);
    double countPotentialEnergy(particle * p, double cell_size);
    double countKineticEnergy();
    void update(bool initial, double dt, double cell_size);
    double countDr2();
};

double smartPow(double x, int16_t n);
