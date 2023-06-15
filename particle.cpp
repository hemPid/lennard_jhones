#include "particle.h"
#include "config_const.h"


double smartPow(double x, int16_t n) {
    double res = 1;
    for (uint16_t i = 0; i < abs(n); i++)
    {
        res *= x;
    }
    if (n >= 0)
    {
        return res;
    }
    
    return 1/res;
}


particle::particle()
{
    count = 0;
}

//getters

vect particle::getPosition() { return position; }

vect particle::getVelocity() { return velocity; }

vect particle::getAceleration() { return aceleration; }

vect particle::getAverageVelocity() { return average_velocity / count;}

void particle::setAceleration(vect & force) {
    //sets aceleration a = F / m
    aceleration = force / mass;
}

//position and velocity setters

void particle::setPosition(vect & pos) {
    position = pos;
}

void particle::setPosition(double x, double y, double z) {
    position.x = x;
    position.y = y;
    position.z = z;
}

void particle::setVelocity(vect & vel) {
    velocity = vel;
}

void particle::setVelocity(double x, double y, double z) {
    velocity.x = x;
    velocity.y = y;
    velocity.z = z;
}

void particle::setRandomVelocity(double max_projection) {
    //generates random velocity vector (|v_x,y,z| <= max_projection)
    velocity.x = rand()-RAND_MAX/2;
    velocity.y = rand()-RAND_MAX/2;
    velocity.z = rand()-RAND_MAX/2;
    velocity *= 2*max_projection/RAND_MAX;
}

void particle::updateAverageVelocity() {
    //updates avg velocity coounter
    count += 1;
    average_velocity += velocity;
}

double particle::countDistanceSqr(particle * p) {
    //counts distance to particle p
    return (position - p->getPosition()).sqr_len();
}

void particle::period(vect & r, double cell_size) {
    //perioding vector to count force and potential energy
    //used for periodic boundary conditions
    double buf = 0;
    r.x += cell_size/2;
    r.y += cell_size/2;
    r.z += cell_size/2;
    periodPosition(r.x, buf, cell_size);
    periodPosition(r.y, buf, cell_size);
    periodPosition(r.z, buf, cell_size);
    r.x -= cell_size/2;
    r.y -= cell_size/2;
    r.z -= cell_size/2;
}

void particle::periodPosition(double & axis, double & pre_axis, double cell_size) {
    //periodic boundary conditions
    //used to pereodize movement of particle
    int64_t q = floor(axis/cell_size);
    axis -= q*cell_size;
    pre_axis -= q*cell_size;
}

vect particle::countForce(particle * p, double cell_size) {
    //count force from particle p (lenard-jhones force)
    vect r;
    vect force;
    r = p->position - position;
    period(r, cell_size);
    double len = r.sqr_len();
    force = r * (24*(smartPow(len,-4) - 2*smartPow(len,-7)));
    return force;
}

double particle::countPotentialEnergy(particle * p, double cell_size) {
    //potential energy of interaction with particle p
    vect r;
    r = p->position - position;
    
    period(r, cell_size);
    double dist = r.sqr_len();
    return 4*(smartPow(dist,-6) - smartPow(dist, -3));
}

double particle::countKineticEnergy() {
    //counts kinetic energy
    return (mass*(velocity.sqr_len()))/2;
}

void particle::update(bool initial, double dt, double cell_size) {
    //updating position using aceleration
    vect dr;
    if (initial)
    {
        //if r(t-dt) is not define
        dr = velocity*dt;
    } else {
        //verlet scheme r(t+dt)-r(t)=r(t)-r(t-dt)+a(dt)^2
        dr = position - pre_position + aceleration*dt*dt;
    }
    pre_position = position;
    position += dr;
    dr_from_zero += dr;
    periodPosition(position.x, pre_position.x, cell_size);
    periodPosition(position.y, pre_position.y, cell_size);
    periodPosition(position.z, pre_position.z, cell_size);
    velocity = (position - pre_position)/dt;
}

double particle::countDr2() {
    //counts square offset from initial position
    return dr_from_zero.sqr_len();
}