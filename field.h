#pragma once
#include "config_const.h"
#include "particle.h"
#include "vector"
#include <cmath>
#include <cstdint>
#include <fstream>

class field
{
private:
    bool free_run_active;
    double cell_size;
    double dt;
    int64_t n;
    particle * parts;
    double potential_enery;
    double kinetic_energy;
    double full_energy;
    void count_interactions();
    void updatePositions();
public:
    field(int64_t n_p, double c_size, double temp);
    ~field();
    void printInfo();
    void makeTicks(int16_t ticks_count);
    void tick();
    double countAverageDr2();
};
