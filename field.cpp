#include "field.h"

field::field(int64_t n_p, double c_size, double temp)
{
    //variables defenition
    n = n_p;
    cell_size = c_size;
    dt = dtime;
    double max_velocity = sqrt(3*temp);
    parts = new particle[n];
    free_run_active = false;
    //crystal structure
    int64_t side_count = ceil(pow(n, 1.0/3));
    double crystal_step = cell_size/side_count;
    vect full_momentum;
    int64_t i = 0;
    int64_t j = 0;
    int64_t k = 0;
    for (size_t p = 0; p < n; ++p)
    {
        parts[p].setPosition(
            i*crystal_step,
            j*crystal_step,
            k*crystal_step
        );
        parts[p].setRandomVelocity(max_velocity);
        parts[p].update(true, dt, cell_size);
        full_momentum += parts[p].getVelocity() * mass;
        ++i;
        if (i >= side_count)
        {
            i = 0;
            ++j;
        }
        if (j >= side_count)
        {
            j = 0;
            ++k;
        }
    }
    //making full_momentum = 0
    full_momentum /= n;
    for (size_t p = 0; p < n; ++p)
    {
        vect d;
        d = parts[p].getVelocity() - full_momentum;
        parts[p].setVelocity(d);
    }
}

field::~field()
{
    delete [] parts;
}

void field::printInfo() {
    for (size_t p = 0; p < n; ++p)
    {
        vect pos = parts[p].getPosition();
        vect vel = parts[p].getVelocity();
        std::cout << p << ":\n";
        std::cout << "Posintion: x=" << pos.x << " | " <<
        "y=" << pos.y << " | " << "z=" << pos.z << "\n";
        std::cout << "Velocity: x=" << vel.x << " | " <<
        "y=" << vel.y << " | " << "z=" << vel.z << "\n";
    }
}

void field::makeTicks(int16_t ticks_count) {
    std::ofstream energy_log, positions, diffusion;
    energy_log.open(energy_path);
    positions.open(pos_path);
    diffusion.open(diffusion_path);
    for (size_t i = 0; i < ticks_count; i++)
    {
        if (!(i % 100))
        {
            std::cout << "Tick: " << i << "/" << ticks_count << "\n";
        }
        tick();
        energy_log << "Tick: " << i << ", energy: " <<
        full_energy << ", penergy: " << potential_enery <<
        ", kenergy: " << kinetic_energy << "\n";
        diffusion << "Tick: " << i << ", diff: " << countAverageDr2() << "\n";
        positions << n << "\n";
        positions << "Frame " << i << "\n";
        for (int64_t j = 0; j < n; ++j)
        {
            positions << 1 << " " << parts[j].getPosition().x << " " <<
            parts[j].getPosition().y << " " <<
            parts[j].getPosition().z << "\n";
            if (i >= 0.9*ticks_count)
            {
                parts[j].updateAverageVelocity();     
            }
            
        }
    }
    energy_log.close();
    positions.close();
    diffusion.close();
    std::ofstream vels_out;
    vels_out.open(velocity_path);
    vect vel;
    vels_out << "n=" << n << "\n";
    for (size_t i = 0; i < n; i++)
    {
        vel = parts[i].getAverageVelocity();
        vels_out << i << " " << sqrt(vel.sqr_len()) << " " << vel.x <<
        " " << vel.y << " " << vel.z << "\n";
    }
    vels_out.close();
}


void field::count_interactions() {
    potential_enery = 0;
    kinetic_energy = 0;
    for (size_t i = 0; i < n; i++)
    {
        //std::cerr << i << "Done\n";
        kinetic_energy += parts[i].countKineticEnergy();
        vect force;
        for (size_t j = 0; j < n; j++)
        {
            if (i == j) continue;
            if (i < j)
            {
                potential_enery += parts[i].countPotentialEnergy(&parts[j], cell_size);
            }
            force += parts[i].countForce(&parts[j], cell_size);
        }
        parts[i].setAceleration(force);
    }
    full_energy = potential_enery + kinetic_energy;
}

void field::updatePositions() {
    for (size_t i = 0; i < n; i++)
    {
        parts[i].update(false, dt, cell_size);
    }
}



void field::tick() {
    count_interactions();
    updatePositions();
}

double field::countAverageDr2() {
    double sum = 0;
    for (size_t i = 0; i < n; i++)
    {
        sum += parts[i].countDr2();
    }
    return sum / n;
}