#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector>
using namespace std;


double stand() {
    return rand()/RAND_MAX;
}
double exp_rand() {
    return -log(1 - stand());
}
int number_of_neutrons () {
    double rnd = stand();

    if(rnd <= 0.025)
        return 0;
    else if(rnd <= 0.855)
        return 1;
    else if(rnd <= 0.925)
        return 2;
    else if(rnd <= 0.975)
        return 3;
    else
        return 4;
}


class Vector {
public:
    double x, y, z;
    Vector(double x, double y, double z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }
    Vector(double r) {
        double R = r * pow(stand(), 0.33);
        double cos_theta = 2*stand() - 1;
        double phi = 2 * 3.14159265359 * stand();
        this->x = R * cos(phi) * cos_theta;
        this->y = R * sin(phi) * cos_theta;
        this->z = R * sqrt(1 - cos_theta*cos_theta);
    }
    Vector() {
        double R = exp_rand() * pow(stand(), 0.33);
        double cos_theta = 2*stand() - 1;
        double phi = 2 * 3.14159265359 * stand();
        this->x = R * cos(phi) * cos_theta;
        this->y = R * sin(phi) * cos_theta;
        this->z = R * sqrt(1 - cos_theta*cos_theta);
    }
    Vector operator+(const Vector & a) {
        return Vector(this->x + a.x, this->y + a.y, this->z + a.z);
    }
    double norma(){
        return x*x + y*y + z*z;
    }
};

class Neutron {
public:
    double r;
    Vector position;
    Neutron(double x, double y, double z, double r) {
        this->position = Vector(x, y, z);
        this->r = r;
    }
    Neutron(double r) {
        this->position = Vector(r);
        this->r = r;
    }
    Neutron(const Neutron & a) {
        this->position = a.position;
        this->r = a.r;
    }
    void move_neutron() {
        this->position = this->position + Vector();
    }
    int check() {
        if(position.norma() < r*r)
            return 1;
        else
            return 0;
    }
};
int task_one(FILE* task);
int task_two(FILE* task);

int main(void) {
    FILE* task = fopen("./task.txt", "w");
    task_one(task);
    task_two(task);
    fclose(task);

    return 0;
}

int task_one(FILE* task) {
    const long time = 100, number_of_experiments = 100;
    double r;
    long exp, good_exp;
    long number_radius, number_exp, exp_time, expl, expl_for_born, new_born, total_new_born, total_born;
    Neutron** neutron_array[time];
    fprintf(task, "First generation: 1 neutron in the center\n\n");

    for(number_radius = 1; number_radius < 6; number_radius++) {
        exp = good_exp = 0;
        r = 2.0*number_radius;
        fprintf(task, "Radius - %f\n", r);
        for(number_exp = 0; number_exp < number_of_experiments; number_exp++){
            exp++;
            printf("Radius - %lf\n", r);
            neutron_array[0][0] = new Neutron(0.0, 0.0, 0.0, r);
            printf("Radius - %lf\n", r);
            total_born = 1;
            for(exp_time = 0; exp_time < time-1; exp_time++) {
                for(expl = 0; expl < total_born; expl++) {
                    neutron_array[exp_time][expl]->move_neutron();
                    new_born = neutron_array[exp_time][expl]->check();
                    for(expl_for_born = 0; expl_for_born < new_born; expl_for_born++) {
                        neutron_array[exp_time+1][total_new_born + expl_for_born] = new Neutron(*(neutron_array[exp_time][expl]));
                    }
                    total_new_born += new_born;
                    new_born = 0;
                    delete neutron_array[exp_time][expl];
                }
                total_born = total_new_born;
                total_new_born = 0;
                if(total_born == 0){
                    good_exp++;
                    break;
                }
            }
        }
        printf("Radius - %lf, probability - %lf\n", r, 1 - (double)good_exp/exp);
        fprintf(task, "Radius - %lf, probability - %lf\n\n", r, 1 - (double)good_exp/exp);
    }
}

int task_two(FILE* task) {

}

