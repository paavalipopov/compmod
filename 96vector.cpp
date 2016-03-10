#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <vector>
using namespace std;


double stand() {
    return 1.0*rand()/RAND_MAX;
}
double exp_rand() {
    return -log(1.0 - stand());
}
int number_of_neutrons() {
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
            return number_of_neutrons();
        else
            return 0;
    }
};
int task_one(FILE* task);
int task_two(FILE* task);

int main(void) {
    srand(time(NULL));
    FILE* task = fopen("./task.txt", "w");
    //task_one(task);
    task_two(task);
    fclose(task);

    return 0;
}

int task_one(FILE* task) {
    const long time = 30, number_of_experiments = 500;
    double r;
    long good_exp;
    long number_radius, number_exp, exp_time, expl, expl_for_born, new_born, total_new_born, total_born;
    vector<Neutron> line1;
    vector<Neutron> line2;
    fprintf(task, "First generation: 1 neutron in the center\n\n");

    for(number_radius = 1; number_radius < 6; number_radius++) {
        good_exp = 0;
        r = 4.0*number_radius;
        cout << "R - " << r << endl;
        //fprintf(task, "Radius \t%f\n", r);
        for(number_exp = 0; number_exp < number_of_experiments; number_exp++){
            //fprintf(task, "Experiment \t%d\n", number_exp);
            new_born = total_new_born = 0;
            line1.push_back(Neutron(0, 0, 0, r));
            total_born = 1;


            for(exp_time = 0; exp_time < time; exp_time++) {

                //cout << "Time -  " << exp_time << endl;
                for(expl = 0; expl < total_born; expl++) {
                    //cout << "Neutron # - " << expl << endl;
                    line1[expl].move_neutron();
                    new_born = line1[expl].check();
                    //cout << "New born neutron for it - " << new_born << endl;
                    for(expl_for_born = 0; expl_for_born < new_born; expl_for_born++) {
                        line2.push_back(Neutron(line1[expl]));
                    }
                    total_new_born += new_born;
                    new_born = 0;


                }

                total_born = total_new_born;
                //cout << "Total born - " << total_born << endl;
                total_new_born = 0;
                //fprintf(task, "Time \t%d \tNeutrons in next generation \t%d\n", exp_time, total_born);
                line1.clear();
                line1 = line2;
                line2.clear();
                if(total_born == 0){
                    good_exp++;
                    //cout << "Good exp - " << good_exp << endl;
                    break;
                }
            }
            line1.clear();
        }
        cout << "Probability - " << (1 - (double)good_exp/number_of_experiments) << endl;
        fprintf(task, "Radius - %f, probability - %f\n\n", r, 1.0 - (double)good_exp/number_of_experiments);
    }
}

int task_two(FILE* task) {
    const long time = 20, number_of_experiments = 500;
    const double r = 12;
    long old_generation, new_generation;
    double ratio_gen = 0;

    long number_exp, exp_time, expl, expl_for_born, new_born, total_new_born, total_born;
    vector<Neutron> line1;
    vector<Neutron> line2;

    cout << "R - " << r << endl;

    for(number_exp = 0; number_exp < number_of_experiments; number_exp++){
            for(int i = 0; i < 50; i++) {
                line1.push_back(Neutron(r));
            }
            new_born = total_new_born = 0;
            new_generation = 0;
            total_born = old_generation = 50;

            for(exp_time = 0; exp_time < time; exp_time++) {
                    for(expl = 0; expl < total_born; expl++) {
                    //cout << "Neutron # - " << expl << endl;
                    line1[expl].move_neutron();
                    new_born = line1[expl].check();
                    //cout << "New born neutron for it - " << new_born << endl;
                    for(expl_for_born = 0; expl_for_born < new_born; expl_for_born++) {
                        line2.push_back(Neutron(line1[expl]));
                    }
                    total_new_born += new_born;
                    new_born = 0;
                }

                new_generation = total_born = total_new_born;
                //cout << "New gen - " << new_generation << endl;
                //cout << "Total born - " << total_born << endl;
                total_new_born = 0;
                //fprintf(task, "Time \t%d \tNeutrons in next generation \t%d\n", exp_time, total_born);
                line1.clear();
                line1 = line2;
                line2.clear();
                //cout << "Old gen - " << old_generation << endl;
                ratio_gen += (double)new_generation / old_generation;
                old_generation = new_generation;
                //cout << "Ratio - " << (ratio_gen) << endl;
                if(total_born == 0){
                    //cout << "Good exp - " << good_exp << endl;
                    break;
                }
            }
            line1.clear();
        }
        ratio_gen /= (time * number_of_experiments);
        cout << "Ratio - " << (ratio_gen) << endl;
        fprintf(task, "Radius - %f, ratio - %f\n\n", r, ratio_gen);
}



