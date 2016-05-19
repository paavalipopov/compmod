#include <iostream>
#include <string>
#include <cmath>
#include <fstream>
#include <vector>
#include "ffta.h"
#include "utilities.h"

using namespace std;
const double TwoPi = 6.283185307179586;

void FFTAnalysis(vector<double> & AVal, vector<double> & FTvl, int Nvl, int Nft) {
  int i, j, n, m, Mmax, Istp;
  double Tmpr, Tmpi, Wtmp, Theta;
  double Wpr, Wpi, Wr, Wi;
  double *Tmvl;

  n = Nvl * 2; Tmvl = new double[n];

  for (i = 0; i < n; i+=2) {
   Tmvl[i] = 0;
   Tmvl[i+1] = AVal[i/2];
  }

  i = 1; j = 1;
  while (i < n) {
    if (j > i) {
      Tmpr = Tmvl[i]; Tmvl[i] = Tmvl[j]; Tmvl[j] = Tmpr;
      Tmpr = Tmvl[i+1]; Tmvl[i+1] = Tmvl[j+1]; Tmvl[j+1] = Tmpr;
    }
    i = i + 2; m = Nvl;
    while ((m >= 2) && (j > m)) {
      j = j - m; m = m >> 1;
    }
    j = j + m;
  }

  Mmax = 2;
  while (n > Mmax) {
    Theta = -TwoPi / Mmax; Wpi = sin(Theta);
    Wtmp = sin(Theta / 2); Wpr = Wtmp * Wtmp * 2;
    Istp = Mmax * 2; Wr = 1; Wi = 0; m = 1;

    while (m < Mmax) {
      i = m; m = m + 2; Tmpr = Wr; Tmpi = Wi;
      Wr = Wr - Tmpr * Wpr - Tmpi * Wpi;
      Wi = Wi + Tmpr * Wpi - Tmpi * Wpr;

      while (i < n) {
        j = i + Mmax;
        Tmpr = Wr * Tmvl[j] - Wi * Tmvl[j-1];
        Tmpi = Wi * Tmvl[j] + Wr * Tmvl[j-1];

        Tmvl[j] = Tmvl[i] - Tmpr; Tmvl[j-1] = Tmvl[i-1] - Tmpi;
        Tmvl[i] = Tmvl[i] + Tmpr; Tmvl[i-1] = Tmvl[i-1] + Tmpi;
        i = i + Istp;
      }
    }

    Mmax = Istp;
  }

  for (i = 0; i < Nft; i++) {
    j = i * 2; FTvl[i] = 2*sqrt(pow(Tmvl[j],2) + pow(Tmvl[j+1],2))/Nvl;
  }

  delete []Tmvl;
}

class Doings {
    vector<double> x;
    vector<double> v;
    double v0;
    double kv[4], kx[4];
    double omega2, period, elementary_step, epsilon;
    ofstream results;
    string path;
    int number_of_steps, period_in_steps;

public:


    Doings(string path, double epsilon = 0.1, double omega2 = 5, double elementary_time = 0.1, int number_of_steps = 1000) {
        if(epsilon > 2*omega2) {
            cout << "Epsilon is greater than 2 * omega squared, automatic reduction" << endl;
            epsilon = omega2 * 1.9;
        }
        x.push_back(0);
        v.push_back(sqrt(4*omega2 - 2*epsilon));
        v0 = v[0];
        this->omega2 = omega2;
        this->elementary_step = elementary_time;
        this->number_of_steps = number_of_steps;
        this->path = path;
        this->epsilon = epsilon;
    }


    void find_elementary_step(int pow_2, double add){
        elementary_step += add;
        calculation();

        if (period_in_steps == pow_2)
            return;
        if (period_in_steps < pow_2 && add < 0) {
            find_elementary_step(pow_2, add);
        };
        if (period_in_steps < pow_2 && add > 0) {
            find_elementary_step(pow_2, -add/2);
        };
        if (period_in_steps > pow_2 && add < 0) {
            find_elementary_step(pow_2, -add/2);
        };
        if (period_in_steps > pow_2 && add > 0) {
            find_elementary_step(pow_2, add);
        };
        return;
    }

    int calculation(int number_of_stepr = 100000) {
        this->number_of_steps = number_of_stepr;
        results.open(path, ios_base::trunc);
        //results << endl << endl << "V0 = " << v[0] << "\tOmega squared = " << this->omega2 << "\t Epsilon = " << epsilon << endl;
        v.clear();
        x.clear();
        v.push_back(v0);
        x.push_back(0);
        period_in_steps = 0;
        results << "Step\t" << "V\t" << "X\t" << endl;
        for(int i = 0; i < number_of_steps; i++) {
            results << i << "\t" << v[i] << "\t" << x[i] << endl;

            kv[0] = formula_for_v(x[i], elementary_step);
            kx[0] = v[i] * elementary_step;
            kv[1] = formula_for_v(x[i] + kx[0]/2, elementary_step);
            kx[1] = (v[i] + kv[0]/2) * elementary_step;
            kv[2] = formula_for_v(x[i] + kx[1]/2, elementary_step);
            kx[2] = (v[i] + kv[1]/2) * elementary_step;
            kv[3] = formula_for_v(x[i] + kx[2], elementary_step);
            kx[3] = (v[i] + kv[2]) * elementary_step;

            v.push_back(v[i] + (kv[0] + 2*kv[1] + 2*kv[2] + kv[3])/6);
            x.push_back(x[i] + (kx[0] + 2*kx[1] + 2*kx[2] + kx[3])/6);
        }

        for(int i = 0; i < number_of_steps-1; i++) {
            if( fabs(x[i]) > fabs(x[i+1]) ) {
                period_in_steps = i*4;
                period = period_in_steps * elementary_step;
                break;
            }
        }

        if(period_in_steps < 1) {
            results << "Can't find period" << endl;
            cout << "Can't find period" << endl;
        }
        else {
           // results << " Period in steps = " << period_in_steps << endl;
            cout << "Period in steps = " << period_in_steps << endl;
        }

        results.close();
        return period_in_steps;
    }

    double formula_for_v(double x, double elementary_step) {
        return -1 * elementary_step * (this->omega2) * sin(x);
    }

    ~Doings() {
        x.clear();
        v.clear();
    };

    double get_vel(int i) {
        return v[i];
    }

    int get_el_step() {
        return elementary_step;
    }
};

int pow_of_2(int number) {
        if (number == 0)
            return 0;
        return ( (pow_of_2(number/2)) + 1);
    }

int pow_2(int power) {
        if ( power == 0)
            return 1;
        return ( 2 * pow_2(power-1));
    }

void print_v_coeffs(string path, vector<double> v_coeffs, int last) {
        ofstream results;
        results.open(path, ios_base::trunc);
        //results << "Numberof harmonica\tCoefficients\tW" << endl;
        results << "Numberof harmonica\tW" << endl;
        for(int i = 0; i < last; i++)
            results << i << "\t" << v_coeffs[i]*v_coeffs[i] << endl;
            //results << i << "\t" << v_coeffs[i] << "\t" << v_coeffs[i]*v_coeffs[i] << endl;
}

int isPowerOfTwo (unsigned int x) {
 return (
   x == 1 || x == 2 || x == 4 || x == 8 || x == 16 || x == 32 ||
   x == 64 || x == 128 || x == 256 || x == 512 || x == 1024 ||
   x == 2048 || x == 4096 || x == 8192 || x == 16384 ||
   x == 32768 || x == 65536 || x == 131072 || x == 262144 ||
   x == 524288 || x == 1048576 || x == 2097152 ||
   x == 4194304 || x == 8388608 || x == 16777216 ||
   x == 33554432 || x == 67108864 || x == 134217728 ||
   x == 268435456 || x == 536870912 || x == 1073741824 ||
   x == 2147483648);
}

int main(void) {
    double elementary_time;
    int period;
    double epsilon, omega2;
    cout << "Omega^2 = ";
    cin >> omega2;
    cout << endl << endl;

    do {
    cout << "Epsilon = ";
    cin >> epsilon;
    cout << "Omega^2 = " << omega2 << endl;

    string yes;
    do {
        cout << "Enter elementary time: ";
        cin >> elementary_time;
        Doings test("test.txt", epsilon, omega2, elementary_time, 100000);
        period = test.calculation();
        cout << "Ok? - ";
        cin >> yes;
    }
    while (yes != "y");

    Doings result("result.txt", epsilon, omega2, elementary_time, 100000);
    period = result.calculation();

    int tmp = pow_2(pow_of_2(period));
    cout << "pow 2 = " << tmp << endl;
    result.find_elementary_step(tmp, -elementary_time/3);
    period = result.calculation( (int) (tmp*1.2) );


    vector<double> velocities;
    for(int i = 0; i < period; i++)
        velocities.push_back(result.get_vel(i));


    vector <complex<double> > vel_comp = convert(velocities);
    vector <complex<double> > v_cf_com = ffta(vel_comp, false);
    vector <double> v_cf = convert(v_cf_com);


    print_v_coeffs("vel_energy.txt", v_cf, period);

    cout << "DONE" << endl << endl;

    }

    while(true);

    return 0;
}

