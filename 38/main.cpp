#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
#include <string>

#define M_PI 3.14159265359

using namespace std;

class Ssalc {
public:
    Ssalc(string path, int N, double length) {
        this->path = path;
        this->N = N;
        this->length = length;

        calculations();

        printResults();
    }

    ~Ssalc() {
        fourierF.clear();
        fourierU.clear();
        uAnalytical.clear();
    }

    vector<double> getU() {
        return u;
    }

    vector<double> getUAnalytical() {
        return uAnalytical;
    }

private:
    int N;
    vector<double> fourierF;
    vector<double> fourierU;
    vector<double> u;
    vector<double> uAnalytical;

    double length;
    double step;
    ofstream results;
    string path;

    void calculations() {
        fourierF = vector<double>(N);
        fourierU = vector<double>(N);
        u = vector<double>(N);
        uAnalytical = vector<double>(N);

        step = length / N;

        for(int m = 1; m < N-1; m++)
            findFourierFAndU(m);

        findU();
    }

    void findFourierFAndU(int& m) {
        double Fm = 0;
        double Um = 0;

        for(int j = 1; j < N; ++j)
            Fm += heatFunction(step * j) * sin(1.0 * M_PI * m  * j / N);

        Fm *= 2.0 / N;
        Um = -1.0 * (length * length) * Fm / (M_PI * M_PI * m * m);

        fourierF[m] = Fm;
        fourierU[m] = Um;
    }

    void findU() {
        for(int j = 1; j < N; ++j) {
            for(int m = 1; m < N; ++m)
                u[j] += fourierU[m] * sin(1.0 * M_PI * m * j / N);
            uAnalytical[j] = exp( -1.0 / (1 - pow(j * step - 1, 2)));
        }
    }

    double heatFunction(const double& x) {
        const double y = 1 - x;

        const double f = exp( -1.0 / (1 - y*y) ) * (6 * pow(y, 4) - 2) / (pow(1 - y*y, 4));

        return f;
    }

    void printResults() {
        results.open(path, ios_base::trunc);

        results << "Step\tU\tU Analytical\tDifference" << endl;

        for(int i  = 0; i < N; i++)
            results << i << "\t" << u[i] << "\t" << uAnalytical[i] << "\t" << fabs(u[i] - uAnalytical[i]) << endl;

        results.close();
    }
};


int main()
{
    int N;
    double length = 2;
    string path;
    string yes;
    do {
        cout << "N = ";
        cin >> N;
        path = "results.txt";

        Ssalc experiment(path, N, length);

        cout << "Done. Repeat? ";
        cin >> yes;
    }
    while(yes == "y" || yes == "yes");

    return 0;
}
