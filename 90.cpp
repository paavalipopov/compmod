#include <iostream>
#include <ctime>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <string>

using namespace std;

double stand() {
    return 1.0 * rand()/RAND_MAX;
}


class Chamber {
    int n, n_start;
    int time_av, flag;
    double dispersion, n_av, n_sq_av;
    fstream results;

public:
    Chamber(int n, string fout) {
        this->n = n;
        n_start = n;
        dispersion = n_av = n_sq_av = flag = 0;
        results.open(fout);
        results << "Number of molecules - " << this->n << endl << endl;
        results << "Step \tn" << endl;
    }

    Chamber() {};

    ~Chamber() {
        results << endl << endl;
        results.close();
    };

    void make_step(int step) {
        if(stand() <= 1.0*n/n_start)
			n--;
		else
			n++;

        if(n <= n_start/2 && flag == 0) {
            flag = 1;
            time_av = step;
            n_av = 0;
        }
        if(flag == 1) {
            n_av += n;
            n_sq_av += n*n;
        }
        results << ' ' << step << "\t" << n << endl;
    }



    void average_calculation(int final_step) {
        n_av = n_av / (final_step - time_av);
        n_sq_av = n_sq_av / (final_step - time_av);
        dispersion = n_sq_av - n_av * n_av;
    }

    double get_dispersion() {
        return dispersion;
    }

    double get_n_av() {
        return n_av;
    }

    double get_n_sq_av() {
        return n_sq_av;
    }

    void print_final_state() {
        cout << "Number of molecules - " << n_start << endl;
        cout << "n average - " << n_av << endl;
        cout << "Dispersion - " << dispersion << endl;
        results << endl << "n average" << n_av << "\t dispersion - " << dispersion << endl << endl;
    }
};

int main(void) {
    srand(time(0));

    const int number_of_molecules = 160;
    const int number_of_ensembles = 100;

    double dispersion = 0;
    double n_av = 0;
    double n_sq_av = 0;


	fstream fs("results.txt", ios::out);
	fs.close();


    for(int i = 0; i < number_of_ensembles; i++) {
        cout << endl << endl << "Experiment #" << i <<endl;
        Chamber cmb(number_of_molecules, "results.txt");

        for(int t = 0; t < 10000; t++) {
            cmb.make_step(t);
        }
        cmb.average_calculation(9999);
        cmb.print_final_state();

        n_av += cmb.get_n_av();
        n_sq_av += cmb.get_n_sq_av();

        cout << "dispersion / n average = " << sqrt(cmb.get_dispersion()) / cmb.get_n_av() << endl;
        cout << "1 / sqrt(n) = " << 1.0 / sqrt(number_of_molecules) << endl;

    }


    n_av /= number_of_ensembles;
    n_sq_av /= number_of_ensembles;
    dispersion  = n_sq_av - n_av * n_av;

    cout << endl << endl << "Final results for ensembles" << endl << "n average = " << n_av << endl;
    cout << "n sq average = " << n_sq_av << endl;
    cout << "dispersion = " << dispersion << endl;
    cout << "sqrt(dispersion) / n average = " << sqrt(dispersion)/n_av << endl;
    cout << "1 / sqrt(n) = " << 1.0 / sqrt(number_of_molecules) << endl;

    return 0;
}
