#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <stdbool.h>

struct chamber {
    int n, m;
    int time_av;
    double dispersion, n_av, n_sq_av;
};

int main(void) {
    FILE* results;
    const int N = 160;
    const int number_of_chambers = 3;
    int time_end = 1000;
    int j, t, flag;
    bool molecule[N];
    struct chamber chamber[number_of_chambers];

    results = fopen("./results.txt", "w");
    srand(time(NULL));

    for(int i = 0; i < number_of_chambers; i++) {
        flag = 1;
        chamber[i].n = N;
        chamber[i].n_av = chamber[i].n_sq_av = 0;
        fprintf(results, "Cycle %d\t N = %d\n", i+1, N);
        for(j = 0; j < N; j++){
            molecule[j] = false;
        }

        fprintf(results, "Time\t n\n");
        for(t = 0; t < time_end; t++) {
            j = (rand() * rand()) % N;
            if(molecule[j]) {
                molecule[j] = false;
                chamber[i].n++;
            }
            else {
                molecule[j] = true;
                chamber[i].n--;
            }
            if(chamber[i].n < N/1.8 && flag == 1){
                chamber[i].time_av = t;
                flag = 0;
            }
            if(flag == 0) {
                chamber[i].n_av += chamber[i].n;
                chamber[i].n_sq_av += chamber[i].n*chamber[i].n;
            }
            fprintf(results, "%d\t %d\n", t, chamber[i].n);
        }
        chamber[i].n_av = chamber[i].n_av / (double)(time_end - chamber[i].time_av);
        chamber[i].n_sq_av = chamber[i].n_sq_av / (double)(time_end - chamber[i].time_av);
        chamber[i].dispersion = chamber[i].n_sq_av - chamber[i].n_av * chamber[i].n_av;
        fprintf(results, "n average\t dispersion\t relative fluctuation\t relative fluctuation - 1/sqrt(N)\n");
        fprintf(results, "%f\t %f\t %f\t %f\n", chamber[i].n_av, chamber[i].dispersion, chamber[i].dispersion / chamber[i].n_av, chamber[i].dispersion / chamber[i].n_av - 1.0/sqrt(N));
        printf("Cycle %d, N = %d\nn average\t\t\t %0.2f\ndispersion\t\t\t %0.2f\nrelative fluctuation\t\t %0.2f\nrelative fluctuation - 1/sqrt(N) %0.2f\n\n\n", i+1, N, chamber[i].n_av, chamber[i].dispersion, chamber[i].dispersion / chamber[i].n_av, chamber[i].dispersion / chamber[i].n_av - 1.0/sqrt(N));
        fprintf(results, "\n\n\n");
    }
    fclose(results);
    return 0;
}
