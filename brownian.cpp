#include<iostream>
#include<cstdlib>
#include<Windows.h>
#include<ctime>
#include<iomanip>
#include<string>
#include<random>
#include<cmath>
#include <fstream>
#include <vector>
#include <utility>
using std::cout;
using std::endl;

#define PI 3.14159265

int main()

{

    double origin = 0.0;
    int timesteps = 262;
    int N = 500;

    double R1 = 1.0;
    double boxlims = 500;
    int i, j, n, p;
    double dt = 0.1;
    double std_dev = 10;

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::default_random_engine generator;
    std::normal_distribution<double> normal(0.0,std_dev);

    std::vector<double> x_pos(N);
    std::vector<double> y_pos(N);   

    for (i = 0; i < N; i++)
        {
            x_pos[i] = boxlims/2;
            y_pos[i] = boxlims/2;
        }

    std::ofstream File3("brownianinfo.csv");
    File3 << "timesteps" << "," << "dt" << "," << "boxlims" << "," << "N" << "," << "R1" << endl;
    File3 << timesteps << "," << dt << "," << boxlims << "," << N << "," << R1 << endl;
    File3.close();

    std::ofstream File1("brownian_positions.csv");

    std::string labelx = "x_part";
    std::string labely = "y_part";
    for (i = 0; i < x_pos.size() - 1; i++)
    {       
        std::string iteration1 = std::to_string(i+1);
        File1 << (labelx + iteration1) << ", " << (labely + iteration1) << ", ";
    }
    std::string iteration2 = std::to_string(x_pos.size());
    File1 << (labelx + iteration2) << ", " << (labely + iteration2) << endl;

    for (i = 0; i < x_pos.size() - 1; i++)
    {
            File1 << x_pos[i] << "," << y_pos[i] << ",";
    }
    File1 << x_pos[x_pos.size() - 1] << "," << y_pos[x_pos.size() - 1] << endl;


    std::vector<double> t(timesteps);

    for (i = 0; i < timesteps; i++)
    {

        t[i] = i*dt;

        std::vector<double> shift_x(N); 
        std::vector<double> shift_y(N);

        for (j = 0; j < N; j++)
        {
            double x_rand = normal(generator);
            double y_rand = normal(generator);
            

            shift_x[j] = x_rand*sqrt(dt);
            shift_y[j] = y_rand*sqrt(dt);
        
        }

        for (j = 0; j < N; j++)
        {
            if (x_pos[j] + shift_x[j] > boxlims)
            {
                x_pos[j] = x_pos[j] + shift_x[j] - boxlims;
            }

            else if (x_pos[j] + shift_x[j] < 0)
            {
                x_pos[j] = x_pos[j] + shift_x[j] + boxlims;
            }

            else
            {
                x_pos[j] = x_pos[j] + shift_x[j];
            }

            if (y_pos[j] + shift_y[j] > boxlims)
            {
                y_pos[j] = y_pos[j] + shift_y[j] - boxlims;
            }

            else if (y_pos[j] + shift_y[j] < 0)
            {
                y_pos[j] = y_pos[j] + shift_y[j] + boxlims;
            }

            else
            {
                y_pos[j] = y_pos[j] + shift_y[j];
            }
        }

        for (j = 0; j < x_pos.size() - 1; j++)
        {
            File1 << x_pos[j] << "," << y_pos[j] << ",";
        }
        File1 << x_pos[x_pos.size() - 1] << "," << y_pos[x_pos.size() - 1] << endl;

    }

    File1.close();

    return 0;

}