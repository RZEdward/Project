#include "matplotlibcpp.h"
#include<iostream>
#include<cstdlib>
#include<Windows.h>
#include<ctime>
#include<iomanip>
#include<string>
#include<random>
#include<cmath>
namespace plt = matplotlibcpp;
using std::cout;
using std::endl;

#define PI 3.14159265

/*length function
random start
timestep
force to each box
notes
github*/

double length(double x1, double x2, double y1, double y2)
{
    double distance = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
    return distance;
}

int main()
{

    double boxlims = 1000.0; //Define The Size Of The Box
    int intboxlims = boxlims; //int version to add to string labels
    double origin = 0.0;
    int timesteps = 10000;
    int num_particles = 100;
    double start_pos = boxlims/2;
    
    double l_o = boxlims/20; //Maximum Distance Threshold For Spring To Take Effect
    int int_l_o = l_o;

    int i; //Loop
    int j; //Loop
    int n; //Loop
    double l; //Distance Between Particles
    double k = 2.0; //Constant for repulsion force F = k(l-l_o)
    double theta; //Angle Between Particles

    double u = 1; //Drag coefficient

    double dt = 1.0/40; //kind of 40fps

    double mean_motion = 0.0;
    double std_dev_motion = 25.0; //RNG Parameters
    double mean_start = start_pos;
    double std_dev_start = boxlims/20;

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(mean_start, std_dev_start); // RNG start position

    std::normal_distribution<double> distribute(mean_motion, std_dev_motion); // RNG motion

    std::vector<double> x_pos(num_particles);
    std::vector<double> y_pos(num_particles);

    for (i = 0; i < num_particles; i++)
    {
        x_pos[i] = distribution(generator);
        y_pos[i] = distribution(generator);
    }

    std::vector<int> empty = {}; //used to remove axis ticks later

    std::vector<int> t(timesteps); //mean square displacement graph x-axis
    std::vector<double> mean_square_displacement(timesteps); //mean square displacement graph y-axis

    std::string one = "Box Dimensions: ";
    std::string two = std::to_string(intboxlims);
    std::string three = " by ";
    std::string four = std::to_string(intboxlims);
    std::string five = "\nParticles Interact Within ";
    std::string six = std::to_string(int_l_o);
    std::string seven = " Units";
    std::string label = one + two + three + four + five + six + seven;

    std::string a = "Brownian Motion + Simple Harmonic Repulsion - ";
    std::string b = std::to_string(num_particles);
    std::string c = " Particles";
    std::string title = a + b + c;

    for (i = 0; i < timesteps; i++)
    {
        t[i] = i*dt;

        if (i == 0 || i == 10 || i % 100 == 0) //only show every 10th figure
        {

            plt::figure();
            plt::title(title);
            plt::xlim(origin,boxlims);
            plt::ylim(origin,boxlims);
            plt::xticks(empty);
            plt::yticks(empty);
            plt::xlabel(label);
            plt::plot(x_pos, y_pos, "ko");
            std::string filename = "C:/Users/ross/OneDrive/Desktop/Uni Work/Project/projectslides/slide";
            std::string num = std::to_string(i+1);
            std::string filetype = ".png";
            std::string savepoint = filename + num + filetype;
            plt::save(savepoint);

        }
        
        std::vector<double> shift_x(num_particles); //Total shift considering random + repulsion components
        std::vector<double> shift_y(num_particles);
        std::vector<double> x_rand(num_particles); //Randomly Generated Motion
        std::vector<double> y_rand(num_particles);
        std::vector<double> F_x(num_particles); //Repulsion Force
        std::vector<double> F_y(num_particles);
        
        for (j = 0; j < num_particles; j++)
        {

            x_rand[j] = distribute(generator);
            y_rand[j] = distribute(generator);

            for (n = 0; n < num_particles; n++)
            {

                if (n == j) //Particle can't interact with itself
                {
                    continue;
                }

                double x_n = x_pos[n];
                double y_n = y_pos[n];

                double x_j = x_pos[j];
                double y_j = y_pos[j];

                l = length(x_j, x_n, y_j, y_n);

                if (length(x_j, x_n + boxlims, y_j, y_n) < l) //right box
                {
                    l = length(x_j, x_n + boxlims, y_j, y_n);
                    x_n = x_n + boxlims;
                }

                else if (length(x_j, x_n - boxlims, y_j, y_n) < l) //left box
                {
                    l = length(x_j, x_n - boxlims, y_j, y_n);
                    x_n = x_n - boxlims;
                }

                else if (length(x_j, x_n, y_j, y_n + boxlims) < l) //top box
                {
                    l = length(x_j, x_n, y_j, y_n + boxlims);
                    y_n = y_n + boxlims;
                }

                else if (length(x_j, x_n, y_j, y_n - boxlims) < l) //bottom box
                {
                    l = length(x_j, x_n, y_j, y_n - boxlims);
                    y_n = y_n - boxlims;
                }

                else if (length(x_j, x_n + boxlims, y_j, y_n + boxlims) < l) //top-right box
                {
                    l = length(x_j, x_n + boxlims, y_j, y_n + boxlims);
                    x_n = x_n + boxlims;
                    y_n = y_n + boxlims;
                }

                else if (length(x_j, x_n - boxlims, y_j, y_n + boxlims) < l) //top-left box
                {
                    l = length(x_j, x_n - boxlims, y_j, y_n + boxlims);
                    x_n = x_n - boxlims;
                    y_n = y_n + boxlims;
                }

                else if (length(x_j, x_n - boxlims, y_j, y_n - boxlims) < l) //bottom-left box
                {
                    l = length(x_j, x_n - boxlims, y_j, y_n - boxlims);
                    x_n = x_n - boxlims;
                    y_n = y_n - boxlims;
                }

                else if (length(x_j, x_n + boxlims, y_j, y_n - boxlims) < l) //bottom-right box
                {
                    l = length(x_j, x_n + boxlims, y_j, y_n - boxlims);
                    x_n = x_n + boxlims;
                    y_n = y_n - boxlims;
                }
                
                if (l > l_o) //only allow close particles to interact
                {
                    continue;
                }

                double arg = (y_pos[j] - y_n)/(x_pos[j] - x_n);
                theta = atan(arg);
        
                if (x_pos[j] - x_n > 0)
                {
                    F_x[j] += abs(cos(theta)*k*(l - l_o)); //abs() used to ensure correct force direction
                }
                    
                else if (y_pos[j] - y_n > 0)
                {
                    F_y[j] += abs(sin(theta)*k*(l - l_o));
                }
                
                else
                {
                    F_x[j] += cos(theta)*k*(l - l_o);
                    F_y[j] += sin(theta)*k*(l - l_o);
                }
                
            }

            shift_x[j] += x_rand[j]*sqrt(dt) + (1/u)*F_x[j]*dt; //still don't fully understand use of sqrt(dt)
            shift_y[j] += y_rand[j]*sqrt(dt) + (1/u)*F_y[j]*dt;

        }

        double count = 0.0; //used to count total distance for all particles to box center

        for (j = 0; j < num_particles; j++)
        {
            if (x_pos[j] + shift_x[j] > boxlims)
            {
                x_pos[j] = x_pos[j] + shift_x[j] - boxlims; //periodic boundary conditions
            }

            else if (x_pos[j] + shift_x[j] < 0)
            {
                x_pos[j] = x_pos[j] + shift_x[j] + boxlims;
            }

            else if (y_pos[j] + shift_y[j] > boxlims)
            {
                y_pos[j] = y_pos[j] + shift_y[j] - boxlims;
            }

            else if (y_pos[j] + shift_y[j] < 0)
            {
                y_pos[j] = y_pos[j] + shift_y[j] + boxlims;
            }

            else
            {
                x_pos[j] = x_pos[j] + shift_x[j];
                y_pos[j] = y_pos[j] + shift_y[j];
            }

            double x_dist = abs(start_pos - x_pos[j]);
            double y_dist = abs(start_pos - y_pos[j]);

            double sqr_disp = x_dist*x_dist + y_dist*y_dist;

            count += sqr_disp;
        }

        double m_s_d = count/num_particles;

        mean_square_displacement[i] = m_s_d;
    }

    std::vector<double> model_msd(timesteps);
    for (i = 0; i < timesteps; i++)
    {
        model_msd[i] = 2*i*std_dev_motion*std_dev_motion; //model for mean square displacement
    } 

    plt::figure();
    plt::title("Mean Square Displacement of All Particles");
    plt::xlabel("Time (s)");
    plt::ylabel("Mean Square Displacement");
    //plt::plot(t, model_msd, "b--");
    plt::plot(t, mean_square_displacement, "k-");
    plt::save("C:/Users/ross/OneDrive/Desktop/Uni Work/Project/projectslides/meansquaredisplacement.png");

    return 0;

}