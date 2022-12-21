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
#define k 1.0 //Constant for repulsion force F = k(l-l_o)
#define u 1 //Drag coefficient


//easier to plot v_rms here - don't have to produce any movie slides/ visually represent anything


double length(double x1, double x2, double y1, double y2)
{
    double distance = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
    return distance;
}

int main()
{
    
    double boxlims;
    double origin = 0.0;
    int timesteps = 2000; //would do more but a) computer would die over 10000, b)particles reach equilibrium very quickly
    int num_particles = 20*20; //requires square number to evenly distribute particles in the box at the start

    int runs = 5; //tried to run from phi = 0.78 to 1 in steps of 0.02, but it bugs out after 9 iterations (phi = 0.94)

    std::vector<double> packing_fraction(runs);
    std::vector<double> boxlimits(runs);
    double phi;
    double r1 = 1.0;
    double r2 = 1.4*r1;
    double d1 = 2*r1;
    double d2 = 2*r2;

    int i; //Loop
    int j; //Loop
    int n; //Loop
    int p; //Loop
    int g; //Loop
    
    double l; //Distance Between Particles
    double l_o; //Maximum Distance Threshold For Spring To Take Effect
    double theta; //Angle Between Particles

    double dt = 1.0/10;

    std::vector<double> x_pos(num_particles);
    std::vector<double> y_pos(num_particles);
    std::vector<double> radii(num_particles);
    
    int numb_particles = num_particles; //Solution to use of odd number of particles
    if (num_particles % 2 == 1)
    {
        numb_particles = num_particles + 2;
    }

    std::string title1 = "Log_10 Root Mean Square Velocity - "; //label making
    std::string title2 = std::to_string(num_particles);
    std::string title3 = " Particles";
    std::string title = title1 + title2 + title3;

    plt::figure(1);
    plt::title(title);
    plt::xlabel("Time (s)");
    plt::ylabel("Log_10 $\\sqrt{<v^{2}>}$");
    std::vector<std::string> color = {"black", "grey", "red", "orange", "yellow", "green", "blue", "purple", "pink"};

    for (i = 0; i < packing_fraction.size(); i++) //first choose the packing fractions - then compute radii
    {
        packing_fraction[i] = 0.84 + i*0.02; //= N*pi*r^2/L^2
        boxlimits[i] = sqrt(PI*(numb_particles/2*r1*r1 + num_particles/2*r2*r2)/packing_fraction[i]);
    }

    for (p = 0; p < packing_fraction.size(); p++) //run all the code from repulsion.cpp - but for each packing fraction (and without displaying plots)
    {
        
        phi = packing_fraction[p];
        boxlims = boxlimits[p];
        cout << phi << endl; //just testing everything looks ok
        cout << boxlims << endl;

        double mean_motion = 0.0;
        double std_dev_motion = 2*boxlims; //RNG Parameters - make particles collide at the start
        std::random_device rd{};
        std::mt19937 gen{rd()};
        std::default_random_engine generator;
        std::normal_distribution<double> distribute(mean_motion, std_dev_motion); // RNG motion
        std::uniform_real_distribution<double> uniform(0.0,boxlims);

        for (i = 0; i < num_particles; i++) //define radius of each particle
        {

            if (i % 2 == 0)
            {
                radii[i] = r1;
            }
            else
            {
                radii[i] = r2;
            }

            /*double a = sqrt(num_particles); //defining uniformly distributed start positions
            int A = int(a);
            double b = 0;
            double z = 0;

            for (j = 0; j < A; j++)
            {
                if (i >= (A + (j-1)*A) && i < (A + j*A))
                {
                    b = j;
                }
            }

            x_pos[i] = boxlims/(2*sqrt(num_particles)) + b*boxlims/sqrt(num_particles);

            for (j = 0; j < A; j++)
            {
                if (i % int(A) == j)
                {
                    z = j;
                }
            }

            y_pos[i] = boxlims/(2*sqrt(num_particles)) + z*boxlims/sqrt(num_particles);*/
        
        }

        //uniformly random start positions
        for (i = 0; i < num_particles; i++)
        {
            x_pos[i] = uniform(generator);
            y_pos[i] = uniform(generator);
        }

        std::vector<double> t(timesteps); //mean square displacement graph x-axis data
        std::vector<double> mean_square_displacement(timesteps); //mean square displacement graph y-axis data
        std::vector<double> rms_velocity(timesteps); //root mean square velocity graph y-axis data
        std::vector<double> log_rms_velocity(timesteps); //log root mean square velocity graph y-axis data

        for (i = 0; i < timesteps; i++) //start of the big loop - updating positions + tracking velocity
        {

            t[i] = i*dt;

            std::vector<double> shift_x(num_particles); //Total shift considering random + repulsion components
            std::vector<double> shift_y(num_particles);
            std::vector<double> x_rand(num_particles); //Randomly Generated Motion
            std::vector<double> y_rand(num_particles);
            std::vector<double> F_x(num_particles); //Repulsion Force
            std::vector<double> F_y(num_particles);
            
            double v_rms_total = 0.0;
            for (j = 0; j < num_particles; j++) //particle of interest
            {

                x_rand[j] = distribute(generator); 
                y_rand[j] = distribute(generator);

                for (n = 0; n < num_particles; n++) //secondary particle that interacts
                {

                    if (n == j) //Particle can't interact with itself
                    {
                        continue;
                    }

                    l_o = radii[j] + radii[n]; //interactions only occur upon particle overlap

                    double x_n = x_pos[n];
                    double y_n = y_pos[n];

                    double x_j = x_pos[j];
                    double y_j = y_pos[j];

                    l = length(x_j, x_n, y_j, y_n);

                    std::vector<int> verify(8);

                    //next set of if statements checks for particles in neighbouring boxes close enough to interact 

                    if (length(x_j, x_n + boxlims, y_j, y_n) < l) //right box
                    {
                        l = length(x_j, x_n + boxlims, y_j, y_n);
                        verify[0] = 1;
                    }

                    if (length(x_j, x_n - boxlims, y_j, y_n) < l) //left box
                    {
                        l = length(x_j, x_n - boxlims, y_j, y_n);
                        verify[0] = 0;
                        verify[1] = 1;
                    }

                    if (length(x_j, x_n, y_j, y_n + boxlims) < l) //top box
                    {
                        l = length(x_j, x_n, y_j, y_n + boxlims);
                        verify[0] = 0;
                        verify[1] = 0;
                        verify[2] = 1;
                    }

                    if (length(x_j, x_n, y_j, y_n - boxlims) < l) //bottom box
                    {
                        l = length(x_j, x_n, y_j, y_n - boxlims);
                        verify[0] = 0;
                        verify[1] = 0;
                        verify[2] = 0;
                        verify[3] = 1;
                    }

                    if (length(x_j, x_n + boxlims, y_j, y_n + boxlims) < l) //top-right box
                    {
                        l = length(x_j, x_n + boxlims, y_j, y_n + boxlims);
                        verify[0] = 0;
                        verify[1] = 0;
                        verify[2] = 0;
                        verify[3] = 0;
                        verify[4] = 1;
                    }

                    if (length(x_j, x_n - boxlims, y_j, y_n + boxlims) < l) //top-left box
                    {
                        l = length(x_j, x_n - boxlims, y_j, y_n + boxlims);
                        verify[0] = 0;
                        verify[1] = 0;
                        verify[2] = 0;
                        verify[3] = 0;
                        verify[4] = 0;
                        verify[5] = 1;
                    }

                    if (length(x_j, x_n - boxlims, y_j, y_n - boxlims) < l) //bottom-left box
                    {
                        l = length(x_j, x_n - boxlims, y_j, y_n - boxlims);
                        verify[0] = 0;
                        verify[1] = 0;
                        verify[2] = 0;
                        verify[3] = 0;
                        verify[4] = 0;
                        verify[5] = 0;
                        verify[6] = 1;
                    }

                    if (length(x_j, x_n + boxlims, y_j, y_n - boxlims) < l) //bottom-right box
                    {
                        l = length(x_j, x_n + boxlims, y_j, y_n - boxlims);
                        verify[0] = 0;
                        verify[1] = 0;
                        verify[2] = 0;
                        verify[3] = 0;
                        verify[4] = 0;
                        verify[5] = 0;
                        verify[6] = 0;
                        verify[7] = 1;
                    }
                    
                    if (l > l_o) //only allow overlapping particles to interact
                    {
                        continue;
                    }

                    if (verify[0] == 1)
                    {
                        x_n += boxlims;
                    }
                    if (verify[1] == 1)
                    {
                        x_n -= boxlims;
                    }
                    if (verify[2] == 1)
                    {
                        y_n += boxlims;
                    }
                    if (verify[3] == 1)
                    {
                        y_n -= boxlims;
                    }
                    if (verify[4] == 1)
                    {
                        x_n += boxlims;
                        y_n += boxlims;
                    }
                    if (verify[5] == 1)
                    {
                        x_n -= boxlims;
                        y_n += boxlims;
                    }
                    if (verify[6] == 1)
                    {
                        x_n -= boxlims;
                        y_n -= boxlims;
                    }
                    if (verify[7] == 1)
                    {
                        x_n += boxlims;
                        y_n -= boxlims;
                    }

                /*double sin_theta = (y_pos[j] - y_n)/l;
                double cos_theta = (x_pos[j] - x_n)/l;
                F_x[j] += cos_theta*k*(l - l_o); //sum all the forces from one particle to every other particle
                F_y[j] += sin_theta*k*(l - l_o);*/

                double arg = (y_pos[j] - y_n)/(x_pos[j] - x_n);
                theta = atan2(arg); //angle between particles
        
                if (x_pos[j] - x_n > 0)
                {
                    F_x[j] += abs((cos(theta))*k*(l - l_o)); //Define strength of repulsion force - abs() used to ensure correct force direction
                }
                else
                {
                    F_x[j] += cos(theta)*k*(l - l_o); 
                }
                    
                if (y_pos[j] - y_n > 0)
                {
                    F_y[j] += abs((sin(theta))*k*(l - l_o));
                }
                else
                {
                     F_y[j] += sin(theta)*k*(l - l_o); //sum all the forces from one particle to every other particle
                }
        
                }

                /*if (i == 0)
                {
                    shift_x[j] += (1/u)*F_x[j]*dt + x_rand[j]*sqrt(dt);
                    shift_y[j] += (1/u)*F_y[j]*dt + y_rand[j]*sqrt(dt);
                }
                else
                {*/
                    shift_x[j] += (1/u)*F_x[j]*dt;
                    shift_y[j] += (1/u)*F_y[j]*dt;
                //}

                v_rms_total += shift_x[j]*shift_x[j] + shift_y[j]*shift_y[j];

            }

            double v_rms = v_rms_total/num_particles;
            rms_velocity[i] = v_rms;
            log_rms_velocity[i] = log10(v_rms);

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

            }

        } //End of i loop - timesteps

    std::vector<double> x(t.size()-1);
    std::vector<double> y(t.size()-1);
    for (g = 0; g < x.size(); g++) //dont want to show first step velocity (has the random component) - tried using t[1:] but doesn't seem to work on c++
    {
        x[g] = t[g+1];
        y[g] = rms_velocity[g+1];
    }

    std::string phii = std::to_string(packing_fraction[p]);
    plt::figure(1);
    plt::plot(x, y, {{"color", color[p]}, {"linestyle", "-"}, {"label", phii}});
    
    } //End of p loop - different packing fractions

    plt::legend({{"loc", "upper right"}, {"title", "$\\phi$"}});
    plt::save("C:/Users/ross/OneDrive/Desktop/Uni Work/Project/projectslides/rootmeansquarevelocity.png"); //save final figure

    return 0;

}