#include<iostream>
#include<cstdlib>
#include<ctime>
#include<iomanip>
#include<string>
#include<random>
#include<cmath>
#include <fstream>
#include <vector>
#include <utility>
#include <math.h>
using std::cout;
using std::endl;

#define PI 3.14159265

int main()

{

    double x1,x2,y1,y2,distance;
    double k = 1.0;
    double u = 1.0;
    double G = 1.0;
    double origin = 0.0;
    int timesteps = 1000;
    int num_particles = 20*20;
    int numb_particles = num_particles;
    double kB = 1.380649*pow(10,-23);
    double eps = 1.0;
    double sig = 1.7818;
    int I, J;

    if (num_particles % 2 == 1)
    {
        numb_particles = num_particles + 2;
    }

    double R1 = 1.0;
    double R2 = 1.0;
    double phi = 0.25;
    double boxlims = sqrt(PI*(numb_particles/2*R1*R1 + num_particles/2*R2*R2)/phi);
    double r, D;
    int i, j, n, p;
    double x_particle_n, y_particle_n, x_n, y_n, x_j, y_j;
    double dy, dx, theta;
    double dt = 0.025;
    double a = 0.5;

    double std_dev = 0.1*boxlims/sqrt(num_particles);

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::default_random_engine generator;
    std::uniform_real_distribution<double> uniform(0.0,boxlims);
    std::normal_distribution<double> distribution(0.0, std_dev);

    std::vector<double> x_pos(num_particles);
    std::vector<double> y_pos(num_particles);   
    std::vector<double> x_posR1(numb_particles/2);
    std::vector<double> y_posR1(numb_particles/2);
    std::vector<double> x_posR2(num_particles/2);
    std::vector<double> y_posR2(num_particles/2);
    std::vector<double> radii(num_particles);

    for (i = 0; i < num_particles; i++)
    {
        if (i % 2 == 0)
        {
            radii[i] = R1;
        }
        else
        {
            radii[i] = R2;
        }
    }

    for (i = 0; i < num_particles; i++)
        {
            //x_pos[i] = uniform(generator);
            //y_pos[i] = uniform(generator); //want no initial contact
           

            for (j = 0; j < sqrt(num_particles); j++)
            {
                if ( j*sqrt(num_particles) <= i && i < (j+1)*sqrt(num_particles) )
                {
                    I = j; 
                    continue;
                }
            }
            
            for (j = 0; j < sqrt(num_particles); j++)
            {
                if ( j*sqrt(num_particles) <= i && i < (j+1)*sqrt(num_particles) )
                {
                    J = i - j*sqrt(num_particles); 
                    continue;
                }
            }
                
            x_pos[i] = boxlims/(2*sqrt(num_particles)) + I*boxlims/sqrt(num_particles)  +  distribution(generator);

            y_pos[i] = boxlims/(2*sqrt(num_particles)) + J*boxlims/sqrt(num_particles)  +  distribution(generator);

        }



    for (i = 0; i < num_particles/2; i++) 
        {
            x_posR1[i] = x_pos[2*i];
            y_posR1[i] = y_pos[2*i];

            x_posR2[i] = x_pos[2*i+1];
            y_posR2[i] = y_pos[2*i+1];
        }
    
    if (num_particles % 2 == 1)
    {
        x_posR1[x_posR1.size()-1] = x_pos[x_pos.size()-1];
        y_posR1[x_posR1.size()-1] = y_pos[x_pos.size()-1];  
    }

    int numR1 = 2*x_posR1.size();
    int numR2 = 2*x_posR2.size();

    std::ofstream File3("LJproperties.csv");
    File3 << "timesteps" << "," << "boxlims" << "," << "numR1" << "," << "numR2" << "," << "R1" << "," << "R2" << "," << "phi" << endl;
    File3 << timesteps << "," << boxlims << "," << numR1 << "," << numR2 << "," << R1 << "," << R2 << "," << phi << endl;
    File3.close();

    std::ofstream File1("LJR1.csv");
    std::ofstream File2("LJR2.csv");

    std::string labelx = "x_part";
    std::string labely = "y_part";
    for (i = 0; i < x_posR1.size() - 1; i++)
    {       
        std::string iteration1 = std::to_string(i+1);
        File1 << (labelx + iteration1) << ", " << (labely + iteration1) << ", ";
    }
    std::string iteration2 = std::to_string(x_posR1.size());
    File1 << (labelx + iteration2) << ", " << (labely + iteration2) << endl;

    for (i = 0; i < x_posR2.size() - 1; i++)
    {       
        std::string iteration3 = std::to_string(i+1);
        File2 << (labelx + iteration3) << ", " << (labely + iteration3) << ", ";
    }
    std::string iteration4 = std::to_string(x_posR2.size());
    File2 << (labelx + iteration4) << ", " << (labely + iteration4) << endl;

    /*for (i = 0; i < x_posR1.size() - 1; i++)
    {
            File1 << x_posR1[i] << "," << y_posR1[i] << ",";
    }
    File1 << x_posR1[x_posR1.size() - 1] << "," << y_posR1[x_posR1.size() - 1] << endl;

    for (i = 0; i < x_posR2.size() - 1; i++)
    { 
        File2 << x_posR2[i] << "," << y_posR2[i] << ",";
    }
    File2 << x_posR2[x_posR2.size() - 1] << "," << y_posR2[x_posR2.size() - 1] << endl;*/


    std::vector<double> dts(timesteps);
    dts[0] = dt;

    for (i = 0; i < timesteps; i++)
    {

        // find the max velocity, then update dt by dt = a/max_v

        std::vector<double> shift_x(num_particles); 
        std::vector<double> shift_y(num_particles);
        std::vector<double> x_rand(num_particles);
        std::vector<double> y_rand(num_particles);
        std::vector<double> F_x(num_particles);
        std::vector<double> F_y(num_particles);

        double max_vel = 0;
        for (j = 0; j < num_particles; j++)
        {
            for (n = 0; n < num_particles; n++) 
            {
                if (n == j) 
                {
                    continue;
                }

                D = radii[j] + radii[n];

                x_n = x_pos[n];
                y_n = y_pos[n];

                x_j = x_pos[j];
                y_j = y_pos[j];

                std::vector<std::vector<double>> coords = {{x_j, x_n + boxlims, y_j, y_n + boxlims},{x_j, x_n + boxlims, y_j, y_n},{x_j, x_n + boxlims, y_j, y_n - boxlims},{x_j, x_n, y_j, y_n - boxlims},{x_j, x_n - boxlims, y_j, y_n - boxlims},{x_j, x_n - boxlims, y_j, y_n},{x_j, x_n - boxlims, y_j, y_n + boxlims},{x_j, x_n, y_j, y_n + boxlims}};
                r = sqrt((x_j - x_n)*(x_j - x_n) + (y_j - y_n)*(y_j - y_n));

                x_particle_n = x_n;
                y_particle_n = y_n;
                for (p = 0; p < coords.size(); p++)
                {
                    x1 = coords[p][0];
                    x2 = coords[p][1];
                    y1 = coords[p][2];
                    y2 = coords[p][3];
                    distance = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
                    if (distance < r)
                    {
                        r = distance;
                        x_particle_n = (coords[p])[1];
                        y_particle_n = (coords[p])[3];
                    }
                }

                x_n = x_particle_n;
                y_n = y_particle_n;
            

                dx = x_j - x_n;
                dy = y_j - y_n;
                theta = atan2(dy,dx);

                double R = radii[j];

                F_x[j] += cos(theta)*(24*eps*pow(sig,6)*(2*pow(sig,6)-pow(r,6)))/pow(r,13);         
                F_y[j] += sin(theta)*(24*eps*pow(sig,6)*(2*pow(sig,6)-pow(r,6)))/pow(r,13);
            }

            double vel = sqrt((1/u)*F_x[j]*(1/u)*F_x[j] + (1/u)*F_y[j]*(1/u)*F_y[j]);
            if (vel > max_vel)
            {
                max_vel = vel;
            }

        }

        double max_dt = 0.4;
        double base_dt = 0.025;

        double displacement = sqrt((1/u)*F_x[j]*dt*(1/u)*F_x[j]*dt + (1/u)*F_y[j]*dt*(1/u)*F_y[j]*dt);

        if (displacement > R1)
        {
            //dt = (R1/displacement)*base_dt;
            dt = base_dt;
        }
        else
        {
            dt = base_dt;
        }


        //dt = a/max_vel;
        if (i != 0)
        {
            dts[i] = dt;
        }

        for (j = 0; j < num_particles; j++)
        {
            //cout << "(" << F_x[j] << ", " << F_y[j] << ")    "; 
            shift_x[j] = (1/u)*F_x[j]*dt;
            shift_y[j] = (1/u)*F_y[j]*dt;
        
        }

        for (j = 0; j < num_particles/2; j++)
        {
            x_posR1[j] = x_pos[2*j];
            y_posR1[j] = y_pos[2*j];
            x_posR2[j] = x_pos[2*j+1];
            y_posR2[j] = y_pos[2*j+1];
        }

        if (num_particles % 2 == 1)
        {
            x_posR1[x_posR1.size()-1] = x_pos[x_pos.size()-1];
            y_posR1[x_posR1.size()-1] = y_pos[x_pos.size()-1];
        }

        for (j = 0; j < x_posR1.size() - 1; j++)
        {
            File1 << x_posR1[j] << "," << y_posR1[j] << ",";
        }
        File1 << x_posR1[x_posR1.size() - 1] << "," << y_posR1[x_posR1.size() - 1] << endl;

        for (j = 0; j < x_posR2.size() - 1; j++)
        { 
            File2 << x_posR2[j] << "," << y_posR2[j] << ",";
        }
        File2 << x_posR2[x_posR2.size() - 1] << "," << y_posR2[x_posR2.size() - 1] << endl;

        for (j = 0; j < num_particles; j++)
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

        cout << i+1 << "/" << timesteps << endl;

    }

    std::ofstream File4("dts_LJ.csv");
    for (i = 0; i < timesteps - 1; i++)
    {
        File4 << dts[i] << ",";
    }
    File4 << dts[timesteps - 1] << endl;

    File4.close();


    File1.close();
    File2.close();

    return 0;

}