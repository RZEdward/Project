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

    double x1,x2,y1,y2,distance;
    double k = 1.0;
    double u = 1.0;
    double G = 1.0;
    double boxlims;
    double origin = 0.0;
    int timesteps = 3000;
    int num_particles = 300;
    int numb_particles = num_particles;

    if (num_particles % 2 == 1)
    {
        numb_particles = num_particles + 2;
    }

    double R1 = 1.0;
    double R2 = 1.4;
    double r, D;
    int i, j, n, p, m;
    double x_particle_n, y_particle_n, x_n, y_n, x_j, y_j;
    double dy, dx, theta;
    double dt = 0.4;
    double t;
    int runs = 12;
    double phi;
    int numR1, numR2;

    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::default_random_engine generator;

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
    
    std::ofstream File_rms_v("rms_v.csv");
    std::ofstream File_packingfractions("packingfractions.csv");

    for (m = 1; m < timesteps; m++)
        {   
            int iteration = m;
            File_rms_v << iteration << ", ";
        }
        File_rms_v << timesteps << endl;

    for (m = 0; m < runs; m++)
    {

        phi = 0.78 + m*0.02;
        if (m != runs - 1)
        {
            File_packingfractions << phi << ",";
        }
        else
        {
            File_packingfractions << phi << endl;
        }
        boxlims = sqrt(PI*(numb_particles/2*R1*R1 + num_particles/2*R2*R2)/phi);

        std::uniform_real_distribution<double> uniform(0.0,boxlims);

        for (i = 0; i < num_particles; i++)
            {
                x_pos[i] = uniform(generator);
                y_pos[i] = uniform(generator);
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

        numR1 = 2*x_posR1.size();
        numR2 = 2*x_posR2.size();

        for (i = 0; i < timesteps; i++)
        {

            std::vector<double> shift_x(num_particles); 
            std::vector<double> shift_y(num_particles);
            std::vector<double> x_rand(num_particles);
            std::vector<double> y_rand(num_particles);
            std::vector<double> F_x(num_particles);
            std::vector<double> F_y(num_particles);

            double v2_total = 0.0;
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
                    
                    if (r > D) 
                    {
                        continue;
                    }

                    dx = x_j - x_n;
                    dy = y_j - y_n;
                    theta = atan2(dy,dx);
                    /* 
                    
                    V(r) = GR**3((1-r)/D)**2 (D - r) - if r < D

                    D = R_j + R_n
                    r = |R_j - R_n| = l
                    G = constant = particle modulus = 1
                    R = R_j

                    */

                    F_x[j] += cos(theta)*2*G*R1*R1*R1*(D - r)/D;          
                    F_y[j] += sin(theta)*2*G*R1*R1*R1*(D - r)/D; 
                }

                shift_x[j] = (1/u)*F_x[j]*dt;
                shift_y[j] = (1/u)*F_y[j]*dt;

                double v_x = (1/u)*F_x[j];
                double v_y = (1/u)*F_y[j];

                v2_total += v_x*v_x + v_y*v_y;
            
            }

            double v2_average = v2_total/num_particles;
            double rms_v = sqrt(v2_average); //this is the rms_v at this timestep at this packing fraction
            if (i != timesteps - 1)
            {
                File_rms_v << rms_v << ",";
            }
            else
            {
                File_rms_v << rms_v << endl;
            }

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

        }

    }

    std::ofstream File_info("info.csv");
    File_info << "timesteps" << "," << "dt" << "," << "runs" << "," << "num_particles" << endl;
    File_info << timesteps << "," << dt << "," << runs << "," << num_particles << endl;
    File_info.close();

    File_rms_v.close();

    return 0;

}