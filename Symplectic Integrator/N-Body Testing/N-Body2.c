// Written by: Nishant Mishra, Research Intern @ PPPL

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

double G = 10000;

double rand_gen() 
{
   // return a uniformly distributed random value
   return ( (double)(rand()) + 1. )/( (double)(RAND_MAX) + 1. );
}

double normalRandom() 
{
   // return a normally distributed random value
   double v1=rand_gen();
   double v2=rand_gen();
   return cos(2*3.14*v2)*sqrt(-2.*log(v1));
}



double radius(double coord1x, double coord1y, double coord1z, double coord2x, double coord2y, double coord2z)
{
    double first_term = pow((coord2x - coord1x), 2);
    double second_term = pow((coord2y - coord1y), 2);
    double third_term = pow((coord2z - coord1z), 2);

    double final_term = pow((first_term + second_term + third_term), (1.00/2.00));
    
    return final_term; 
} 

double force(double mass, int index, int num_bodies, double* radii, double* masses, double* coord) 
{
    double net_Force = 0;

    int z = 0;
    for (int a = 1; a <= num_bodies; a++)
    {   
        for (int b = a + 1; b <= num_bodies; b++)
        {

            if (a == index + 1)
            {
                net_Force += (-G * mass * masses[b-1])/pow(radii[z],2) * (coord[a-1]-coord[b-1])/(radii[z]);
            }

            if (b == index + 1)
            {
                net_Force += (-G * mass * masses[a-1])/pow(radii[z],2) * (coord[b-1]-coord[a-1])/(radii[z]);
            }

            z++;
        }
    }
    return net_Force;
}

double energy(double mass, int index, int num_bodies, 
              double x_velocity, double y_velocity, double z_velocity, 
              double* masses, double* radii) 
{
    double kinetic_energy = (1.00 / 2.00) * mass * (pow(x_velocity,2) + pow(y_velocity,2) + pow(z_velocity,2));
    double potential_energy = 0;

    int c = 0;
    for (int d = 1; d <= num_bodies; d++)
    {   
        for (int e = d + 1; e <= num_bodies; e++)
        {
            
            if (d == index + 1 || e == index + 1)
            {
                potential_energy += (-G * masses[d-1] * masses[e-1])/radii[c];
            }
            c++;
        }
    }
    return (kinetic_energy + potential_energy);
}

//CHEKCS IF EXCEED LIM CONDITION IS MET - 1 is limit is exceeded, 0 if not
int exceedLim(double* x_A_array, double* y_A_array, double* z_A_array, 
              double* x_array,   double* y_array,   double* z_array, 
              int num_bodies,    double lim,        int max_min)
{
    double x_diff, y_diff, z_diff;

    if (max_min == 1)
    {
        for (int v = 0; v < num_bodies; v++)
        {
            x_diff = fabs(x_A_array[v] - x_array[v]);
            y_diff = fabs(y_A_array[v] - y_array[v]);
            z_diff = fabs(z_A_array[v] - z_array[v]);

            if (x_diff > lim || y_diff > lim || z_diff > lim)
            {
                return 1;
            }
        }
        return 0;
    }

    else 
    {
        for (int y = 0; y < num_bodies; y++)
        {
            x_diff = fabs(x_A_array[y] - x_array[y]);
            y_diff = fabs(y_A_array[y] - y_array[y]);
            z_diff = fabs(z_A_array[y] - z_array[y]);

            if (x_diff > lim || y_diff > lim || z_diff > lim)
            {
                return 0;
            }
        }
        return 1;
    }
}

//RETURNS A  RANDOM INTEGER 
double returnRandoms(int lower, int upper) 
{ 
    double num = (rand() % (upper - lower + 1)) + lower; 
    return num; 
} 

int integrate(double time_elapsed, double dt, int num_bodies, double min_lim, double max_lim, int skip_save, 
              double* masses, double* init_x_pos, double* init_y_pos, double* init_z_pos, double* init_x_veloc, 
              double* init_y_veloc, double* init_z_veloc) 
{

    int num_radii = (num_bodies * (num_bodies-1))/2;

    //PHYSICAL PARAMTER INITIATION
    double x_pos[num_bodies];
    double y_pos[num_bodies];
    double z_pos[num_bodies];

    double x_veloc[num_bodies];
    double y_veloc[num_bodies];
    double z_veloc[num_bodies];

    double energies[num_bodies];

    double radii[num_radii];
    double ET = 0; 
    double time = 0;

    //FILE POINTER ARRAY INITILIAZATION
    FILE* write_x[num_bodies];
    FILE* write_y[num_bodies];
    FILE* write_z[num_bodies];

    FILE* append_x[num_bodies];
    FILE* append_y[num_bodies];
    FILE* append_z[num_bodies];

    FILE* write_E[num_bodies];
    FILE* append_E[num_bodies];

    int fileCount;

    //STRING ARRAY INITIALIZITION
    char x_string[50];
    char y_string[50];
    char z_string[50];
    char E_string[50];

    char hold[50];

    //printf("Before: %d\n", num_bodies);

    //STRING GENERATOR + FILE POINTER ALLOCATION
    for (int w = 0; w < num_bodies; w++)
    {
        fileCount = w+1;

        itoa(fileCount, hold, 2); 

        sprintf(x_string, "x%d.txt", fileCount);
        sprintf(y_string, "y%d.txt", fileCount);
        sprintf(z_string, "z%d.txt", fileCount);

        sprintf(E_string, "E%d.txt", fileCount);

        write_x[w] = fopen(x_string, "w");
        write_y[w] = fopen(y_string, "w");
        write_z[w] = fopen(z_string, "w");

        write_E[w] = fopen(E_string, "w");

        append_x[w] = fopen(x_string, "a");
        append_y[w] = fopen(y_string, "a");
        append_z[w] = fopen(z_string, "a");

        append_E[w] = fopen(E_string, "a");
    }

    FILE* ETFILE = fopen("ET.txt", "w");
    FILE* timeFILE = fopen("time.txt", "w");

    FILE* ETFILEA = fopen("ET.txt", "a");
    FILE* timeFILEA = fopen("time.txt", "a");

    //INITIAL PARAMETER ASSIGNMENT/FILE SAVING
    for (int f = 0; f < num_bodies; f++)
    {
        x_pos[f] = init_x_pos[f];
        y_pos[f] = init_y_pos[f];
        z_pos[f] = init_z_pos[f];
        
        fprintf(write_x[f], "%Lf\n", x_pos[f]);
        fprintf(write_y[f], "%Lf\n", y_pos[f]);
        fprintf(write_z[f], "%Lf\n", z_pos[f]);

        fclose(write_x[f]);
        fclose(write_y[f]);
        fclose(write_z[f]);
        
        x_veloc[f] = init_x_veloc[f];
        y_veloc[f] = init_y_veloc[f];
        z_veloc[f] = init_z_veloc[f];
    }

    //INITIAL RADII ARRAY ASSIGNMENT
    int count = 0; 
    int starter = 1;
    for (int g = 0; g < num_bodies; g++) 
    {          
        for (int h = starter; h < num_bodies; h++) 
        {

            radii[count] = radius(x_pos[g], y_pos[g], z_pos[g], x_pos[h], y_pos[h], z_pos[h]);
            count++;
        }
        starter++;
    }
    count = 0;
  
    //INITIAL ENERGY ARRAY ASSIGNMENT
    ET = 0;
    for (int i = 0; i < num_bodies; i++)
    {
        energies[i] = energy(masses[i], i, num_bodies, x_veloc[i], y_veloc[i], z_veloc[i], masses, radii);
        ET += energies[i];

        fprintf(write_E[i], "%Lf\n", energies[i]);
        fclose(write_E[i]);
    }
        
    fprintf(ETFILE, "%Lf\n", ET);
    fprintf(timeFILE, "%Lf\n", time);   
    fclose(ETFILE);
    fclose(timeFILE);    

    //HALF, HALFH, HALFD ARRAY INTIALIZATION FOR ADAPTIVE TIME-STEPPING 
    double x_half_pos[num_bodies];
    double y_half_pos[num_bodies];
    double z_half_pos[num_bodies];

    double x_halfH_pos[num_bodies];
    double y_halfH_pos[num_bodies];
    double z_halfH_pos[num_bodies];

    double x_halfD_pos[num_bodies];
    double y_halfD_pos[num_bodies];
    double z_halfD_pos[num_bodies];
    
    double radii_half[num_radii];

    int loop_count = 1; 

    int adp_Time;

    //MAIN INTEGRATION LOOP
    do
    {   
        loop_count++;
        //HALF, HALFH, HALFD COMPUTATION
        for (int k = 0; k < num_bodies; k++)
        {
            x_half_pos[k] = x_veloc[k] * (dt/2.00) + x_pos[k];
            y_half_pos[k] = y_veloc[k] * (dt/2.00) + y_pos[k];
            z_half_pos[k] = z_veloc[k] * (dt/2.00) + z_pos[k];

            x_halfH_pos[k] = x_veloc[k] * (dt/4.00) + x_pos[k];
            y_halfH_pos[k] = y_veloc[k] * (dt/4.00) + y_pos[k];
            z_halfH_pos[k] = z_veloc[k] * (dt/4.00) + z_pos[k]; 

            x_halfD_pos[k] = x_veloc[k] * (dt) + x_pos[k];
            y_halfD_pos[k] = y_veloc[k] * (dt) + y_pos[k];
            z_halfD_pos[k] = z_veloc[k] * (dt) + z_pos[k]; 
        }  

        //ADAPTIVE TIME STEPPING MAX LIM CALCULATION
        adp_Time = exceedLim(x_halfH_pos, y_halfH_pos, z_halfH_pos, x_half_pos, y_half_pos, z_half_pos, num_bodies, max_lim, 1);
        while (adp_Time == 1)
        {
            if (dt < .000005)
            {
                break;
            }
            
            dt = dt / 2.00;
           
            for (int l = 0; l < num_bodies; l++)
            {
                x_half_pos[l] = x_veloc[l] * (dt/2.00) + x_pos[l];
                y_half_pos[l] = y_veloc[l] * (dt/2.00) + y_pos[l];
                z_half_pos[l] = z_veloc[l] * (dt/2.00) + z_pos[l];

                x_halfH_pos[l] = x_veloc[l] * (dt/4.00) + x_pos[l];
                y_halfH_pos[l] = y_veloc[l] * (dt/4.00) + y_pos[l];
                z_halfH_pos[l] = z_veloc[l] * (dt/4.00) + z_pos[l]; 
            }
            adp_Time = exceedLim(x_halfH_pos, y_halfH_pos, z_halfH_pos, x_half_pos, y_half_pos, z_half_pos, num_bodies, max_lim, 1);
        }
        
        //ADAPTIVE TIME STEPPING MIN LIM CALCULATION
        adp_Time = exceedLim(x_halfD_pos, y_halfD_pos, z_halfD_pos, x_half_pos, y_half_pos, z_half_pos, num_bodies, min_lim, 0);
        while (adp_Time == 1)
        {
            if (dt > 0.5)
            {
                break;
            }

            dt = 2 * dt;

            for (int m = 0; m < num_bodies; m++)
            {
                x_half_pos[m] = x_veloc[m] * (dt/2.00) + x_pos[m];
                y_half_pos[m] = y_veloc[m] * (dt/2.00) + y_pos[m];
                z_half_pos[m] = z_veloc[m] * (dt/2.00) + z_pos[m];

                x_halfD_pos[m] = x_veloc[m] * (dt) + x_pos[m];
                y_halfD_pos[m] = y_veloc[m] * (dt) + y_pos[m];
                z_halfD_pos[m] = z_veloc[m] * (dt) + z_pos[m]; 
            }
            adp_Time = exceedLim(x_halfD_pos, y_halfD_pos, z_halfD_pos, x_half_pos, y_half_pos, z_half_pos, num_bodies, max_lim, 0);
        }

        //RADII HALF LOOP COMPUTATION
        starter = 1;
        for (int n = 0; n < num_bodies; n++) 
        {
        
            for (int o = starter; o < num_bodies; o++) 
            {
                radii_half[count] = radius(x_half_pos[n], y_half_pos[n], z_half_pos[n], x_half_pos[o], y_half_pos[o], z_half_pos[o]);
                count++;
            }
            starter++;
        }
        count = 0;

        //LEAPFROG VELOCITY+POSITION COMPUTATION
        for (int p = 0; p < num_bodies; p++)
        {
            x_veloc[p] = x_veloc[p] + dt * force(masses[p], p, num_bodies, radii_half, masses, x_half_pos)/masses[p];
            y_veloc[p] = y_veloc[p] + dt * force(masses[p], p, num_bodies, radii_half, masses, y_half_pos)/masses[p];
            z_veloc[p] = z_veloc[p] + dt * force(masses[p], p, num_bodies, radii_half, masses, z_half_pos)/masses[p];

            x_pos[p] = x_half_pos[p] + (dt / 2.00) * x_veloc[p];
            y_pos[p] = y_half_pos[p] + (dt / 2.00) * y_veloc[p];
            z_pos[p] = z_half_pos[p] + (dt / 2.00) * z_veloc[p];

            if (loop_count % skip_save == 0)
            {
                fprintf(append_x[p], "%Lf\n", x_pos[p]);
                fprintf(append_y[p], "%Lf\n", y_pos[p]);
                fprintf(append_z[p], "%Lf\n", z_pos[p]);                
            }
        }
        
        //RADII COMPUTATION
        starter = 1;
        for (int q = 0; q < num_bodies; q++) 
        {
            for (int r = starter; r < num_bodies; r++) 
            {
                radii[count] = radius(x_pos[q], y_pos[q], z_pos[q], x_pos[r], y_pos[r], z_pos[r]);
                count++;
            }
            starter++;
        }
        count = 0;

        //ENERGY COMPUTATION
        ET = 0;
        for (int s = 0; s < num_bodies; s++)
        {
            energies[s] = energy(masses[s], s, num_bodies, x_veloc[s], y_veloc[s], z_veloc[s], masses, radii);
            ET += energies[s];        
            
            if (loop_count % skip_save == 0) 
            {
                fprintf(append_E[s], "%Lf\n", energies[s]);
            }
        }

        if (loop_count % skip_save == 0)
        {
            fprintf(ETFILEA, "%Lf\n", ET);     
            fprintf(timeFILEA,"%Lf\n", time);
        }
        
        time += dt; 
    }
    while (time < time_elapsed); 

    //FILE POINTER CLOSER
    for (int ww = 0; ww < num_bodies; ww++)
    {
        fclose(append_x[ww]); 
        fclose(append_y[ww]);
        fclose(append_z[ww]);

        fclose(append_E[ww]);
    }

    fclose(ETFILEA);
    fclose(timeFILEA);

    return 1; 
}

int main(double time_elapsed, double dt, double min_lim, double max_lim, int num_bodies, int skip_save, int max_file, double sigma, double Mi) 
{

    printf("Integration Code (C) Initiation: Success\n\n");
    printf("%d\n", num_bodies);
    double masses[num_bodies];
    
    double init_x_pos[num_bodies];
    double init_y_pos[num_bodies];
    double init_z_pos[num_bodies];
    
    double init_x_veloc[num_bodies];
    double init_y_veloc[num_bodies];
    double init_z_veloc[num_bodies];
    
    double coord_array[3];

    //RANDOM PARAMETER ASSIGNMENT TO EACH BODY IN SYSTEM
    for (int counter = 0; counter < num_bodies; counter++)
    {
        masses[counter] = normalRandom()*sigma+Mi;//returnRandoms(1, 50);

        init_x_pos[counter] = normalRandom()*sigma+Mi; //returnRandoms(1, 50);
        init_y_pos[counter] = normalRandom()*sigma+Mi; //returnRandoms(1, 50);
        init_z_pos[counter] = normalRandom()*sigma+Mi; //returnRandoms(1, 50);

        init_x_veloc[counter] = 1; //normalRandom()*sigma+Mi; //returnRandoms(1, 50);
        init_y_veloc[counter] = 1; //normalRandom()*sigma+Mi; //returnRandoms(1, 50);
        init_z_veloc[counter] = 1; //normalRandom()*sigma+Mi; //returnRandoms(1, 50);
    }

    _setmaxstdio(max_file);

    integrate(time_elapsed, dt, num_bodies, min_lim, max_lim, skip_save, masses, 
              init_x_pos, init_y_pos, init_z_pos, init_x_veloc, init_y_veloc, init_z_veloc);
    
    printf("Integration Code (C) Termination: Success\n\n");
    
    return 1;
}