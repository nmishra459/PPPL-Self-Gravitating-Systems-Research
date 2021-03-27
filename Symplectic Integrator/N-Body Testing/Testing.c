
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

double G = 100;

double force(double mass, int index, int num_bodies, double* radii, double* masses, double* coord) {
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

double energy(double mass, int index, int num_bodies, double x_velocity, double y_velocity, double z_velocity, double* masses, double* radii) {
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

double radius(double coord1x, double coord1y, double coord1z, double coord2x, double coord2y, double coord2z)
{
    double first_term = pow((coord2x - coord1x), 2);
    double second_term = pow((coord2y - coord1y), 2);
    double third_term = pow((coord2z - coord1z), 2);

    double final_term = pow((first_term + second_term + third_term), (1.00/2.00));
    
    return final_term; 
} 

int exceedLim(double* x_A_array, double* z_A_array, double* y_A_array, double* x_array, double* y_array, double* z_array, int num_bodies, double lim, int max_min)
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

int integrate(double time_elapsed, int num_bodies, double* masses, double* init_x_pos, double* init_y_pos, double* init_z_pos, double* init_x_veloc, double* init_y_veloc, double* init_z_veloc) {

    int num_radii = (num_bodies * (num_bodies-1))/2;

    double x_pos[num_bodies];
    double y_pos[num_bodies];
    double z_pos[num_bodies];

    double x_veloc[num_bodies];
    double y_veloc[num_bodies];
    double z_veloc[num_bodies];

    double energies[num_bodies];

    double radii[num_radii];
    double ET, time;

    FILE* write_x[num_bodies];
    FILE* write_y[num_bodies];
    FILE* write_z[num_bodies];

    FILE* append_x[num_bodies];
    FILE* append_y[num_bodies];
    FILE* append_z[num_bodies];

    FILE* write_E[num_bodies];
    FILE* append_E[num_bodies];

    int fileCount;

    char x_string[50];
    char y_string[50];
    char z_string[50];
    char E_string[50];

    char hold[50];

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

        append_x[w] = fopen(x_string, "a");
        append_y[w] = fopen(y_string, "a");
        append_z[w] = fopen(z_string, "a");

        write_E[w] = fopen(E_string, "w");
        append_E[w] = fopen(E_string, "a");
    }

    FILE* ETFILE = fopen("ET.txt", "w");
    FILE* ETFILEA = fopen("ET.txt", "a");
    FILE* timeFILE = fopen("time.txt", "w");
    FILE* timeFILEA = fopen("time.txt", "a");

    for (int f = 0; f < num_bodies; f++)
    {
        x_pos[f] = init_x_pos[f];
        y_pos[f] = init_y_pos[f];
        z_pos[f] = init_z_pos[f];

        x_veloc[f] = init_x_veloc[f];
        y_veloc[f] = init_y_veloc[f];
        z_veloc[f] = init_z_veloc[f];
    }

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
    
    for (int i = 0; i < num_bodies; i++)
    {
        energies[i] = energy(masses[i], i, num_bodies, x_veloc[i], y_veloc[i], z_veloc[i], masses, radii);
    }
    
    for (int j = 0; j < num_bodies; j++) 
    {

        ET += energies[j];
    }

    for (int x = 0; x < num_bodies; x++)
    {
        fprintf(write_x[x], "%lf\n", x_pos[x]);
        fprintf(write_y[x], "%lf\n", y_pos[x]);
        fprintf(write_z[x], "%lf\n", z_pos[x]);

        fprintf(write_E[x], "%lf\n", energies[x]);
    }

    fprintf(ETFILE, "%lf\n", ET);
    fprintf(timeFILE, "%lf\n", time);
    
    time = 0;
    ET = 0;
    int i = 0;
    double totalTime = time_elapsed;
    double dt = 0.00001;              //initial time step
    double minLim = 0.0000000000005;  //lower bound for distance traveled
    double maxLim = 0.001;            //upper bound for distance traveled

    double x_half_pos[num_bodies];
    double y_half_pos[num_bodies];
    double z_half_pos[num_bodies];

    double x_halfH_pos[num_bodies];
    double y_halfH_pos[num_bodies];
    double z_halfH_pos[num_bodies];

    double x_halfD_pos[num_bodies];
    double y_halfD_pos[num_bodies];
    double z_halfD_pos[num_bodies];
    
    double radii_half[num_bodies];

    int adp_Time;
    int loop_count = 1;
    int divisor = 10000;
  
    do
    {   
        loop_count++;
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

        // ADAPTIVE TIME STEPPING
        //////////////////////////////////////////////////////////////////////////////////////
        adp_Time = exceedLim(x_halfH_pos, y_halfH_pos, z_halfH_pos, x_half_pos, y_half_pos, z_half_pos, num_bodies, maxLim, 1);
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
            adp_Time = exceedLim(x_halfH_pos, y_halfH_pos, z_halfH_pos, x_half_pos, y_half_pos, z_half_pos, num_bodies, maxLim, 1);
        }
        
        adp_Time = exceedLim(x_halfD_pos, y_halfD_pos, z_halfD_pos, x_half_pos, y_half_pos, z_half_pos, num_bodies, minLim, 0);
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
            adp_Time = exceedLim(x_halfD_pos, y_halfD_pos, z_halfD_pos, x_half_pos, y_half_pos, z_half_pos, num_bodies, maxLim, 0);
        }
        //////////////////////////////////////////////////////////////////////////////////////

        //LEAPFROG ALGORITHIM COMPUTATION
        
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

        for (int p = 0; p < num_bodies; p++)
        {
            x_veloc[p] = x_veloc[p] + dt * force(masses[p], p, num_bodies, radii_half, masses, x_half_pos)/masses[p];
            y_veloc[p] = y_veloc[p] + dt * force(masses[p], p, num_bodies, radii_half, masses, y_half_pos)/masses[p];
            z_veloc[p] = z_veloc[p] + dt * force(masses[p], p, num_bodies, radii_half, masses, z_half_pos)/masses[p];

            x_pos[p] = x_half_pos[p] + (dt / 2.00) * x_veloc[p];
            y_pos[p] = y_half_pos[p] + (dt / 2.00) * y_veloc[p];
            z_pos[p] = z_half_pos[p] + (dt / 2.00) * z_veloc[p];
        }
        
        if (loop_count % divisor == 0) {
            for (int xx = 0; xx < num_bodies; xx++)
            {
                fprintf(append_x[xx], "%lf\n", x_pos[xx]);
                fprintf(append_y[xx], "%lf\n", y_pos[xx]);
                fprintf(append_z[xx], "%lf\n", z_pos[xx]);
            }
        }

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

        for (int s = 0; s < num_bodies; s++)
        {
            energies[s] = energy(masses[s], s, num_bodies, x_veloc[s], y_veloc[s], z_veloc[s], masses, radii);
        }

        for (int t = 0; t < num_bodies; t++) 
        {
            ET = ET + energies[t];
        }

        //TEXT FILE APPENDING
        if (loop_count % divisor == 0)
        {
            for (int yy = 0; yy < num_bodies; yy++)
            {
                fprintf(append_E[yy], "%lf\n", energies[yy]);
            }
        
            fprintf(ETFILEA, "%lf\n", ET);     
            fprintf(timeFILEA,"%lf\n", time);
        }
        
        time += dt; //Time Increment   
        ET = 0;
    }
    while (time < totalTime); 

    for (int w = 0; w < num_bodies; w++)
    {
        fclose(write_x[w]); 
        fclose(write_y[w]);
        fclose(write_z[w]);

        fclose(append_x[w]); 
        fclose(append_y[w]);
        fclose(append_z[w]);

        fclose(write_E[w]);
        fclose(append_E[w]);
    }
    
    return 1; 
}

int main() 
{
    printf("Integration Code (C) Initiation: Success\n\n");

/*    //BODY 1 PARAMETERS
    double m1 = 1;
    double x10 = 1;
    double y10 = 1;
    double z10 = 0;    
    double vx10 = -0.5;
    double vy10 = -0.5;
    double vz10 = 0;
   
    //BODY 2 PARAMETERS
    double m2 = 1;
    double x20 = 2;
    double y20 = 1;
    double z20 = 0;
    double vx20 = 0.5;
    double vy20 = -0.5;
    double vz20 = 0;

    //BODY 3 PARAMETERS
    double m3 = 1;
    double x30 = 3;
    double y30 = 1;
    double z30 = 0;
    double vx30 = -0.5;
    double vy30 = -0.5;
    double vz30 = 0;

    double m4 = 1;
    double x40 = 4;
    double y40 = 1;
    double z40 = 0;
    double vx40 = 0.5;
    double vy40 = 0.5;
    double vz40 = 0;

    double m5 = 1;
    double x50 = 1;
    double y50 = 2;
    double z50 = 0;
    double vx50 = -0.5;
    double vy50 = 0.5;
    double vz50 = 0;

    double m6 = 1;
    double x60 = 2;
    double y60 = 2;
    double z60 = 0;
    double vx60 = 0.5;
    double vy60 = 0.5;
    double vz60 = 0;

    double m7 = 1;
    double x70 = 3;
    double y70 = 2;
    double z70 = 0;
    double vx70 = -0.5;
    double vy70 = 0.5;
    double vz70 = 0;

    double m8 = 1;
    double x80 = 4;
    double y80 = 2;
    double z80 = 0;
    double vx80 = 0.5;
    double vy80 = 0.5;
    double vz80 = 0;*/

    double time_elapsed = 40;
    int num_bodies = 5;
 
    double masses[num_bodies];
    
    double init_x_pos[num_bodies];
    double init_y_pos[num_bodies];
    double init_z_pos[num_bodies];
    
    double init_x_veloc[num_bodies];
    double init_y_veloc[num_bodies];
    double init_z_veloc[num_bodies];
    
    for (int counter = 0 ; counter < num_bodies; counter++)
    {
        masses[counter] = rand() % 100;

        init_x_pos[counter] = rand() % 100;
        init_y_pos[counter] = rand() % 100;
        init_z_pos[counter] = rand() % 100;

        init_x_veloc[counter] = rand() % 100;
        init_y_veloc[counter] = rand() % 100;
        init_z_veloc[counter] = rand() % 100;
    }

    /*double masses[] =  {m1,m2,m3,m4,m5,m6, m7, m8};
    
    double init_x_pos[] = {x10, x20, x30, x40, x50,x60, x70, x80};
    double init_y_pos[] = {y10, y20, y30, y40, y50,y60, y70, y80};
    double init_z_pos[] = {z10, z20, z30, z40, z50,z60, z70, z80};
    
    double init_x_veloc[] = {vx10, vx20, vx30, vx40, vx50, vx60, vx70, vx80};
    double init_y_veloc[] = {vy10, vy20, vy30, vy40, vy50, vy60, vy70, vy80};
    double init_z_veloc[] = {vz10, vz20, vz30, vz40, vz50, vz60, vz70, vz80};*/

    integrate(time_elapsed, num_bodies, masses, init_x_pos, init_y_pos, init_z_pos, init_x_veloc, init_y_veloc, init_z_veloc);
    
    printf("Integration Code (C) Termination: Success\n\n");
    
    return 1;
}

