#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#define PI 3.14159265

double G = 1;

double drand (double low, double high)
{
    return ((double)rand() * (high - low)) / (double)RAND_MAX + low;
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

// ghost_force(masses[p], p, num_bodies, ghost_radii_half, masses, ghost_x_half_pos);            

double ghost_force(double mass, int index, int num_bodies, double* radii, double* masses, double* coord) 
{
    double net_Force = 0;

    int z = (index+1) * num_bodies - 1;
    int b = num_bodies - 1;

    //printf("%d\n::::", z);
    for (int a = 0; a < num_bodies; a++)
    {   
        if (index == a || radii[z-b] == 0)
        {
            net_Force += 0;
        }
        else
        {
            net_Force += (-G * mass * masses[a])/(pow(radii[z-b],2) * (coord[index]-coord[a])/(radii[z-b]));
        }
        b--;
    }
    return net_Force;
}

double kinetic_energy(double mass, double x_velocity, double y_velocity, double z_velocity)
{
    return ((1.00 / 2.00) * mass * (pow(x_velocity,2) + pow(y_velocity,2) + pow(z_velocity,2)));
}

double potential_energy(int num_bodies, int index, double* masses, double* radii)
{
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
    return (potential_energy);    
}

double energy(double mass,       int index,         int num_bodies, 
              double x_velocity, double y_velocity, double z_velocity, 
              double* masses,    double* radii) 
{
    double kin_energy = kinetic_energy(mass, x_velocity, y_velocity, z_velocity);
    double poten_energy = potential_energy(num_bodies, index, masses, radii);
    return (kin_energy + poten_energy);
}

double radius(double coord1x, double coord1y, double coord1z, double coord2x, double coord2y, double coord2z)
{
    double first_term = pow((coord2x - coord1x), 2);
    double second_term = pow((coord2y - coord1y), 2);
    double third_term = pow((coord2z - coord1z), 2);

    double final_term = pow((first_term + second_term + third_term), (1.00/2.00));
    
    return final_term; 
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

int integrate(double time_elapsed, double dt, int num_bodies, double min_lim, double max_lim, int skip_save, 
              double* masses, double* init_x_pos, double* init_y_pos, double* init_z_pos, double* init_x_veloc, 
              double* init_y_veloc, double* init_z_veloc, double* ghost_init_x_pos, double* ghost_init_y_pos, 
              double* ghost_init_z_pos, double* ghost_init_x_veloc, double* ghost_init_y_veloc, double* ghost_init_z_veloc) 
{

    int num_radii = (num_bodies * (num_bodies-1))/2;
    int num_ghost_radii = num_bodies * num_bodies;

    //PHYSICAL PARAMTER INITIATION
    double x_pos[num_bodies];
    double y_pos[num_bodies];
    double z_pos[num_bodies];

    double x_veloc[num_bodies];
    double y_veloc[num_bodies];
    double z_veloc[num_bodies];

    double kinetic_energies[num_bodies];
    double potential_energies[num_bodies];
    double energies[num_bodies];

    double radii[num_radii];

    double ghost_x_pos[num_bodies];
    double ghost_y_pos[num_bodies];
    double ghost_z_pos[num_bodies];

    double ghost_x_veloc[num_bodies];
    double ghost_y_veloc[num_bodies];
    double ghost_z_veloc[num_bodies];

    double ghost_radii[num_ghost_radii];

    double ET = 0; 
    double KET = 0;
    double PET = 0;
    double time = 0;

    //FILE POINTER ARRAY INITILIAZATION
    FILE* write_x[num_bodies];
    FILE* write_y[num_bodies];
    FILE* write_z[num_bodies];

    FILE* append_x[num_bodies];
    FILE* append_y[num_bodies];
    FILE* append_z[num_bodies];

    FILE* write_KE[num_bodies];
    FILE* append_KE[num_bodies];

    FILE* write_PE[num_bodies];
    FILE* append_PE[num_bodies];

    FILE* write_E[num_bodies];
    FILE* append_E[num_bodies];

    FILE* write_ghost_x[num_bodies];
    FILE* write_ghost_y[num_bodies];
    FILE* write_ghost_z[num_bodies];

    FILE* append_ghost_x[num_bodies];
    FILE* append_ghost_y[num_bodies];
    FILE* append_ghost_z[num_bodies];

    int fileCount;

    //STRING ARRAY INITIALIZITION
    char x_string[50];
    char y_string[50];
    char z_string[50];
    char E_string[50];
    char KE_string[50];
    char PE_string[50];

    char ghost_x_string[50];
    char ghost_y_string[50];
    char ghost_z_string[50];

    char hold[50];

    //printf("Before: %d\n", num_bodies);

    //STRING GENERATOR + FILE POINTER ALLOCATION
    for (int w = 0; w < num_bodies; w++)
    {
        fileCount = w+1;

        sprintf(hold, "%d", fileCount); 

        sprintf(x_string, "x%d.txt", fileCount);
        sprintf(y_string, "y%d.txt", fileCount);
        sprintf(z_string, "z%d.txt", fileCount);

        sprintf(E_string, "E%d.txt", fileCount);
        sprintf(E_string, "E%d.txt", fileCount);
        
        sprintf(KE_string, "KE%d.txt", fileCount);
        sprintf(KE_string, "KE%d.txt", fileCount);
       
        sprintf(PE_string, "PE%d.txt", fileCount);        
        sprintf(PE_string, "PE%d.txt", fileCount);

        sprintf(ghost_x_string, "gx%d.txt", fileCount);
        sprintf(ghost_y_string, "gy%d.txt", fileCount);
        sprintf(ghost_z_string, "gz%d.txt", fileCount);

        write_x[w] = fopen(x_string, "w");
        write_y[w] = fopen(y_string, "w");
        write_z[w] = fopen(z_string, "w");

        write_E[w] = fopen(E_string, "w");
        write_KE[w] = fopen(KE_string, "w");
        write_PE[w] = fopen(PE_string, "w");

        append_x[w] = fopen(x_string, "a");
        append_y[w] = fopen(y_string, "a");
        append_z[w] = fopen(z_string, "a");

        append_E[w] = fopen(E_string, "a");
        append_KE[w] = fopen(KE_string, "a");
        append_PE[w] = fopen(PE_string, "a");

        write_ghost_x[w] = fopen(ghost_x_string, "w");
        write_ghost_y[w] = fopen(ghost_y_string, "w");
        write_ghost_z[w] = fopen(ghost_z_string, "w");

        append_ghost_x[w] = fopen(ghost_x_string, "a");
        append_ghost_y[w] = fopen(ghost_y_string, "a");
        append_ghost_z[w] = fopen(ghost_z_string, "a");

    }

    FILE* ETFILE = fopen("ET.txt", "w");
    FILE* KETFILE = fopen("KET.txt", "w");
    FILE* PETFILE = fopen("PET.txt", "w");
    FILE* timeFILE = fopen("time.txt", "w");

    FILE* ETFILEA = fopen("ET.txt", "a");
    FILE* KETFILEA = fopen("KET.txt", "a");
    FILE* PETFILEA = fopen("PET.txt", "a");
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

        ghost_x_pos[f] = ghost_init_x_pos[f];
        ghost_y_pos[f] = ghost_init_y_pos[f];
        ghost_z_pos[f] = ghost_init_z_pos[f];
        
        fprintf(write_ghost_x[f], "%Lf\n", ghost_x_pos[f]);
        fprintf(write_ghost_y[f], "%Lf\n", ghost_y_pos[f]);
        fprintf(write_ghost_z[f], "%Lf\n", ghost_z_pos[f]);

        fclose(write_ghost_x[f]);
        fclose(write_ghost_y[f]);
        fclose(write_ghost_z[f]);
        
        ghost_x_veloc[f] = ghost_init_x_veloc[f];
        ghost_y_veloc[f] = ghost_init_y_veloc[f];
        ghost_z_veloc[f] = ghost_init_z_veloc[f];
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
  
    //INITIAL GHOST RADII ARRAY ASSIGNMENT
    for (int cz = 0; cz < num_bodies; cz++) 
    {          
        for (int ha = 0; ha < num_bodies; ha++) 
        {
            if (cz == ha)
            {
                ghost_radii[count] = 0;
                count++;
            }
            else
            {
                ghost_radii[count] = radius(ghost_x_pos[cz], ghost_y_pos[cz], ghost_z_pos[cz], x_pos[ha], y_pos[ha], z_pos[ha]);
                count++;  
            }
            
        }
    }
    count = 0;

    //INITIAL ENERGY ARRAY ASSIGNMENT
    ET = 0;
    KET = 0;
    PET = 0;
    for (int i = 0; i < num_bodies; i++)
    {
        energies[i] = energy(masses[i], i, num_bodies, x_veloc[i], y_veloc[i], z_veloc[i], masses, radii);
        kinetic_energies[i] = kinetic_energy(masses[i], x_veloc[i], y_veloc[i], z_veloc[i]);
        potential_energies[i] = potential_energy(num_bodies, i, masses, radii);
        
        ET += energies[i];
        KET += kinetic_energies[i];
        PET += potential_energies[i];

        fprintf(write_E[i], "%Lf\n", energies[i]);
        fclose(write_E[i]);

        fprintf(write_KE[i], "%Lf\n", kinetic_energies[i]);
        fclose(write_KE[i]);

        fprintf(write_PE[i], "%Lf\n", potential_energies[i]);
        fclose(write_PE[i]);
    }
        
    fprintf(ETFILE, "%Lf\n", ET);
    fprintf(KETFILE, "%Lf\n", KET);
    fprintf(PETFILE, "%Lf\n", PET);
    fprintf(timeFILE, "%Lf\n", time);   
 
    fclose(ETFILE);
    fclose(KETFILE);
    fclose(PETFILE);
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

    double ghost_x_half_pos[num_bodies];
    double ghost_y_half_pos[num_bodies];
    double ghost_z_half_pos[num_bodies];

    double ghost_x_halfH_pos[num_bodies];
    double ghost_y_halfH_pos[num_bodies];
    double ghost_z_halfH_pos[num_bodies];

    double ghost_x_halfD_pos[num_bodies];
    double ghost_y_halfD_pos[num_bodies];
    double ghost_z_halfD_pos[num_bodies];
    
    double ghost_radii_half[num_ghost_radii];

    int loop_count = 1; 

    int adp_Time;
    int ghost_adp_Time;
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

            ghost_x_half_pos[k] = ghost_x_veloc[k] * (dt/2.00) + ghost_x_pos[k];
            ghost_y_half_pos[k] = ghost_y_veloc[k] * (dt/2.00) + ghost_y_pos[k];
            ghost_z_half_pos[k] = ghost_z_veloc[k] * (dt/2.00) + ghost_z_pos[k];

            ghost_x_halfH_pos[k] = ghost_x_veloc[k] * (dt/4.00) + ghost_x_pos[k];
            ghost_y_halfH_pos[k] = ghost_y_veloc[k] * (dt/4.00) + ghost_y_pos[k];
            ghost_z_halfH_pos[k] = ghost_z_veloc[k] * (dt/4.00) + ghost_z_pos[k]; 

            ghost_x_halfD_pos[k] = ghost_x_veloc[k] * (dt) + ghost_x_pos[k];
            ghost_y_halfD_pos[k] = ghost_y_veloc[k] * (dt) + ghost_y_pos[k];
            ghost_z_halfD_pos[k] = ghost_z_veloc[k] * (dt) + ghost_z_pos[k];
        }  


        //ADAPTIVE TIME STEPPING MAX LIM CALCULATION
        adp_Time = exceedLim(x_halfH_pos, y_halfH_pos, z_halfH_pos, x_half_pos, y_half_pos, z_half_pos, num_bodies, max_lim, 1);
        ghost_adp_Time = exceedLim(ghost_x_halfH_pos, ghost_y_halfH_pos, ghost_z_halfH_pos, ghost_x_half_pos, ghost_y_half_pos, ghost_z_half_pos, num_bodies, max_lim, 1);

        while (adp_Time == 1 || ghost_adp_Time == 1)
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
            
                ghost_x_half_pos[l] = ghost_x_veloc[l] * (dt/2.00) + ghost_x_pos[l];
                ghost_y_half_pos[l] = ghost_y_veloc[l] * (dt/2.00) + ghost_y_pos[l];
                ghost_z_half_pos[l] = ghost_z_veloc[l] * (dt/2.00) + ghost_z_pos[l];

                ghost_x_halfH_pos[l] = ghost_x_veloc[l] * (dt/4.00) + ghost_x_pos[l];
                ghost_y_halfH_pos[l] = ghost_y_veloc[l] * (dt/4.00) + ghost_y_pos[l];
                ghost_z_halfH_pos[l] = ghost_z_veloc[l] * (dt/4.00) + ghost_z_pos[l];     

            }

            adp_Time = exceedLim(x_halfH_pos, y_halfH_pos, z_halfH_pos, x_half_pos, y_half_pos, z_half_pos, num_bodies, max_lim, 1);
            ghost_adp_Time = exceedLim(ghost_x_halfH_pos, ghost_y_halfH_pos, ghost_z_halfH_pos, ghost_x_half_pos, ghost_y_half_pos, ghost_z_half_pos, num_bodies, max_lim, 1);

        }
        
        //ADAPTIVE TIME STEPPING MIN LIM CALCULATION
        adp_Time = exceedLim(x_halfD_pos, y_halfD_pos, z_halfD_pos, x_half_pos, y_half_pos, z_half_pos, num_bodies, min_lim, 0);
        ghost_adp_Time = exceedLim(ghost_x_halfD_pos, ghost_y_halfD_pos, ghost_z_halfD_pos, ghost_x_half_pos, ghost_y_half_pos, ghost_z_half_pos, num_bodies, min_lim, 0);

        while (adp_Time == 1 && ghost_adp_Time == 1)
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

                ghost_x_half_pos[m] = ghost_x_veloc[m] * (dt/2.00) + ghost_x_pos[m];
                ghost_y_half_pos[m] = ghost_y_veloc[m] * (dt/2.00) + ghost_y_pos[m];
                ghost_z_half_pos[m] = ghost_z_veloc[m] * (dt/2.00) + ghost_z_pos[m];

                ghost_x_halfD_pos[m] = ghost_x_veloc[m] * (dt) + ghost_x_pos[m];
                ghost_y_halfD_pos[m] = ghost_y_veloc[m] * (dt) + ghost_y_pos[m];
                ghost_z_halfD_pos[m] = ghost_z_veloc[m] * (dt) + ghost_z_pos[m]; 
            }
            adp_Time = exceedLim(x_halfD_pos, y_halfD_pos, z_halfD_pos, x_half_pos, y_half_pos, z_half_pos, num_bodies, min_lim, 0);
            ghost_adp_Time = exceedLim(ghost_x_halfD_pos, ghost_y_halfD_pos, ghost_z_halfD_pos, ghost_x_half_pos, ghost_y_half_pos, ghost_z_half_pos, num_bodies, min_lim, 0);


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

        for (int na = 0; na < num_bodies; na++) 
        {          
            for (int oa = 0; oa < num_bodies; oa++) 
            {
                if (na == oa)
                {
                    ghost_radii_half[count] = 0;
                }
                else
                { 
                    ghost_radii_half[count] = radius(ghost_x_half_pos[na], ghost_y_half_pos[na], ghost_z_half_pos[na], x_half_pos[oa], y_half_pos[oa], z_half_pos[oa]);
                }
                count++;
            }
        }
        //printf("%lf:::\n ", ghost_radii_half[1]);

        //printf("%d\n", count);
        count = 0;    

        //LEAPFROG VELOCITY+POSITION COMPUTATION
        for (int p = 0; p < num_bodies; p++)
        {
            x_veloc[p] = x_veloc[p] + dt * force(masses[p], p, num_bodies, radii_half, masses, x_half_pos)/masses[p];
            y_veloc[p] = y_veloc[p] + dt * force(masses[p], p, num_bodies, radii_half, masses, y_half_pos)/masses[p];
            z_veloc[p] = z_veloc[p] + dt * force(masses[p], p, num_bodies, radii_half, masses, z_half_pos)/masses[p];
            //printf("%Lf\n", x_veloc[p]);

            x_pos[p] = x_half_pos[p] + (dt / 2.00) * x_veloc[p];
            y_pos[p] = y_half_pos[p] + (dt / 2.00) * y_veloc[p];
            z_pos[p] = z_half_pos[p] + (dt / 2.00) * z_veloc[p];

            ghost_x_veloc[p] = ghost_x_veloc[p] + dt * ghost_force(masses[p], p, num_bodies, ghost_radii_half, masses, ghost_x_half_pos)/masses[p];            
            //printf("%Lf\n", ghost_x_veloc[p]);

            ghost_y_veloc[p] = ghost_y_veloc[p] + dt * ghost_force(masses[p], p, num_bodies, ghost_radii_half, masses, ghost_y_half_pos)/masses[p];
            ghost_z_veloc[p] = ghost_z_veloc[p] + dt * ghost_force(masses[p], p, num_bodies, ghost_radii_half, masses, ghost_z_half_pos)/masses[p];

            ghost_x_pos[p] = ghost_x_half_pos[p] + (dt / 2.00) * ghost_x_veloc[p];
            ghost_y_pos[p] = ghost_y_half_pos[p] + (dt / 2.00) * ghost_y_veloc[p];
            ghost_z_pos[p] = ghost_z_half_pos[p] + (dt / 2.00) * ghost_z_veloc[p];

            if (loop_count % skip_save == 0)
            {
                fprintf(append_x[p], "%Lf\n", x_pos[p]);  
                fprintf(append_y[p], "%Lf\n", y_pos[p]);
                fprintf(append_z[p], "%Lf\n", z_pos[p]); 

                fprintf(append_ghost_x[p], "%Lf\n", ghost_x_pos[p]);  
                fprintf(append_ghost_y[p], "%Lf\n", ghost_y_pos[p]);
                fprintf(append_ghost_z[p], "%Lf\n", ghost_z_pos[p]);                
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

        starter = 1;
        for (int qa = 0; qa < num_bodies; qa++) 
        {
            for (int ra = starter; ra < num_bodies; ra++) 
            {
                ghost_radii[count] = radius(ghost_x_pos[qa], ghost_y_pos[qa], ghost_z_pos[qa], x_pos[ra], y_pos[ra], z_pos[ra]);
                count++;
            }
            starter++;
        }
        count = 0;

        //ENERGY COMPUTATION
        ET = 0;
        KET = 0;
        PET = 0;
        for (int s = 0; s < num_bodies; s++)
        {
            energies[s] = energy(masses[s], s, num_bodies, x_veloc[s], y_veloc[s], z_veloc[s], masses, radii);
            kinetic_energies[s] = kinetic_energy(masses[s], x_veloc[s], y_veloc[s], z_veloc[s]);
            potential_energies[s] = potential_energy(num_bodies, s, masses, radii);
            
            ET += energies[s];        
            KET += kinetic_energies[s];
            PET += potential_energies[s];    

            if (loop_count % skip_save == 0) 
            {
                fprintf(append_E[s], "%Lf\n", energies[s]);
                fprintf(append_KE[s], "%Lf\n", kinetic_energies[s]);
                fprintf(append_PE[s], "%Lf\n", potential_energies[s]);

            }
        }

        if (loop_count % skip_save == 0)
        {
            fprintf(ETFILEA, "%Lf\n", ET);     
            fprintf(KETFILEA, "%Lf\n", KET);     
            fprintf(PETFILEA, "%Lf\n", PET);     
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

        fclose(append_ghost_x[ww]); 
        fclose(append_ghost_y[ww]);
        fclose(append_ghost_z[ww]);

        fclose(append_E[ww]);
        fclose(append_KE[ww]);
        fclose(append_PE[ww]);

    }

    fclose(ETFILEA);
    fclose(KETFILEA);
    fclose(PETFILEA);

    fclose(timeFILEA);

    return 1; 
}

int main() 
{

    printf("Integration Code (C) Initiation: Success\n\n");

    double time_elapsed = 6500;         //total time passed during simulation
    double dt = 0.00001;                //initial time step
    double min_lim = 0.0000000000005;   //lower bound for distance traveled
    double max_lim = 0.1;               //upper bound for distance traveled
    int num_bodies = 3;                 //number of bodies in the system
    int skip_save = 10000;              //the number of loops that need to pass before saving to a txt file

    int max_file = 8000;                //maximum number of files that can be open at once

    double preturb = 1.0001;

    double masses[num_bodies];
    
    double init_x_pos[num_bodies];
    double init_y_pos[num_bodies];
    double init_z_pos[num_bodies];
    
    double init_x_veloc[num_bodies];
    double init_y_veloc[num_bodies];
    double init_z_veloc[num_bodies];

    double ghost_masses[num_bodies];
    
    double ghost_init_x_pos[num_bodies];
    double ghost_init_y_pos[num_bodies];
    double ghost_init_z_pos[num_bodies];
    
    double ghost_init_x_veloc[num_bodies];
    double ghost_init_y_veloc[num_bodies];
    double ghost_init_z_veloc[num_bodies];
    
    double v, theta, phi, r;

    //RANDOM PARAMETER ASSIGNMENT TO EACH BODY IN SYSTEM
    for (int counter = 0; counter < num_bodies; counter++)
    {
        v = drand(0,1);
        theta = drand(0, 2*PI);
        phi = acos((2*v)-1);
        r = pow(drand(0,1), 1/3);

        masses[counter] = 1;
        ghost_masses[counter] = 1;

        init_x_pos[counter] = 100 * r * sin(phi) * cos(theta);
        init_y_pos[counter] = 100 * r * sin(phi) * sin(theta);
        init_z_pos[counter] = 100 * r * cos(phi); 
       
        init_x_veloc[counter] = drand(-.1,.1);
        init_y_veloc[counter] = drand(-.1,.1);
        init_z_veloc[counter] = drand(-.1,.1);

        ghost_init_x_pos[counter] = init_x_pos[counter];
        ghost_init_y_pos[counter] = init_y_pos[counter];
        ghost_init_z_pos[counter] = init_z_pos[counter] ;
       
        ghost_init_x_veloc[counter] = init_x_veloc[counter] * preturb;
        ghost_init_y_veloc[counter] = init_y_veloc[counter]; //* preturn
        ghost_init_z_veloc[counter] = init_z_veloc[counter]; //* preturb

    }

    _setmaxstdio(max_file);

    integrate(time_elapsed, dt, num_bodies, min_lim, max_lim, skip_save, masses, init_x_pos, init_y_pos, init_z_pos, init_x_veloc, init_y_veloc, init_z_veloc, ghost_init_x_pos, 
        ghost_init_y_pos, ghost_init_z_pos, ghost_init_x_veloc, ghost_init_y_veloc, ghost_init_z_veloc);
    
    printf("Integration Code (C) Termination: Success\n\n");
    
    return 1;
}
