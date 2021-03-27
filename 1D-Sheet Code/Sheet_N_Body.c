#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.14159265
#define G 1.00

double drand(double low, double high);
double force(int index, double sheet_count, double* coord);
double radius(double coord1x, double coord2x);
int integrate(double* x_pos_0, double* x_veloc_0, double sheet_count, double sheet_mass, double time_elapsed, int skip_save);
int compareDouble(const void *pa, const void *pb);


typedef struct 
{
    double position;
    double velocity;
}
Particle;


int main(void)
{
    // initial Physical Constants
    double sheet_count = 100.00;
    int max_file = 8000;
    double dt = 0.05;
    double time_elapsed = 10.00;//0.000079;
    int skip_save = 1;
    double sheet_mass = 1.00/sheet_count;

    double x_pos_0[(int)sheet_count];
    double x_veloc_0[(int)sheet_count];

    printf("Integration Code (C) Initiatiated \n\n");

    //RANDOM PARAMETER ASSIGNMENT TO EACH BODY IN SYSTEM
    for (int i = 0; i < sheet_count; i++)
    {
        x_pos_0[i] = drand(-1, 1);
        x_veloc_0[i] = drand(-1.6, 1.6);
    }

    integrate(x_pos_0, x_veloc_0, sheet_count, sheet_mass, time_elapsed, skip_save);

    printf("Integration Code (C) Termination: Success\n\n");

    return 1;
}

int comparePosition(const void *pa, const void *pb)
{
    Particle *p1 = (Particle *) pa;
    Particle *p2 = (Particle *) pb;

    if (p1->position < p2->position) return -1;
    if (p1->position > p2->position) return 1;    
    return 0;
}

double drand(double low, double high)
{
    return ((double)rand() * (high - low) / (double)RAND_MAX) + low;
}


double quad_solver(double a, double b, double c) {
    double discriminant, root1, root2, realPart, imagPart;

    discriminant = b * b - 4.00 * a * c;

    // condition for real and different roots
    if (discriminant > 0.00) 
    {
        //printf("%lf\n", a);
        root1 = (-b + sqrt(discriminant)) / (2.00 * a);
        root2 = (-b - sqrt(discriminant)) / (2.00 * a);
        if (root1 > 0 && root2 > 0)
        {
            if(root1 < root2)
            {
                return root1;
            }
            else 
            {
                return root2;
            }
        }
        else if (root1 > 0)
        {
            return root1;
        }
        else if (root2 > 0)
        {
            return root2;
        }
    }
    // condition for real and equal roots
    else if (discriminant == 0) 
    {
        return -b / (2.00 * a);
    }
    // if roots are not real
    else 
    {
        return -1.00;
    }

    return -1.00;
} 

int integrate(double* x_pos_0, double* x_veloc_0, double sheet_count, double sheet_mass, double time_elapsed, int skip_save)
{
    int num_sheets = (int)sheet_count;
    //PHYSICAL PARAMTER INITIATION
    Particle particles[num_sheets];

    double total_time = 0.00;
    double time = -1.00;
    double x_half_pos[num_sheets];

    //INITIAL PARAMETER ASSIGNMENT/FILE SAVING
    for (int i = 0; i < num_sheets; i++)
    {
        particles[i].position = x_pos_0[i];
        particles[i].velocity = x_veloc_0[i];
    }

    qsort(particles, num_sheets, sizeof(Particle), comparePosition);
    
    int loop_count = 1;
    int flag = -1;
    int new_flag = -1;
    double accel_diff = 0.5 * (1.00/sheet_mass);
    double dt = 0.5;
    do
    {        
        //printf("%d\n", flag);

        if (flag == num_sheets - 1 || flag == num_sheets - 2)
        {
            dt = -1;
            int counter = 0;
            while (dt < 0)
            {
                dt = quad_solver(accel_diff, particles[counter].velocity - particles[counter+1].velocity, particles[counter].position - particles[counter+1].position);
                counter++;

            }
        }
        else
        {
            dt = -1;
            int counter = flag+1;
            while (dt < 0)
            {
                dt = quad_solver(accel_diff, particles[counter].velocity - particles[counter+1].velocity, particles[counter].position - particles[counter+1].position);
                counter++;

            }            
        }

        for (int i = 0; i < num_sheets - 1; i++)
        {
            if (i != flag && i + 1 != flag)
            {
                time = quad_solver(accel_diff, particles[i].velocity - particles[i+1].velocity, particles[i].position - particles[i+1].position);
                //printf("%lf\n", time);

                if (time < dt && time > 0) 
                {
                    dt = time;
                    new_flag = i;
                    //printf("%d\n", flag);

                }
            }
            else if ((i+1 == flag) && (i+1 != num_sheets - 1)) 
            {
                time = quad_solver(accel_diff, particles[i].velocity - particles[i+2].velocity, particles[i].position - particles[i+2].position);
                if (time < dt && time > 0) 
                {
                    dt = time;
                    new_flag = i;
                }                
            }

        }

        //printf("%d\n", flag);

        for (int i = 0; i < num_sheets; i++)
        {
            double index = (double) i;
            double acceleration = ((sheet_count - 1.00) - 2.00*index)/(sheet_mass);
            particles[i].velocity = particles[i].velocity + dt * (acceleration);
            particles[i].position = particles[i].position + particles[i].velocity * dt + (1.00/2.00) * acceleration * pow(dt, 2.00);
            if (i == new_flag)
            {
                printf("%Lf\n",particles[i].position);
                printf("%Lf\n",particles[i+1].position);

            }


        }

        //SWAPPING - MIGHT NOT WORK?? USE * Maybe
        Particle temp = particles[new_flag];
        particles[new_flag] = particles[new_flag+1];
        particles[new_flag+1] = temp;

        flag = new_flag;
        //printf("%lf\n", total_time);


        total_time += dt;
    }
    while (total_time < time_elapsed);

    // for (int i = 0; i < num_sheets; i++)
    // {
    //    printf("%lf\n", particles[i].position);
    // }

    //FILE POINTER ARRAY INITILIAZATION
    FILE* write_x[num_sheets];
    FILE* write_v[num_sheets];

    int fileCount;

    //STRING ARRAY INITIALIZITION
    char x_string[50];
    char v_string[50];
    char hold[50];

    //STRING GENERATOR + FILE POINTER ALLOCATION
    for (int i = 0; i < num_sheets; i++)
    {
        fileCount = i + 1;
        sprintf(hold, "%d", fileCount);
 
        sprintf(x_string, "x%d.txt", fileCount);
        write_x[i] = fopen(x_string, "w");
        fprintf(write_x[i], "%Lf\n", particles[i].position);
        fclose(write_x[i]);

        sprintf(v_string, "v%d.txt", fileCount);
        write_v[i] = fopen(v_string, "w");
        fprintf(write_v[i], "%Lf\n", particles[i].velocity);
        fclose(write_v[i]);        
    }

    return 1;
}
