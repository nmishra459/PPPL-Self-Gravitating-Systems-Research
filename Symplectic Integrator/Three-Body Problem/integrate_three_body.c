
#include <stdio.h>
#include <math.h>
#include <stdlib.h>

double G = 1.00;

double force(double mass1, double mass2, double mass3, double radius1, double radius2, double coord1, double coord2, double coord3) {
    double first_term = (G * mass1 * mass2) / (pow(radius1,2)) * (coord2 - coord1) / (radius1);
    double second_term = (G * mass1 * mass3) / pow(radius2,2) * (coord3 - coord1) / (radius2);
    return (first_term + second_term);
}

double energy(double mass1, double mass2, double mass3, double xvelocity, double yvelocity, double zvelocity, double radius1, double radius2) {
    double first_term = (1.00 / 2.00) * mass1 * (pow(xvelocity,2) + pow(yvelocity,2) + pow(zvelocity,2));
    double second_term = ((-G * mass1 * mass2) / radius1);
    double third_term = ((-G * mass1 * mass3) / radius2);
    return (first_term + second_term + third_term);
}

double radius(double coord1x, double coord1y, double coord1z, double coord2x, double coord2y, double coord2z) {
    double first_term = pow((coord2x - coord1x),2);
    double second_term = pow((coord2y - coord1y),2);
    double third_term = pow((coord2z - coord1z),2);
    double final_term = pow((first_term + second_term + third_term), (1.00/2.00));
    return final_term; 
} 

int integrate(double m1, double m2, double m3, double x10, double y10, double z10, double vx10, double vy10, double vz10, double x20, double y20, double z20, double vx20, double vy20, double vz20, double x30, double y30, double z30, double vx30, double vy30, double vz30) {

    double r12, r23, r13, x1, y1, z1, x2, y2, z2, x3, y3, z3, vx1, vy1, vz1, vx2, vy2, vz2, vx3, vy3, vz3, E1, E2, E3, ET, time;

    //File: Write
    FILE* x1FILE = fopen("x1.txt", "w");
    FILE* y1FILE = fopen("y1.txt", "w");
    FILE* z1FILE = fopen("z1.txt", "w");

    FILE* x2FILE = fopen("x2.txt", "w");
    FILE* y2FILE = fopen("y2.txt", "w");
    FILE* z2FILE = fopen("z2.txt", "w");

    FILE* x3FILE = fopen("x3.txt", "w");
    FILE* y3FILE = fopen("y3.txt", "w");
    FILE* z3FILE = fopen("z3.txt", "w");

    FILE* E1FILE = fopen("E1.txt", "w");
    FILE* E2FILE = fopen("E2.txt", "w");
    FILE* E3FILE = fopen("E3.txt", "w");
    FILE* ETFILE = fopen("ET.txt", "w");

    FILE* timeFILE = fopen("time.txt", "w");

    //File: Append
    FILE* x1FILEA = fopen("x1.txt", "a");
    FILE* y1FILEA = fopen("y1.txt", "a");
    FILE* z1FILEA = fopen("z1.txt", "a");

    FILE* x2FILEA = fopen("x2.txt", "a");
    FILE* y2FILEA = fopen("y2.txt", "a");
    FILE* z2FILEA = fopen("z2.txt", "a");

    FILE* x3FILEA = fopen("x3.txt", "a");
    FILE* y3FILEA = fopen("y3.txt", "a");
    FILE* z3FILEA = fopen("z3.txt", "a");

    FILE* E1FILEA = fopen("E1.txt", "a");
    FILE* E2FILEA = fopen("E2.txt", "a");
    FILE* E3FILEA = fopen("E3.txt", "a");
    FILE* ETFILEA = fopen("ET.txt", "a");

    FILE* timeFILEA = fopen("time.txt", "a");

    x1 = x10;
    y1 = y10;
    z1 = z10;

    x2 = x20;
    y2 = y20;
    z2 = z20;

    x3 = x30;
    y3 = y30;
    z3 = z30;

    r12 = radius(x1, y1, z1, x2, y2, z2);
    r23 = radius(x2, y2, z2, x3, y3, z3);
    r13 = radius(x1, y1, z1, x3, y3, z3);

    vx1 = vx10;
    vy1 = vy10;
    vz1 = vz10;

    vx2 = vx20;
    vy2 = vy20;
    vz2 = vz20;

    vx3 = vx30;
    vy3 = vy30;
    vz3 = vz30;

    E1 = energy(m1, m2, m3, vx1, vy1, vz1, r12, r13);  
    E2 = energy(m2, m1, m3, vx2, vy2, vz2, r12, r23);
    E3 = energy(m3, m1, m2, vx3, vy3, vz3, r13, r23);
    ET = E1 + E2 + E3;

    time = 0;

    fprintf(x1FILE, "%lf\n", x1);
    fprintf(y1FILE, "%lf\n", y1);
    fprintf(z1FILE, "%lf\n", z1);

    fprintf(x2FILE, "%lf\n", x2);
    fprintf(y2FILE, "%lf\n", y2);
    fprintf(z2FILE, "%lf\n", z2);

    fprintf(x3FILE, "%lf\n", x3);
    fprintf(y3FILE, "%lf\n", y3);
    fprintf(z3FILE, "%lf\n", z3);

    fprintf(E1FILE, "%lf\n", E1);
    fprintf(E2FILE, "%lf\n", E2);
    fprintf(E3FILE, "%lf\n", E3);
    fprintf(ETFILE, "%lf\n", ET);

    fprintf(timeFILE, "%lf\n", time);

    int i = 0;
    double totalTime = 3;
    double dt = 0.00001;  //initial time step
    double minLim = 0.0000000000005;  //lower bound for distance traveled
    double maxLim = 0.001;  //upper bound for distance traveled

    double x1half, y1half, z1half;
    double x2half, y2half, z2half; 
    double x3half, y3half, z3half;
    double x1halfH, y1halfH, z1halfH;
    double x2halfH, y2halfH, z2halfH;
    double x3halfH, y3halfH, z3halfH;
    double x1halfD, y1halfD, z1halfD;
    double x2halfD, y2halfD, z2halfD;
    double x3halfD, y3halfD, z3halfD;
    
    do
    {   

        x1half = vx1 * (dt / 2.00) + x1;
        y1half = vy1 * (dt / 2.00) + y1;
        z1half = vz1 * (dt / 2.00) + z1;

        x2half = vx2 * (dt / 2.00) + x2;
        y2half = vy2 * (dt / 2.00) + y2;
        z2half = vz2 * (dt / 2.00) + z2;

        x3half = vx3 * (dt / 2.00) + x3;
        y3half = vy3 * (dt / 2.00) + y3;
        z3half = vz3 * (dt / 2.00) + z3;

        x1halfH = vx1 * (dt / 4.00) + x1;
        y1halfH = vy1 * (dt / 4.00) + y1;
        z1halfH = vz1 * (dt / 4.00) + z1;

        x2halfH = vx2 * (dt / 4.00) + x2;
        y2halfH = vy2 * (dt / 4.00) + y2;
        z2halfH = vz2 * (dt / 4.00) + z2;

        x3halfH = vx3 * (dt / 4.00) + x3;
        y3halfH = vy3 * (dt / 4.00) + y3;
        z3halfH = vz3 * (dt / 4.00) + z3;

        x1halfD = vx1 * (dt) + x1;
        y1halfD = vy1 * (dt) + y1;
        z1halfD = vz1 * (dt) + z1;

        x2halfD = vx2 * (dt) + x2;
        y2halfD = vy2 * (dt) + y2;
        z2halfD = vz2 * (dt) + z2;

        x3halfD = vx3 * (dt) + x3;
        y3halfD = vy3 * (dt) + y3;
        z3halfD = vz3 * (dt) + z3;


        // ADAPTIVE TIME STEPPING
        //////////////////////////////////////////////////////////////////////////////////////
        while (fabs((x1halfH - x1half)) > maxLim || fabs((x2halfH - x2half)) > maxLim
               || fabs((x3halfH - x3half)) > maxLim || fabs((y1halfH - y1half)) > maxLim
               || fabs((y2halfH - y2half)) > maxLim || fabs((y3halfH - y3half)) > maxLim
               || fabs((z1halfH - z1half)) > maxLim || fabs((z2halfH - z2half)) > maxLim
               || fabs((z3halfH - z3half)) > maxLim)
        {
        
            if (dt < .000005){
                break;
            }
            
            dt = dt / 2.00;
           
            x1half = vx1 * (dt / 2.00) + x1;
            y1half = vy1 * (dt / 2.00) + y1;
            z1half = vz1 * (dt / 2.00) + z1;

            x2half = vx2 * (dt / 2.00) + x2;
            y2half = vy2 * (dt / 2.00) + y2;
            z2half = vz2 * (dt / 2.00) + z2;

            x3half = vx3 * (dt / 2.00) + x3;
            y3half = vy3 * (dt / 2.00) + y3;
            z3half = vz3 * (dt / 2.00) + z3;
            
            x1halfH = vx1 * (dt / 4.00) + x1;
            y1halfH = vy1 * (dt / 4.00) + y1;
            z1halfH = vz1 * (dt / 4.00) + z1;

            x2halfH = vx2 * (dt / 4.00) + x2;
            y2halfH = vy2 * (dt / 4.00) + y2;
            z2halfH = vz2 * (dt / 4.00) + z2;

            x3halfH = vx3 * (dt / 4.00) + x3;
            y3halfH = vy3 * (dt / 4.00) + y3;
            z3halfH = vz3 * (dt / 4.00) + z3;


        }
        
        while (fabs((x1halfD - x1half)) < minLim && fabs((x2halfD - x2half)) < minLim
               && fabs((x3halfD - x3half)) < minLim && fabs((y1halfD - y1half)) < minLim
               && fabs((y2halfD - y2half)) < minLim && fabs((y3halfD - y3half)) < minLim
               && fabs((z1halfD - z1half)) < minLim && fabs((z2halfD - z2half)) < minLim
               && fabs((z3halfD - z3half)) < minLim) 
        {

            if (dt > 0.5) {
                break;
            }

            dt = 2 * dt;

            x1half = vx1 * (dt / 2.00) + x1;
            y1half = vy1 * (dt / 2.00) + y1;
            z1half = vz1 * (dt / 2.00) + z1;

            x2half = vx2 * (dt / 2.00) + x2;
            y2half = vy2 * (dt / 2.00) + y2;
            z2half = vz2 * (dt / 2.00) + z2;

            x3half = vx3 * (dt / 2.00) + x3;
            y3half = vy3 * (dt / 2.00) + y3;
            z3half = vz3 * (dt / 2.00) + z3;

            x1halfD = vx1 * (dt) + x1;
            y1halfD = vy1 * (dt) + y1;
            z1halfD = vz1 * (dt) + z1;

            x2halfD = vx2 * (dt) + x2;
            y2halfD = vy2 * (dt) + y2;
            z2halfD = vz2 * (dt) + z2;

            x3halfD = vx3 * (dt) + x3;
            y3halfD = vy3 * (dt) + y3;
            z3halfD = vz3 * (dt) + z3;


        }
        //////////////////////////////////////////////////////////////////////////////////////

        //LEAPFROG ALGORITHIM COMPUTATION
        double rhalf12 = radius(x1half, y1half, z1half, x2half, y2half, z2half);
        double rhalf23 = radius(x2half, y2half, z2half, x3half, y3half, z3half);
        double rhalf13 = radius(x1half, y1half, z1half, x3half, y3half, z3half);

        vx1 = vx1 + dt * force(m1, m2, m3, rhalf12, rhalf13, x1half, x2half, x3half) / m1;
        vx2 = vx2 + dt * force(m2, m1, m3, rhalf12, rhalf23, x2half, x1half, x3half) / m2;
        vx3 = vx3 + dt * force(m3, m1, m2, rhalf13, rhalf23, x3half, x1half, x2half) / m3;

        vy1 = vy1 + dt * force(m1, m2, m3, rhalf12, rhalf13, y1half, y2half, y3half) / m1;
        vy2 = vy2 + dt * force(m2, m1, m3, rhalf12, rhalf23, y2half, y1half, y3half) / m2;
        vy3 = vy3 + dt * force(m3, m1, m2, rhalf13, rhalf23, y3half, y1half, y2half) / m3;

        vz1 = vz1 + dt * force(m1, m2, m3, rhalf12, rhalf13, z1half, z2half, z3half) / m1;
        vz2 = vz2 + dt * force(m2, m1, m3, rhalf12, rhalf23, z2half, z1half, z3half) / m2;
        vz3 = vz3 + dt * force(m3, m1, m2, rhalf13, rhalf23, z3half, z1half, z2half) / m3;

        x1 = x1half + dt / 2.00 * vx1;
        y1 = y1half + dt / 2.00 * vy1;
        z1 = z1half + dt / 2.00 * vz1;

        x2 = x2half + dt / 2.00 * vx2;
        y2 = y2half + dt / 2.00 * vy2;
        z2 = z2half + dt / 2.00 * vz2;

        x3 = x3half + dt / 2.00 * vx3;
        y3 = y3half + dt / 2.00 * vy3;
        z3 = z3half + dt / 2.00 * vz3;

        r12 = radius(x1, y1, z1, x2, y2, z2);
        r23 = radius(x2, y2, z2, x3, y3, z3);
        r13 = radius(x1, y1, z1, x3, y3, z3);
        
        E1 = energy(m1, m2, m3, vx1, vy1, vz1, r12, r13);
        E2 = energy(m2, m1, m3, vx2, vy2, vz2, r12, r23);
        E3 = energy(m3, m1, m2, vx3, vy3, vz3, r13, r23);
        ET = E1 + E2 + E3;

        //TEXT FILE APPENDING
        fprintf(x1FILEA, "%lf\n", x1);
        fprintf(y1FILEA, "%lf\n", y1);
        fprintf(z1FILEA, "%lf\n", z1);

        fprintf(x2FILEA, "%lf\n", x2);
        fprintf(y2FILEA, "%lf\n", y2);
        fprintf(z2FILEA, "%lf\n", z2);

        fprintf(x3FILEA, "%lf\n", x3);
        fprintf(y3FILEA, "%lf\n", y3);
        fprintf(z3FILEA, "%lf\n", z3);

        fprintf(E1FILEA, "%lf\n", E1);
        fprintf(E2FILEA, "%lf\n", E2);
        fprintf(E3FILEA, "%lf\n", E3);
        fprintf(ETFILEA, "%lf\n", ET);     

        fprintf(timeFILEA,"%lf\n", time);

        time += dt; //Time Increment   

    }
    while (time < totalTime);

    //FILE CLOSING
    fclose(x1FILE);
    fclose(y1FILE);
    fclose(z1FILE);
    
    fclose(x2FILE);
    fclose(y2FILE);
    fclose(z2FILE);
    
    fclose(x3FILE);
    fclose(y3FILE);
    fclose(z3FILE);

    fclose(E1FILE);
    fclose(E2FILE);
    fclose(E3FILE);
    fclose(ETFILE);

    fclose(timeFILE);

    fclose(x1FILEA);
    fclose(y1FILEA);
    fclose(z1FILEA);
    
    fclose(x2FILEA);
    fclose(y2FILEA);
    fclose(z2FILEA);
    
    fclose(x3FILEA);
    fclose(y3FILEA);
    fclose(z3FILEA);
    
    fclose(E1FILEA);
    fclose(E2FILEA);
    fclose(E3FILEA);
    fclose(ETFILEA);
       
    fclose(timeFILEA);   

    return 1; 
}

int main() {
    printf("Integration Code (C) Initiation: Success\n\n");

    double m1 = 1;
    double x10 = -0.97000436;
    double y10 = 0.24308753;
    double z10 = 0;    
    double vx10 = -0.5*0.93240737;
    double vy10 = -0.5*0.86473146;
    double vz10 = 0;
   
    //BODY 2 PARAMETERS
    double m2 = 1;
    double x20 = 0.97000436;
    double y20 = -0.24308753;
    double z20 = 0;
    double vx20 = -0.5*0.93240737;
    double vy20 = -0.5*0.86473146;
    double vz20 = 0;

    //BODY 3 PARAMETERS
    double m3 = 1;
    double x30 = 0;
    double y30 = 0;
    double z30 = 0;
    double vx30 = 0.93240737;
    double vy30 = 0.86473146;
    double vz30 = 0;


    integrate(m1, m2, m3, x10, y10, z10, vx10, vy10, vz10, x20, y20, z20, vx20, vy20, vz20, x30, y30, z30, vx30, vy30, vz30);
    
    printf("Integration Code (C) Termination: Success\n\n");
    
    return 1;
}
