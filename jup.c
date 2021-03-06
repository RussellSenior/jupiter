/*
**     J U P I T E R - I O - P R O B E   S I M U L A T I O N
**
**     Written by:  Russell S. Senior
**     Date:  09 December 1987
**
**     This program solve the trajectory of a probe through a 
**     gravitational environment including a stationary planet
**     and a moon moving in a circular orbit.  The method used
**     is the Runge-Kutta.  The step size will be variable.
**
**     Version 11.0
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct { double x,y,r,t; } VECTOR;

#define  PI       3.141592653589793
#define  PERIOD   152859.0
#define  RADIUS   4.218070627e8
#define  PULLJUP  1.2675868e17
#define  PULLIO   5.2477160e12
#define  RADJUP   7.137e7
#define  RADIO    1.726e6
#define  MASS     1.138e3
#define  sq(x)    ((x)*(x))

VECTOR iopos0;    /* position of io at start of interval    */
VECTOR iopos1;    /* position of io at mid-interval         */
VECTOR iopos2;    /* position of io at end of interval      */
VECTOR p0;        /* position of probe at start of interval */
VECTOR p1;        /* predicted position at mid-interval     */
VECTOR p2;        /* predicted position at end of interval  */
VECTOR da;        /* change in acceleration over test int.  */
VECTOR a;         /* average acceleration over interval     */
VECTOR a0;        /* acceleration at p0                     */
VECTOR a1;        /* acceleration at first p1               */
VECTOR a2;        /* acceleration at second p1              */
VECTOR a3;        /* acceleration at p2                     */
VECTOR v;         /* velocity of probe at start of interval */

VECTOR io (double timep);
VECTOR acc (VECTOR *p,VECTOR *i,int useio);

double dt;        /* time interval for calculations         */
double dtest;     /* test time interval                     */
double deviation; /* fraction of allowable change in acc.   */
double timep;     /* elapsed time in simulation             */
double fint;      /* time interval between dumps to file    */
double sint;      /* time interval between dumps to screen  */
double fdump;     /* time at which file dump occurs         */
double sdump;     /* time at which screen dump occurs       */
double phi;       /* offset angle for io                    */
double maxtime;   /* time at which simulation stops         */
double e;         /* total energy of probe                  */
double l;         /* angular momentum of probe              */
double al;        /* alpha - parameter of the trajectory    */
double ep;        /* epsilon - parameter of the trajectory  */
double r;         /* predicted closest approach to jupiter  */

double energy (VECTOR *pos,VECTOR *vel,VECTOR *io,int useio);
double angular (VECTOR *,VECTOR *);

int    iothere;   /* boolean to indicate presence of io     */

char *probe_time(double);

FILE *in,*out;
 
int
main (argc,argv)
    int argc;
    char *argv[];
{
    if (argc != 3)
    {
        fprintf(stderr,"\nUsage: %s infile outfile\n",argv[0]);
        exit(0);
    }
    in = fopen(argv[1],"r");
    out = fopen(argv[2],"w");

    fscanf(in,"%lf %lf %lf",&deviation,&fint,&sint);
    fscanf(in,"%lf %lf %lf %lf",&p0.x,&p0.y,&v.x,&v.y);
    fscanf(in,"%d %lf %lf",&iothere,&phi,&maxtime);
    phi *= PI / 180.0;

    p0.r = sqrt(sq(p0.x) + sq(p0.y));
    iopos0 = io(0.0);
    dtest = 1.0;

    fdump = 0.0;
    sdump = 0.0;

    while (timep < maxtime)
    {
        a0 = acc (&p0,&iopos0,iothere);
        a0.r = sqrt(sq(a0.x) + sq(a0.y));
        p1.x = p0.x + v.x*dtest + a0.x*sq(dtest)/2.0;
        p1.y = p0.y + v.y*dtest + a0.y*sq(dtest)/2.0;
        iopos1 = io (timep + dtest);
        a1 = acc (&p1,&iopos1,iothere);
        a1.r = sqrt(sq(a1.x) + sq(a1.y));
        da.x = a1.x - a0.x;
        da.y = a1.y - a0.y;
        da.r = sqrt(sq(da.x) + sq(da.y));
        dt = dtest * deviation * (a0.r + a1.r) / (2 * da.r);
        if (timep >= fdump || timep >= sdump)
        {
            p0.r = sqrt(sq(p0.x) + sq(p0.y));
            e = energy(&p0,&v,&iopos0,iothere);
            l = angular(&p0,&v);
            al = sq(l) / (sq(MASS) * PULLJUP);
            ep = sqrt (1 + (2.0 * e * sq(l/(MASS*PULLJUP)) / MASS));
            r = al / (1 + ep);

            if (timep >= fdump)
            {
                fprintf(out,
                    "\n%.0f,%18.10e,%18.10e,%18.10e,%18.10e,%18.10e,%18.10e,"
                    "%18.10e,%18.10e,%18.10e,%18.10e,%18.10e",
                    timep,p0.x,p0.y,v.x,v.y,iopos0.x,iopos0.y,e,l,al,ep,r);
                fdump += fint;
            }
            if (timep >= sdump)
            {
                fprintf(stdout,
                    "%s %12.5e %12.5e %12.5e %12.5e\n"
                    "     %12.5e %12.5e %12.5e %12.5e %12.5e\n",
                    probe_time(timep),p0.x,p0.y,v.x,v.y,e,l,al,ep,r);
                sdump += sint;
            }
        }
        iopos1 = io (timep + dt/2.0);
        iopos2 = io (timep + dt);
        p1.x = p0.x + v.x*dt/2.0 + a0.x*sq(dt/2.0)/2.0;
        p1.y = p0.y + v.y*dt/2.0 + a0.y*sq(dt/2.0)/2.0;
        a1 = acc (&p1,&iopos1,iothere);
        p1.x = p0.x + v.x*dt/2.0 + a1.x*sq(dt/2.0)/2.0;
        p1.y = p0.y + v.y*dt/2.0 + a1.y*sq(dt/2.0)/2.0;
        a2 = acc (&p1,&iopos1,iothere);
        p2.x = p0.x + v.x*dt + a2.x*sq(dt)/2.0;
        p2.y = p0.y + v.y*dt + a2.y*sq(dt)/2.0;
        a3 = acc (&p2,&iopos2,iothere);

        a.x = (a0.x + 2.0*a1.x + 2.0*a2.x + a3.x) / 6.0;
        a.y = (a0.y + 2.0*a1.y + 2.0*a2.y + a3.y) / 6.0;
        p0.x += v.x*dt + a.x*sq(dt)/2.0;
        p0.y += v.y*dt + a.y*sq(dt)/2.0;
        v.x += a.x*dt;
        v.y += a.y*dt;
        iopos0.x = iopos2.x;
        iopos0.y = iopos2.y;
        timep += dt;
    }
    fclose (in);
    fclose (out);
}

/*************************************/

VECTOR acc (p,i,useio)
    VECTOR *p;
    VECTOR *i;
    int useio;
{
    VECTOR total;
    VECTOR ptoi;
    
    p->r = sqrt(sq(p->x) + sq(p->y));
    total.r = - PULLJUP / sq(p->r);
    total.x = total.r * (p->x) / (p->r);
    total.y = total.r * (p->y) / (p->r);
    if (useio)
    {
        ptoi.x = (p->x) - (i->x);
        ptoi.y = (p->y) - (i->y);
        ptoi.r = sqrt(sq(ptoi.x) + sq(ptoi.y));
        total.r = -PULLIO / sq(ptoi.r);
        total.x += total.r * ptoi.x / ptoi.r;
        total.y += total.r * ptoi.y / ptoi.r;
    }
    return total;
}

/*************************************/

VECTOR io (timep)
    double timep;
{
    VECTOR pos;

    pos.r = RADIUS;
    pos.t = (timep * 2.0 * PI / PERIOD) + phi;
    pos.x = pos.r * cos(pos.t);
    pos.y = pos.r * sin(pos.t);
    return pos;
}

/*************************************/

double energy(pos,vel,io,useio)
    VECTOR *pos;
    VECTOR *vel; 
    VECTOR *io;
    int useio;
{
    double distance,total;
    distance = sqrt(sq((pos->x)-(io->x)) + sq((pos->y)-(io->y)));
    (vel->r) = sqrt(sq(vel->x) + sq(vel->y));
    total = (0.5 * sq(vel->r)) - (PULLJUP / (pos->r));
    if (useio)
        total -= (PULLIO / distance);
    return (total * MASS);
}

/*************************************/

double angular(pos,vel)
    VECTOR *pos;
    VECTOR *vel;
{
    double a,b;
    a = (pos->x) / (pos->r);
    b = (pos->y) / (pos->r);
    return ((pos->r) * MASS * ((a * (vel->y)) - (b * (vel->x))));
}   
 
/*************************************/

char *probe_time (seconds)
    double seconds;
{
    static char buffer[20];
    int hrs=0,min=0;

    if (seconds < 0.0)
    {
        seconds *= -1;
        buffer[0] = '-';
    }
    else
        buffer[0] = ' ';
    while (seconds >= 3600.0)
    {
        seconds -= 3600.0;
        hrs++;
    }
    while (seconds >= 60.0)
    {
        seconds -= 60.0;
        min++;
    }
    sprintf(buffer+1,"%03d:%02d:%05.2f",hrs,min,seconds);
    return buffer;
}
