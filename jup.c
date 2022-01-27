/*************************************************************************/
/* JUP5.C:  A program to calculate the trajectory of a spacecraft which  */
/*          is encountering the planet Jupiter, in the presense of the   */
/*          moon IO.  This version also includes several bug fixes.      */
/*************************************************************************/

#include <stdio.h>
#include <math.h>

typedef struct {double x,y,r,t;} VECTOR;

#define RADIUS   4.218070627e8
#define PERIOD   152859.0
#define pullio   5.2477160e12
#define pulljup  1.2675868e17
#define radJUP   7.137e7
#define radIO    1.726e6
#define mass     1.138e3

double INTER;

/************************************************/
/*            M A I N   P R O G R A M           */
/************************************************/

FILE *out;

main ()
{
   VECTOR SCacc,SCvel,SCpos,temp,temp2,IOpos,
          IO(),ACC(),INITpos(),INITvel();
   double dt,time = 0.0,l,a,e,R,E,phi,maxtime,
          energy(),L(),alpha(),epsilon(),Rmin(),
          INTERVAL(),IOphase(),sq(),MAXtime(),DT();
   int i = 1,j = 1,skip,skip2,logic;
   out = fopen("GALILEO.DAT","a");
   INTER = INTERVAL();
   skip  = INITskip(1);
   skip2 = INITskip(2);
   SCpos = INITpos();
   SCvel = INITvel();
   logic = IOon();
   phi   = IOphase();
   IOpos = IO(time,phi);
   maxtime = MAXtime();
   dump(0,time,SCpos,SCvel);
   dump(1,time,SCpos,IOpos);
   while (time <= maxtime)
   {
      SCacc = ACC(SCpos,IOpos,logic);
      dt = DT(SCpos,SCvel);
      temp.x = SCpos.x + (SCvel.x * dt) + (SCacc.x * sq(dt));
      temp.y = SCpos.y + (SCvel.y * dt) + (SCacc.y * sq(dt));
      temp.r = sqrt(sq(temp.x) + sq(temp.y));
      temp.t = atan2(temp.y,temp.x);
      IOpos = IO(time+dt,phi);
      temp2 = ACC(temp,IOpos,logic);
      SCacc.x = (SCacc.x + temp2.x) / 2.0;
      SCacc.y = (SCacc.y + temp2.y) / 2.0;
      SCpos.x = SCpos.x + (SCvel.x * dt) + (SCacc.x * sq(dt));
      SCpos.y = SCpos.y + (SCvel.y * dt) + (SCacc.y * sq(dt));
      SCpos.r = sqrt(sq(SCpos.x) + sq(SCpos.y));
      SCpos.t = atan2(SCpos.y,SCpos.x);
      SCvel.x = SCvel.x + (SCacc.x * dt);
      SCvel.y = SCvel.y + (SCacc.y * dt);
      SCvel.r = sqrt(sq(SCvel.x) + sq(SCvel.y));
      SCvel.t = atan2(SCvel.y,SCvel.x);
      time += dt;
      if (i == skip)
      {
         dump(0,time,SCpos,SCvel);
         i = 1;
         E = energy (SCpos,SCvel);
         l = L (SCpos,SCvel);
         a = alpha (l);
         e = epsilon (E,l);
         R = Rmin (a,e);
         printf("     %12.5e %12.5e %12.5e %12.5e %12.5e\n",E,l,a,e,R);
      }  else i++;
      if (j == skip2)
      {
         dump(1,time,SCpos,IOpos);
         j = 1;
         E = energy (SCpos,SCvel);
         l = L (SCpos,SCvel);
         a = alpha (l);
         e = epsilon (E,l);
         R = Rmin (a,e);
         fprintf(out,"%12.5e,%12.5e,%12.5e,%12.5e,%12.5e\n",E,l,a,e,R);
      }  else j++;
   }
   fclose(out);
}

/************************************************/

double INTERVAL ()
{
   double m;
   printf("Please enter the length of standard interval (in seconds) ...");
   scanf ("%lf",&m);
   return m;
}

/************************************************/

INITskip (k)
   short k;
{
   int j;
   if (k == 1)
   {
      printf("Please enter the number of loops between output lines... ");
      scanf ("%d", &j);
   }
   if (k == 2)
   {
      printf("Please enter the number of loops between outfile lines... ");
      scanf ("%d", &j);
   }
   return j;
}

/************************************************/

VECTOR INITpos ()
{
   VECTOR pos;
   printf("Please enter the initial position (x,y) of the SC...");
   scanf ("%lf %lf",&pos.x,&pos.y);
   pos.r = sqrt(sq(pos.x) + sq(pos.y));
   pos.t = atan2(pos.y,pos.x);
   return pos;
}

/************************************************/

VECTOR INITvel ()
{
   VECTOR vel;
   printf("Please enter the initial velocity (x,y) of the SC...");
   scanf ("%lf %lf",&vel.x,&vel.y);
   vel.r = sqrt(sq(vel.x) + sq(vel.y));
   vel.t = atan2(vel.y,vel.x);
   return vel;
}

/************************************************/


VECTOR ACC (location,locIO,on)
   VECTOR location,locIO;
   int on;
{
   VECTOR toJUP,toIO,posIO,total;
   toJUP.r = - pulljup / sq(location.r);
   toJUP.x = toJUP.r * location.x / location.r;
   toJUP.y = toJUP.r * location.y / location.r;
   if (on != 0)
   {
      posIO.x = location.x - locIO.x;
      posIO.y = location.y - locIO.y;
      posIO.r = sqrt(sq(posIO.x) + sq(posIO.y));
      toIO.r  = - pullio  / sq(posIO.r);
      toIO.x  = toIO.r * posIO.x / posIO.r;
      toIO.y  = toIO.r * posIO.y / posIO.r;
      total.x = toJUP.x + toIO.x;
      total.y = toJUP.y + toIO.y;
      total.r = sqrt(sq(total.x) + sq(total.y));
      total.t = atan2(total.y,total.x);
      return total;
   }
   else
   {
      toJUP.t = atan2(toJUP.y,toJUP.x);
      return toJUP;
   }
}

/************************************************/

IOon ()
{
   int logical;
   printf("Shall IO exert a gravitational force (0 for no, else yes)? ... ");
   scanf ("%d",&logical);
   return logical;
}

/************************************************/

double IOphase ()
{
   double phi;
   printf("Please enter the beginning position of IO (in degrees)... ");
   scanf ("%lf", &phi);
   phi = phi * PI / 180.0;
   return phi;
}

/************************************************/

double MAXtime ()
{
   double t;
   printf("Please enter the ending time in seconds... ");
   scanf ("%lf",&t);
   return t;
}

/************************************************/

VECTOR IO (loctime,phi)
   double loctime,phi;
{
   VECTOR pos;
   pos.r = RADIUS;
   pos.t = (loctime * 2.0 * PI / PERIOD) + phi;
   pos.x = pos.r * cos(pos.t);
   pos.y = pos.r * sin(pos.t);
   return pos;
}

/************************************************/

double sq (k)
   double k;
{
   return (k * k);
}

/************************************************/

TIME (second)
   double second;
{
   short hour=0,min=0,logic;
   if (second < 0.0)
   {
      second = - second;
      logic = 1;
      printf("-");
   }
   else
      printf(" ");
   while (second >= 3600.0)
   {  second = second - 3600.0; hour++; }
   while (second >= 60)
   {  second = second - 60.0; min++;    }
   if (hour < 10) printf("0");
   printf("%d:",hour);
   if (min < 10) printf("0");
   printf("%d:",min);
   if (second < 10.0) printf("0");
   printf("%.2f",second);
   return ;
}

/************************************************/

dump (k,l,m,n)
   short k;
   VECTOR m,n;
   double l;
{
   if (k == 0)
   {   TIME (l);
      printf("  %12.5e %12.5e %12.5e %12.5e\n",m.x,m.y,n.x,n.y);
   }
   if (k == 1)
   {
      fprintf(out,"%7.0f,%12.5e,%12.5e,%12.5e,%12.5e,",l,m.x,m.y,n.x,n.y);
   }
   return ;
}

/************************************************/

double energy (pos,vel)
   VECTOR pos,vel;
{
   return (0.5 * mass * sq(vel.r)) - (pulljup * mass / pos.r);
}

/************************************************/

double L (pos,vel)
   VECTOR pos,vel;
{
   double m,n;
   m = pos.x / pos.r;
   n = pos.y / pos.r;
   return (pos.r * mass * ((m * vel.y) - (n * vel.x)));
}

/************************************************/

double alpha (l)
   double l;
{
   return (sq(l) / (sq(mass) * pulljup));
}

/************************************************/

double epsilon (E,l)
   double E,l;
{
   return sqrt(1 + (2.0 * E * sq(l / (mass * pulljup)) / mass));
}

/************************************************/

double Rmin (a,e)
   double a,e;
{
   return (a / (1 + e));
}

/************************************************/

double DT (pos,vel)
   VECTOR pos,vel;
{
   VECTOR temp;
   double temp2,temp3,interval;
   temp.x = pos.x / pos.r;
   temp.y = pos.y / pos.r;
   temp2 = (temp.x * vel.x) + (temp.y * vel.y);
   temp3 = sqrt(sq(vel.r) - sq(temp2));
        if (temp3 < 100.0)  interval = INTER;
   else if (temp3 < 500.0)  interval = INTER / 2.0;
   else if (temp3 < 1000.0) interval = INTER / 4.0;
   else if (temp3 < 2000.0) interval = INTER / 10.0;
   else                     interval = INTER / 20.0;
   return interval;
}
