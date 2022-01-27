L/*************************************************************************/
/* JUP6.C:  A program to calculate the trajectory of a spacecraft which  */
/*          is encountering the planet Jupiter, in the presense of the   */
/*          moon IO.  This version reads input data from a textfile      */
/*          specified on the command line.                               */
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

double INTER,dt;

/************************************************/
/*            M A I N   P R O G R A M           */
/************************************************/

FILE *out,*in;

main (argc,argv)
   int argc;
   char *argv[];
{
   VECTOR SCacc,SCvel,SCpos,temp,temp2,IOpos,
          IO(),ACC();
   double time = 0.0,l,a,e,R,E,phi,maxtime,
          energy(),L(),alpha(),epsilon(),Rmin(),
          sq();
   int i = 1,j = 1,skip,skip2,logic;

   if (argc >= 2)
      in = fopen(argv[1],"r");
   else
      in = fopen("STDIN","r");
   if (argc = 3)
      out = fopen(argv[2],"w");
   else
      out = fopen("GALILEO.DAT","a");
   fscanf (in,"%lf",&INTER);
   dt = INTER;
   fscanf (in,"%d", &skip);
   fscanf (in,"%d", &skip2);
   fscanf (in,"%lf %lf",&SCpos.x,&SCpos.y);
   SCpos.r = sqrt(sq(SCpos.x) + sq(SCpos.y));
   SCpos.t = atan2(SCpos.y,SCpos.x);
   fscanf (in,"%lf %lf",&SCvel.x,&SCvel.y);
   SCvel.r = sqrt(sq(SCvel.x) + sq(SCvel.y));
   SCvel.t = atan2(SCvel.y,SCvel.x);
   fscanf (in,"%d",&logic);
   fscanf (in,"%lf", &phi);
   phi = phi * PI / 180.0;
   fscanf (in,"%lf",&maxtime);
   IOpos = IO(time,phi);
   dump(0,time,SCpos,SCvel);
   dump(1,time,SCpos,IOpos);
   E = energy (SCpos,SCvel);
   l = L (SCpos,SCvel);
   a = alpha (l);
   e = epsilon (E,l);
   R = Rmin (a,e);
   printf("     %12.5e %12.5e %12.5e %12.5e %12.5e\n",E,l,a,e,R);
   fprintf(out,"%.8e,%.8e,%.8e,%.8f,%.8e",E,l,a,e,R);
   while (time <= maxtime)
   {
      SCacc = ACC(SCpos,IOpos,logic);
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
         fprintf(out,"%.8e,%.8e,%.8e,%.8f,%.8e",E,l,a,e,R);
      }  else j++;
   }
   fclose(out);
}

/************************************************/

VECTOR ACC (location,locIO,on)
   VECTOR location,locIO;
   int on;
{
   VECTOR toJUP,toIO,posIO,total;
   toJUP.r = - pulljup / sq(location.r);
   if (location.r < radJUP) printf("\n\nCRASH...CRASHED ON JUPITER!\n\n");
   toJUP.x = toJUP.r * location.x / location.r;
   toJUP.y = toJUP.r * location.y / location.r;
   if (on != 0)
   {
      posIO.x = location.x - locIO.x;
      posIO.y = location.y - locIO.y;
      posIO.r = sqrt(sq(posIO.x) + sq(posIO.y));
      if (posIO.r < radIO) printf("\n\nCRASH...CRASHED ON IO!\n\n");
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
      fprintf(out,"\n%.0f,%.8e,%.8e,%.8e,%.8e,",l,m.x,m.y,n.x,n.y);
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
