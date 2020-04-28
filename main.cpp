#include <stdio.h>
#include <math.h>
#include "winbgi2.h"
#include "rk4.h"
const double g = 10.0; //przyspieszenie ziemskie
const double L = 2.0; //d³ugoœæ promienia pó³kuli
const double n = 0.1; //wspó³czynnik tarcia
void rhs(double x, double y[], double f[])
{
	f[0] = y[1];
	if (f[0] > 0)
	{
		f[1] = (g / L) *(-1 * sin(y[0])) - n / L * (g*cos(y[0]) + L * f[0] * f[0]);
	}
	else if (f[0] < 0)
	{
		f[1] = (g / L) *(-1 * sin(y[0])) + n / L * (g*cos(y[0]) + L * f[0] * f[0]);
	}
}

void main()
{
	FILE*f1 = fopen("Energia.xls", "w");
	FILE*f2 = fopen("Ruch.xls", "w");
	FILE*f3 = fopen("Moment.xls", "w");
	graphics(600, 600);
	scale(0., -3., 8., 3.);
	double t = 0, h = 0.01, tmax = 8.0;
	double y1[2];
	double Ek, Ep, Ec, x, y, x0, y0, Vx, Vy, V0, Ms, T;
	double m = 5;
	int a = 2;
	printf("Podaj wartosci poczatkowe : x0, V0 \n");//Promien pó³kuli jest ustalony, st¹d proœba jedynie o jedna wspolrzedna poczatkowa.
	scanf("%lf %lf", &x0, &V0);
	while (a > 0)
	{
		if ((x0 < -L) || (x0 > L) || ((V0 == 0) && (x0 == 0))) //pkt musi le¿eæ wewn¹trz pó³kuli, oraz w pkt. (0;0) nie moze miec predkosci rownej 0
		{
			printf("Podaj wartosci poczatkowe : x0, V0 \n");
			scanf("%lf %lf", &x0, &V0);
		}
		else
		{
			a = -1;
		}
	}
	y0 = sqrt(L*L - x0 * x0);
	y1[0] = atan((L - y0) / x0);
	y1[1] = V0 / L;
	y1[0] = -1.;
	y1[1] = 0.5;
	Ek = 1 / 2.0*m*y1[1] * y1[1] * L*L;
	Ep = g * m*L - m * g*L*cos(y1[0]);
	Ec = Ek + Ep;
	Vx = cos(y1[0])*y1[1] * L;
	Vy = sin(y1[0])*y1[1] * L;
	Ms = -m * g / L * sin(y1[0]);
	if (y1[1] > 0)
	{
		T = -m * n / L * (g*cos(y1[0]) + L * y1[1] * y1[1]);
	}
	else if (y1[1] < 0)
	{
		T = m * n / L * (g*cos(y1[0]) + L * y1[1] * y1[1]);
	}
	fprintf(f1, "%lf \t %lf \t %lf \t %lf \n ", t, Ek, Ep, Ec);
	fprintf(f2, "%lf \t %lf \t %lf \t %lf \t %lf \t %lf \n ", t, x0, y0, Vx, Vy, y1[0]);
	fprintf(f3, "%lf \t %lf \t %lf \t %lf \t %lf \n ", t, Vx, Vy, T*cos(y1[0]), T*sin(y1[0]));
	for (int i = 0; i < tmax / h; ++i)
	{
		vrk4(t, y1, h, 2, rhs, y1);
		t += h;
		Ek = 1 / 2.0*m*y1[1] * y1[1] * L*L;
		Ep = g * m*L - m * g*L*cos(y1[0]);
		Ec = Ek + Ep;
		x = L * sin(y1[0]);
		y = L - L * cos(y1[0]);
		Vx = cos(y1[0])*y1[1] * L;
		Vy = sin(y1[0])*y1[1] * L;
		printf("%lf \n", y1[0]);
		Ms = m * g / L * (-1 * sin(y1[0]));
		if (y1[1] > 0)
		{
			T = -m * n / L * (g*cos(y1[0]) + L * y1[1] * y1[1]);
		}
		else if (y1[1] < 0)
		{
			T = m * n / L * (g*cos(y1[0]) + L * y1[1] * y1[1]);
		}
		fprintf(f1, "%lf \t %lf \t %lf \t %lf \n ", t, Ek, Ep, Ec);
		fprintf(f2, "%lf \t %lf \t %lf \t %lf \t %lf \t %lf \n ", t, x, y, Vx, Vy, y1[1]);
		fprintf(f3, "%lf \t %lf \t %lf \t %lf \t %lf \n ", t, Vx, Vy, T*cos(y1[0]), T*sin(y1[0]));
		if ((y1[0] > 1.57) || (y1[0] < -1.57))
		{
			break;// warunek na 'wyskoczenie' z pólkuli
		}
		/*setcolor(0.1);
		circle(t,y1[0],1); kat fi
		setcolor(0.99);
		circle(t, y1[1], 1); predkosc katowa*/
		setcolor(0.5);
		circle(t, L*sin(y1[0]), 1);//wspolrzedna x
		setcolor(0.7);
		circle(t, L - L * cos(y1[0]), 1); // wspolrzedna y
	}
	fclose(f1);
	fclose(f2);
	fclose(f3);
	wait();
}
