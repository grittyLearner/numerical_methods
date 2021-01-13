/*============================================================


Name: yprime.h

Author: Louis Horner

Date: January 11, 2021

Description:
	Here we define the class Yprime whose member functions
are the following numerical methods to approximate first-order
differential equations.
 - Euler's method
 - Improved Euler's method
 - 4th-order Runge-Kutta method
 - Adams-Bashforth-Moulton method


============================================================*/
#include <cmath>
const long double e = 2.71828182845904523536;
const long double pi = 3.14159265358979323846;


class Yprime
{
protected:
	// function pointer yprime wich is f(x, y) = dy/dx
	long double (*yprime)(const long double&, const long double&);
	// function pointer Y wich is the exact solution 
	long double (*Y)(const long double&);
	// gx is x_n, and gy is y_n
	long double gx;
	long double gy;
	// X is the x-coordinate of the point we wish to approximate
	long double X;
public:
	/* constructors:  Yprime (yprime, [Y], x, y, X)
	where
	x and y are the initial conditions
	h is the step size */
	Yprime (long double (*f)(const long double&, const long double&), const long double& cx, const long double& cy, const long double& cX) : yprime(f), gx(cx), gy(cy), X(cX) {}
	Yprime (long double (*f)(const long double&, const long double&), long double (*cY)(const long double&), const long double& cx, const long double& cy, const long double& cX) : yprime(f), Y(cY), gx(cx), gy(cy), X(cX) {}

	// Eulers method
	void eulers (const long double& h)
	{
		std::cout << "Euler's method:\n";
		long double lx = gx;
		long double ly = gy;
		for (;X > lx;) 
		{   
			ly += h*(*yprime)(lx, ly); 
			lx += h;
			std::cout << "(" << lx << ", " << ly << "); ";
			std::cout << "absolute error: " << ( ((*Y)(lx) - ly > 0 )? ((*Y)(lx) - ly) : -1*((*Y)(lx) - ly)) << "\n";
		}   
		std::cout << "\n";
	}
	
	// Improved Eulers method
	void ieulers (long double h)
	{
		std::cout << "Improved Euler's method:\n";
		long double lx = gx;
		long double ly = gy;
		long double ys = gy;
		for (;X > lx;)
		{
			ys = ly + h*(*yprime)(lx, ly);
			ly += (((*yprime)(lx, ly) + (*yprime)( (lx+h), ys))/2)*h;
			lx += h;
			std::cout << "(" << lx << ", " << ly << "); ";
			std::cout << "absolute error: " << ( ((*Y)(lx) - ly > 0 )? ((*Y)(lx) - ly) : -1*((*Y)(lx) - ly)) << "\n";
		}
		std::cout << "\n";
	}
	// Fourth order Runge-Kutta method
	void rk4 (long double h)
	{
		std::cout << "Fourth Order Runge-Kutta method:\n";
		long double lx = gx;
		long double ly = gy;
		for (;X > lx;)
		{
			long double k1 = (*yprime)(lx, ly);
			long double k2 = (*yprime)(lx + h/2, ly + (h/2)*k1);
			long double k3 = (*yprime)(lx + h/2, ly + (h/2)*k2);
			long double k4 = (*yprime)(lx + h, ly + h*k3);

			ly += (h/6)*(k1 + 2*k2 + 2*k3 + k4);
			lx += h;
			std::cout << "(" << lx << ", " << ly << "); ";
			std::cout << "absolute error: " << ( ((*Y)(lx) - ly > 0 )? ((*Y)(lx) - ly) : -1*((*Y)(lx) - ly)) << "\n";
		}
		std::cout << "\n";
	}
	// rk4, one step only
protected:
	long double RK4o (long double h)
	{
		long double k1 = (*yprime)(gx, gy);
		long double k2 = (*yprime)(gx + h/2, gy + (h/2)*k1);
		long double k3 = (*yprime)(gx + h/2, gy + (h/2)*k2);
		long double k4 = (*yprime)(gx + h, gy + h*k3);
		gy += (h/6)*(k1 + 2*k2 + 2*k3 + k4);
		gx += h;

		return gy;
	}
public:
	// Adams–Bashforth–Moulton method
	void abm (long double h)
	{
		long double YY [4];
		long double * y = YY;
		// y_n values: [n-3, n-2, n-1, n]

		long double YP [4];
		long double * yp = YP;
		// y_n' values: [n-3', n-2', n-1', n']

		// Using the 4th order Runge-Kutta methods to find y1,y2,y3
		*(y) = gy;
		*(yp) = (*yprime)(gx, gy);

		for (int i=1; i<=3; i++)
		{
			*(y+i) = this->RK4o(h);
			*(yp+i) = (*yprime)(gx, gy);
{			// printing
			std::cout << "x_" << i << ": " << gx << "\n";
			std::cout << "y_" << i << ": " << *y << "\n";
			std::cout << "y_" << i << ": " << gy << "\n";
			std::cout << "y_" << i << "': " << *yp << "\n";
			std::cout << "absolute error: " << ( ((*Y)(gx) - gy > 0 )? ((*Y)(gx) - gy) : -1*((*Y)(gx) - gy)) << "\n\n";
}
		}

		// X-(h/2) is so that the next step doesn't run
		for (int i=4; gx < (X-(h/2)) ;i++)
		{
			// ys is y_{n+1}^*
			long double ys = *(y+3) + (h/24)*( 55*(*(yp+3)) - 59*(*(yp+2)) + 37*(*(yp+1)) - 9*(*(yp)) );

{			// printing
			std::cout << "y_" << i+1 << "^*: " << ys << "\n";
}
			gy += (h/24)*( 9*((*yprime)(gx+h, ys)) + 19*(*(yp+3)) - 5*(*(yp+2)) + (*(yp+1)) );

			// shifting y values to clear space for y_{n+1}
			*y = *(y+1);
			*(y+1) = *(y+2);
			*(y+2) = *(y+3);

			// shifting yp values to clear space for yp_{n+1}
			*yp = *(yp+1);
			*(yp+1) = *(yp+2);
			*(yp+2) = *(yp+3);

			*(y+3) = gy;
			gx += h;
			*(yp+3) = (*yprime)(gx, *(y+3));

{			// printing
			std::cout << "y_" << i << ": " << *(y+3) << "\n";
			std::cout << "y_" << i << ": " << gy << "\n";
			std::cout << "x_" << i << ": " << gx << "\n";
			std::cout << "absolute error: " << ( ((*Y)(gx) - gy > 0 )? ((*Y)(gx) - gy) : -1*((*Y)(gx) - gy)) << "\n";
			std::cout << "X: " << X << "\n\n";
}
		}
	}
};

// below are the functions where dy/dx = f(x, y) for the differential equation being approximated
long double f (const long double& x, const long double& y) { return  (x + y - 1); }
long double g (const long double& x, const long double& y) { return  std::pow(e, -1*y); }
long double w (const long double& x, const long double& y) { return  2*x*y; }

/* exact solutions */
long double exY (const long double& x) { return std::pow(e, x*x - 1); }
long double eY (const long double& x) { return (std::pow(e, x) - x); }
