/*============================================================


Name: ode.h

Author: Louis Horner

Date: January 20, 2021

Description:
	Here we define three classes, below, to approximate three
classes of problems using their respective method.

fstOrder
	meant to solve problems of the following form.
		dy/dx = f(x, y)
		y(x_0) = y_0
	Constructor:
		fstOrder (yprime, [Y], x_0, y_0, X)
			X is the value of x at the point we wish to approximate
	Public members:
		eulers(h)
			Euler's method
			h is the step size
			returns y_n
		ieulers(h)
			Improved Euler's method
			h is the step size
			returns y_n
		rk4(h)
			4th order Runge-Kutta method
			h is the step size
			returns y_n
		abm(h)
			Adams–Bashforth–Moulton method
			h is the step size
			returns y_n


sndOrder
	meant to solve problems of the following form.
		d^2y/dx^2 = f(x, y, dy/dx), or equivalently

		y' = u
		u = f(x, y, u)

		y(x_0) = y_0
		y'(x_0) = y_1

	Constructor:
		sndOrder (ydprime, [Y], x, y, u, X)
			x, y, and u are the initial conditions x_0, y_0, y_1
	Public members:
		eulers (h)
			Euler's method
			h is the step size
			returns y_n
		rk4 (h)
			4th order Runge-Kutta method
			h is the step size
			returns y_n


fstOrderSys
	meant to solve problems of the following form.
		dx/dt = f(t, x, y)
		dy/dt = g(t, x, y)

		x(t_0) = x_0
		y(t_0) = y_0

	Constructor:
		fstOrderSys (xprime, yprime, [X], [Y], t_0, x_0, y_0, T) 
			T is the value of t at the point we wish to approximate
	Public members:
		rk4(h)
			4th order Runge-Kutta method
			h is the step size
			returns a pointer to an array: [x_n, y_n, [error in x_n], [error in y_n]]

============================================================*/


class fstOrder
{
protected:
	// function pointer yprime which is f(x, y) = dy/dx
	long double (*yprime)(long double, long double);
	// function pointer Y wich is the exact solution y(x)
	long double (*Y)(long double);
	// gx and gy are x_n and y_n
	long double gx;
	long double gy;
	long double X;

	// value to determine whether an exact solution was given
	bool is_def_Y;
public:
	fstOrder
	(
		long double (*f)(long double, long double), 
		long double cx, 
		long double cy, 
		long double cX
	) : yprime(f), gx(cx), gy(cy), X(cX)
	{
		is_def_Y = 0;
	}
	fstOrder
	(
		long double (*f)(long double, long double), 
		long double (*cY)(long double), 
		long double cx, 
		long double cy, 
		long double cX
	) : yprime(f), Y(cY), gx(cx), gy(cy), X(cX) 
	{
		is_def_Y = 1;
	}

	// Eulers method
	long double eulers (long double h)
	{
		std::cout << "Euler's method applied to a first order differential equation: (h = " << h << ")\n";
		// copies of the variables so that the originals are not modified
		long double lx = gx;
		long double ly = gy;

		for (;(X-(h/2)) > lx;) 
		{   
			ly += h*(*yprime)(lx, ly); 
			lx += h;
		}   

		
		// printing the approximate points and the absolute error (if applicable)
		if (is_def_Y)
		{
			std::cout << "(" << lx << ", " << ly << "), ";
			std::cout << "absolute error: " << ( ((*Y)(lx) - ly > 0 )? ((*Y)(lx) - ly) : -1*((*Y)(lx) - ly)) << "\n\n";
		}
		else
		{
			std::cout << "(" << lx << ", " << ly << ")\n";
		}
		return ly;
	}
	
	// Improved Eulers method
	long double ieulers (long double h)
	{
		std::cout << "Improved Euler's method applied to a first order differential equation: (h = " << h << ")\n";
		// local variables so the class ones are not modified
		long double lx = gx;
		long double ly = gy;
		// y_{n+1}^*
		long double ys = gy;

		for (;(X-(h/2)) > lx;)
		{
			ys = ly + h*(*yprime)(lx, ly);
			ly += (((*yprime)(lx, ly) + (*yprime)( (lx+h), ys))/2)*h;
			lx += h;
		}

		// printing the approximate points and the absolute error (if applicable)
		if (is_def_Y)
		{
			std::cout << "(" << lx << ", " << ly << "), ";
			std::cout << "absolute error: " << ( ((*Y)(lx) - ly > 0 )? ((*Y)(lx) - ly) : -1*((*Y)(lx) - ly)) << "\n\n";
		}
		else
		{
			std::cout << "(" << lx << ", " << ly << ")\n\n";
		}
		return ly;
	}
	// Fourth order Runge-Kutta method
	long double rk4 (long double h)
	{
		std::cout << "Fourth Order Runge-Kutta method applied to a first order differential equation: (h = " << h << ")\n";
		// local variables so the class ones are not modified
		long double lx = gx;
		long double ly = gy;
		for (;(X-(h/2)) > lx;)
		{
			long double k [4];
			k[0] = (*yprime)(lx, ly);
			k[1] = (*yprime)(lx + h/2, ly + (h/2)*k[0]);
			k[2] = (*yprime)(lx + h/2, ly + (h/2)*k[1]);
			k[3] = (*yprime)(lx + h, ly + h*k[2]);

			ly += (h/6)*(k[0] + 2*k[1] + 2*k[2] + k[3]);
			lx += h;
		}

		// printing the approximate points and the absolute error (if applicable)
		if (is_def_Y)
		{
			std::cout << "(" << lx << ", " << ly << "), ";
			std::cout << "absolute error: " << ( ((*Y)(lx) - ly > 0 )? ((*Y)(lx) - ly) : -1*((*Y)(lx) - ly)) << "\n\n";
		}
		else
		{
			std::cout << "(" << lx << ", " << ly << ")\n\n";
		}
		return ly;
	}
protected:
	// rk4, one step only
	long double RK4o (long double h, long double& Lx, long double& Ly)
	{
		long double k [4];
		k[0] = (*yprime)(Lx, Ly);
		k[1] = (*yprime)(Lx + h/2, Ly + (h/2)*k[0]);
		k[2] = (*yprime)(Lx + h/2, Ly + (h/2)*k[1]);
		k[3] = (*yprime)(Lx + h, Ly + h*k[2]);

		Ly += (h/6)*(k[0] + 2*k[1] + 2*k[2] + k[3]);
		Lx += h;
		return Ly;
	}
public:
	// Adams–Bashforth–Moulton method
	long double abm (long double h)
	{
		std::cout << "Adams-Bashforth-Moulton Method applied to a first order differential equation: (h = " << h << ")\n";
		// local variables so the class ones are not modified
		long double lx = gx;
		long double ly = gy;

		long double YY [4];
		long double * y = YY;
		// y_n values: [n-3, n-2, n-1, n]

		long double YP [4];
		long double * yp = YP;
		// y_n' values: [n-3', n-2', n-1', n']

		// Using the 4th order Runge-Kutta methods to find y_1, y_2, and y_3
		*(y) = ly;
		*(yp) = (*yprime)(lx, ly);

		for (int i=1; i<=3; i++)
		{
			*(y+i) = this->RK4o(h, lx, ly);
			*(yp+i) = (*yprime)(lx, ly);
		}

		for (int i=4; lx < (X-(h/2)) ;i++)
		{
			// ys is y_{n+1}^*
			long double ys = *(y+3) + (h/24)*( 55*(*(yp+3)) - 59*(*(yp+2)) + 37*(*(yp+1)) - 9*(*(yp)) );

			ly += (h/24)*( 9*((*yprime)(lx+h, ys)) + 19*(*(yp+3)) - 5*(*(yp+2)) + (*(yp+1)) );

			// shifting y values to clear space for y_{n+1}
			*y = *(y+1);
			*(y+1) = *(y+2);
			*(y+2) = *(y+3);

			// shifting yp values to clear space for yp_{n+1}
			*yp = *(yp+1);
			*(yp+1) = *(yp+2);
			*(yp+2) = *(yp+3);

			*(y+3) = ly;
			lx += h;
			*(yp+3) = (*yprime)(lx, *(y+3));
		}

		if (is_def_Y)
		{
			std::cout << "(" << lx << ", " << *(y+3) << "), ";
			std::cout << "absolute error: " << ( ((*Y)(lx) - ly > 0 )? ((*Y)(lx) - ly) : -1*((*Y)(lx) - ly)) << "\n\n";
		}
		else
		{
			std::cout << "(" << lx << ", " << *(y+3) << ")\n\n";
		}

		return ly;
	}
};


class sndOrder
{
protected:
	// function pointer ydprime which is y'' = u' = f(x, y, u)
	long double (*ydprime)(long double, long double, long double);

	// function pointer Y wich is the exact solution 
	long double (*Y)(long double);

	// x is x_n, and y is y_n, u is u_n
	long double x;
	long double y;
	long double u;

	// X is the x-coordinate of the point we wish to approximate
	long double X;

	// value to determine whether an exact solution was given
	bool is_def_Y;
public:
	sndOrder
	(
		long double (*f)(long double, long double, long double), 
		long double cx, 
		long double cy, 
		long double cu, 
		long double cX
	) : ydprime(f), x(cx), y(cy), u(cu), X(cX)
	{
		is_def_Y = 0;
	}
	sndOrder
	(
		long double (*f)(long double, long double, long double), 
		long double (*cY)(long double), 
		long double cx, 
		long double cy, 
		long double cu, 
		long double cX
	) : ydprime(f), Y(cY), x(cx), y(cy), u(cu), X(cX)
	{
		is_def_Y = 1;
	}

	// Eulers method
	long double eulers (long double h)
	{
		std::cout << "Euler's method applied to a second order differential equation: (h = " << h << ")\n";
		// local variables so the class ones are not modified
		long double lx = x;
		long double ly = y;
		long double lu = u;

		for (;(X-(h/2)) > lx;) 
		{   
			ly += h*lu;
			lu += h*(*ydprime)(lx, ly, lu);
			lx += h;
		}   
		if (is_def_Y)
		{
			std::cout << "(" << lx << ", " << ly << "), ";
			std::cout << "absolute error: " << ( ((*Y)(lx) - ly > 0 )? ((*Y)(lx) - ly) : -1*((*Y)(lx) - ly)) << "\n\n";
		}
		else
		{
			std::cout << "(" << lx << ", " << ly << ")\n\n";
		}
		return ly;
	}

	// 4th order Runge-Kutta method
	long double rk4 (long double h)
	{
		std::cout << "4th order Runge-Kutta method applied to a second order differential equation: (h = " << h << ")\n";
		// local variables so the class ones are not modified
		long double lx = x;
		long double ly = y;
		long double lu = u;

		long double k [4];
		long double m [4];

		for (;(X-(h/2)) > lx;) 
		{   
			k[0] = (*ydprime)(lx, ly, lu);
			m[0] = lu;
			k[1] = (*ydprime)((lx + h/2), (ly + h*m[0]/2), (lu + h*k[0]/2));
			m[1] = (lu + h*k[0]/2);
			k[2] = (*ydprime)(lx + h/2, ly + h*m[1]/2, lu + h*k[1]/2);
			m[2] = lu + h*k[1]/2;
			k[3] = (*ydprime)(lx + h, ly + h*m[2], lu + h*k[2]);
			m[3] = lu + h*k[2];

			ly += h*(m[0] + 2*m[1] + 2*m[2] + m[3])/6;
			lu += h*(k[0] + 2*k[1] + 2*k[2] + k[3])/6;
			lx += h;
		}

		if (is_def_Y)
		{
			std::cout << "(" << lx << ", " << ly << "), ";
			std::cout << "absolute error: " << ( ((*Y)(lx) - ly > 0 )? ((*Y)(lx) - ly) : -1*((*Y)(lx) - ly)) << "\n\n";
		}
		else
		{
			std::cout << "(" << lx << ", " << ly << ")\n\n";
		}

		return ly;
	}
};


class fstOrderSys
{
protected:
	// function pointers yprime and xprime which are y' = g(t, x, y) and x' = f(t, x, y)
	long double (*yprime)(long double, long double, long double);
	long double (*xprime)(long double, long double, long double);

	// function pointers, exact solutions y(t) and x(t)
	long double (*Y)(long double);
	long double (*X)(long double);

	// x is x_n, y is y_n, and t is t_n
	long double x;
	long double y;
	long double t;

	long double T;

	// value to determine whether an exact solution was given
	bool is_def_Y;
public:
	// constructor: fstOrderSys (xprime, yprime, [X], [Y], t_0, x_0, y_0, T) 
	fstOrderSys
	(
		long double (*XX)(long double, long double, long double), 
		long double (*YYY)(long double, long double, long double), 
		long double ct, 
		long double cx, 
		long double cy, 
		long double cT
	) : xprime(XX), yprime(YYY), t(ct), x(cx), y(cy), T(cT)
	{
		is_def_Y = 0;
	}
	fstOrderSys
	(
		long double (*XX)(long double, long double, long double),
		long double (*YYY)(long double, long double, long double), 
		long double (*Xt)(long double), 
		long double (*Yt)(long double), 
		long double ct, 
		long double cx, 
		long double cy, 
		long double cT
	) : xprime(XX), yprime(YYY), X(Xt), Y(Yt), t(ct), x(cx), y(cy), T(cT)
	{
		is_def_Y = 1;
	}
	// 4th order Runge-Kutta method
	long double * rk4 (long double h)
	{
		std::cout << "4th order Runge-Kutta method applied to a system of first order differential equations: (h = " << h << ")\n";
		// local variables so the class ones are not modified
		long double lx = x;
		long double ly = y;
		long double lt = t;
		long double k [4];
		long double m [4];

		// return value: a pointer to an array: [x_n, y_n, [error in x_n], [error in y_n]]
		long double * final_values;
		if (is_def_Y)
		{
			final_values = new long double [4];
		}
		else
		{
			final_values = new long double [2];
		}

		for (;(T-(h/2)) > lt;) 
		{   
			k[0] = (*yprime)(lt, lx, ly);
			m[0] = (*xprime)(lt, lx, ly);
			k[1] = (*yprime)(lt + h/2, lx + h*m[0]/2, ly + h*k[0]/2);
			m[1] = (*xprime)(lt + h/2, lx + h*m[0]/2, ly + h*k[0]/2);
			k[2] = (*yprime)(lt + h/2, lx + h*m[1]/2, ly + h*k[1]/2);
			m[2] = (*xprime)(lt + h/2, lx + h*m[1]/2, ly + h*k[1]/2);
			k[3] = (*yprime)(lt + h, lx + h*m[2], ly + h*k[2]);
			m[3] = (*xprime)(lt + h, lx + h*m[2], ly + h*k[2]);

			ly += (h/6)*(k[0] + 2*k[1] + 2*k[2] + k[3]);
			lx += (h/6)*(m[0] + 2*m[1] + 2*m[2] + m[3]);
			lt += h;
		}   

		// printing the approximate points and the absolute error (if applicable)
		if (is_def_Y)
		{
			std::cout << "(t, x, y): " << "(" << lt << ", " << lx << ", " << ly << ")\n";
			std::cout << "error in x: " << ( ((*X)(lt) - lx > 0 )? ((*X)(lt) - lx ) : -1*((*X)(lt) - lx)) << "\n";
			std::cout << "error in y: " << ( ((*Y)(lt) - ly > 0 )? ((*Y)(lt) - ly ) : -1*((*Y)(lt) - ly)) << "\n\n";
			// return values
			final_values[0] = lx;
			final_values[1] = ly;
			final_values[2] = ( ((*X)(lt) - lx > 0 )? ((*X)(lt) - lx ) : -1*((*X)(lt) - lx));
			final_values[3] = ( ((*Y)(lt) - ly > 0 )? ((*Y)(lt) - ly ) : -1*((*Y)(lt) - ly));
		}
		else
		{
			std::cout << "(t, x, y): " << "(" << lt << ", " << lx << ", " << ly << ")\n\n";
			// return values
			final_values[0] = lx;
			final_values[1] = ly;
		}

		return final_values;
	}
};
