/*============================================================


Name: finite_difference.h

Author: Louis Horner

Date: January 27, 2021

Description:
	Here we define the class finDiff to approximate the 
mesh points of a second-order boundary-value problem using 
the finite difference method, with the central difference 
approximations. It approximates solutions of problems of the
following form.

d^2y/dx^2 + P(x) dy/dx + Q(x) y = f(x)
y(a) = alpha
y(b) = beta

finDiff
	Dependency: matrix.h

	Constructor:
		finDiff (P, Q, f, a, b, alpha, beta)

	Public members:
		solve(n)
			n is the number of mesh points
			prints the values of y_i on the last column of 
			a matrix


============================================================*/


class finDiff
{
protected:
	long double (*P)(long double);
	long double (*Q)(long double);
	long double (*f)(long double);
	long double a; long double b; long double alpha; long double beta;
public:
	finDiff
	(
		long double (*tP)(long double),
		long double (*tQ)(long double),
		long double (*tf)(long double),
		long double ta,
		long double talpha,
		long double tb,
		long double tbeta
	) : P(tP), Q(tQ), f(tf), a(ta), alpha(talpha), b(tb), beta(tbeta)
	{ }

	void solve (int n)
	{
		long double h = ((long double) (b - a))/n;
		matrix M (n-1, n);
		long double x = a;
		int i = 1;
	// filling the matrix
		// first entry
		x += h;
		M.M[0][0] = -2.0 + h*h*(*Q)(x);
		M.M[0][1] = 1.0 + (h/2.0)*(*P)(x);
		M.M[0][n-1] = h*h*(*f)(x) - (1.0 - (h/2.0)*(*P)(x))*alpha;
		// middle entries
		x += h;
		for (;i<n-2; i++)
		{
			M.M[i][i-1] = 1.0 - (h/2.0)*(*P)(x);
			M.M[i][i] = -2.0 + h*h*(*Q)(x);
			M.M[i][i+1] = 1.0 + (h/2.0)*(*P)(x);
			M.M[i][n-1] = h*h*(*f)(x);
			x += h;
		}
		// last entry
		M.M[i][i-1] = 1.0 - (h/2.0)*(*P)(x);
		M.M[i][i] = -2.0 + h*h*(*Q)(x);
		M.M[i][n-1] = h*h*(*f)(x) - (1.0 + (h/2.0)*(*P)(x))*beta;
		// solving and printing the solutions
		M.rref();
		M.print();
	}
};
