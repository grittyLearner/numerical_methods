/*============================================================


Name: wave_eqn.h

Author: Louis Horner

Date: February 9, 2021

Description:
	This file defines the class wave_eqn meant to solve
the wave equation, or problems of the following form:

u(x, t)
c^2 d^2u/dx^2 = d^2u/dt^2
0 < x < X, 0 < t < T
u(x,0) = uxi
u(0,t) = U0
u(X,t) = UX

wave_eqn
	Dependency: matrix.h

	Constructor:
		heat_eqn (u(x,0), c, X, T, U0, UX)

	Public members:
		solve(int n, int m)
			prints the output


============================================================*/
class wave_eqn
{
protected:
	double (*uxi)(double);
	double c;
	double X;
	double T;
	double U0;
	double UX;
public:
	wave_eqn (
		double (*Uxi)(double),
		double C,
		double x,
		double t,
		double u0,
		double uX
	) : uxi(Uxi), c(C), X(x), T(t), U0(u0), UX(uX)
	{}

	void solve (int n, int m)
	{
		// implementing finite difference method
		double h = (double) X/n;
		double k = (double) T/m;
		double lambda = (double) c*k/h;

		matrix M (m, n);
		// time increases by k for each increase in row, x increases by h for each increase in column
		// the leftmost column is time

		for (int row=0; row<m; row++)
		{
			// time
			M.M[row][0] = k*((double)row + 1.0);

			for (int col=0; col<n-1; col++)
			{
				// first row of M
				if (row == 0)
				{
					// first entry
					if (col == 0)
					{
						M.M[row][col+1] = (lambda*lambda/2.0)*((*uxi)(((double)col+2.0)*h) + U0) + (1.0 - lambda*lambda)*(*uxi)(((double)col+1.0)*h);
					}
					// last entry
					else if (col == n-2)
					{
						M.M[row][col+1] = (lambda*lambda/2.0)*(UX + (*uxi)(((double)col)*h)) + (1.0 - lambda*lambda)*(*uxi)(((double)col+1.0)*h);
					}
					// middle entries
					else if (col < n-2)
					{
						M.M[row][col+1] = (lambda*lambda/2.0)*((*uxi)(((double)col+2.0)*h) + (*uxi)(((double)col)*h)) + (1.0 - lambda*lambda)*(*uxi)(((double)col+1.0)*h);
					}
				}
				// second row of M
				else if (row == 1)
				{
					// first entry
					if (col == 0)
					{
						M.M[row][col+1] = lambda*lambda*M.M[row-1][col+2] + 2.0*(1.0 - lambda*lambda)*M.M[row-1][col+1] + lambda*lambda*U0 - (*uxi)(((double)col+1.0)*h);
					}
					// last entry
					else if (col == n-2)
					{
						M.M[row][col+1] = lambda*lambda*UX + 2.0*(1.0 - lambda*lambda)*M.M[row-1][col+1] + lambda*lambda*M.M[row-1][col] - (*uxi)(((double)col+1.0)*h);
					}
					// middle entries
					else 
					{
						M.M[row][col+1] = lambda*lambda*M.M[row-1][col+2] + 2.0*(1.0 - lambda*lambda)*M.M[row-1][col+1] + lambda*lambda*M.M[row-1][col] - (*uxi)(((double)col+1.0)*h);
					}
				}
				// every row after second in M
				else
				{
					// first entry
					if (col == 0)
					{
						M.M[row][col+1] = lambda*lambda*M.M[row-1][col+2] + 2.0*(1.0 - lambda*lambda)*M.M[row-1][col+1] + lambda*lambda*U0 - M.M[row-2][col+1];
					}
					// last entry
					else if (col == n-2)
					{
						M.M[row][col+1] = lambda*lambda*UX + 2.0*(1.0 - lambda*lambda)*M.M[row-1][col+1] + lambda*lambda*M.M[row-1][col] - M.M[row-2][col+1];
					}
					// middle entries
					else 
					{
						M.M[row][col+1] = lambda*lambda*M.M[row-1][col+2] + 2.0*(1.0 - lambda*lambda)*M.M[row-1][col+1] + lambda*lambda*M.M[row-1][col] - M.M[row-2][col+1];
					}
				}
			}
		}

		// output
		M.print();
	}
};
