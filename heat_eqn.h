/*============================================================


Name: heat_eqn.h

Author: Louis Horner

Date: February 9, 2021

Description:
	This file defines the class heat_eqn meant to solve
the heat equation, or problems of the following form:

u(x, t)
c d^2u/dx^2 = du/dt
u(x,0) = uxi
u(0,0) = T0
u(a,0) = Ta

heat_eqn
	Dependency: matrix.h

	Constructor:
		heat_eqn (u(x,0), c, a, T, Ta, T0)

	Public members:
		exp (int n, int m)
			explicit finite difference

		cn (int n, int m)
			Crank-Nicholson Method, an implicit finite difference 
			method


============================================================*/
class heat_eqn
{
protected:
	double (*uxi)(double);
	double c;
	double a;
	double T;
	double Ta;
	double T0;
public:
	heat_eqn (
		double (*UXI)(double),
		double C,
		double A,
		double Tt,
		double Tr,
		double Tl
	) : uxi(UXI), c(C), a(A), T(Tt), Ta(Tr), T0(Tl)
	{}

	void exp (int n, int m)
	{
		double h = (double) a/n;
		double k = (double) T/m;
		double lambda = (double) c*k/(h*h);

		matrix M (m, n);
		// time increases by k for each increase in row, x increases by h for each increase in column
		// the leftmost column is time

		for (int row=0; row<m; row++)
		{
			M.M[row][0] = k*((double)row + 1.0);
			for (int col=0; col<n-1; col++)
			{
				// first row of M
				if (row == 0)
				{
					// first entry
					if (col == 0)
					{
						M.M[row][col+1] = lambda*T0 + (1.0 - 2.0*lambda)*(*uxi)(h) + lambda*(*uxi)(2.0*h);
					}
					// last entry
					else if (col == n-2)
					{
						M.M[row][col+1] = lambda*(*uxi)(((double)col)*h) + (1.0 - 2.0*lambda)*(*uxi)(((double)col+1.0)*h) + lambda*Ta;
					}
					// middle entries
					else 
					{
						M.M[row][col+1] = lambda*(*uxi)(((double)col)*h) + (1.0 - 2.0*lambda)*(*uxi)(((double)col+1.0)*h) + lambda*(*uxi)(((double)col+2.0)*h);
					}
				}
				// every row after first in M
				else
				{
					// first entry
					if (col == 0)
					{
						M.M[row][col+1] = lambda*T0 + (1.0 - 2.0*lambda)*M.M[row-1][col+1] + lambda*M.M[row-1][col+2];
					}
					// last entry
					else if (col == n-2)
					{
						M.M[row][col+1] = lambda*M.M[row-1][col] + (1.0 - 2.0*lambda)*M.M[row-1][col+1] + lambda*Ta;
					}
					// middle entries
					else 
					{
						M.M[row][col+1] = lambda*M.M[row-1][col] + (1.0 - 2.0*lambda)*M.M[row-1][col+1] + lambda*M.M[row-1][col+2];
					}
				}
			}
		}

		// output
		M.print();
	}
	void cn (int n, int m)
	// Crank-Nicholson Method
	{
		double h = (double) a/n;
		double k = (double) T/m;
		double lambda = (double) c*k/(h*h);
		double alpha = 2.0*(1.0 + (1.0/lambda));
		double beta = 2.0*(1.0 - (1.0/lambda));

		matrix M (m, n);
		// time increases by k for each increase in row, x increases by h for each increase in column
		// the leftmost column is time

		for (int row=0; row<m; row++)
		{
			M.M[row][0] = k*((double)row + 1.0);
			matrix F (n-1, n);

			// forcing F.M values to 0 (because a bug happened if not)
			for (int R=0; R<n-1; R++)
			{
				for (int col=0; col<n; col++)
				{
					F.M[R][col] = 0.0;
				}
			}

			if (row == 0)
			// first row of M
			{
				for (int col=0; col<n-1; col++)
				{
					if (col == 0)
					// first row in F
					{
						F.M[col][col] = alpha;
						F.M[col][col+1] = -1.0;
						F.M[col][n-1] = 2.0*T0 + (*uxi)(((double)col+2.0)*h) - beta*(*uxi)(((double)col+1.0)*h);
					}
					else if (col == n-2)
					// last row in F
					{
						F.M[col][col-1] = -1.0;
						F.M[col][col] = alpha;
						F.M[col][n-1] = 2.0*Ta - beta*(*uxi)(((double)col+1.0)*h) + (*uxi)(((double)col)*h);
					}
					else
					// middle rows
					{
						F.M[col][col-1] = -1.0;
						F.M[col][col] = alpha;
						F.M[col][col+1] = -1.0;
						F.M[col][n-1] = (*uxi)(((double)col+2.0)*h) - beta*(*uxi)(((double)col+1.0)*h) + (*uxi)(((double)col)*h);
					}
				}
			}
			else
			// every row of M not the first
			{
				for (int col=0; col<n-1; col++)
				{
					if (col == 0)
					// first row in F
					{
						F.M[col][col] = alpha;
						F.M[col][col+1] = -1.0;
						F.M[col][n-1] = 2.0*T0 + M.M[row-1][col+2] - beta*M.M[row-1][col+1];
					}
					else if (col == n-2)
					// last row in F
					{
						F.M[col][col-1] = -1.0;
						F.M[col][col] = alpha;
						F.M[col][n-1] = 2.0*Ta - beta*M.M[row-1][col+1] + M.M[row-1][col];
					}
					else
					// middle rows
					{
						F.M[col][col-1] = -1.0;
						F.M[col][col] = alpha;
						F.M[col][col+1] = -1.0;
						F.M[col][n-1] = M.M[row-1][col+2] - beta*M.M[row-1][col+1] + M.M[row-1][col];
					}
				}
			}

			// solving for the next row
			F.rref();

			// transferring the next row from F to M
			for (int i=0; i<n-1; i++)
			{
				M.M[row][i+1] = F.M[i][n-1];
			}
		}

		// output
		M.print();
	}
};
