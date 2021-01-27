/*============================================================


Name: finite_difference.h

Author: Louis Horner

Date: January 27, 2021

Description:
	Here we define the necessary classes to approximate the 
mesh points of a second-order boundary-value problem using 
the finite difference method, with the central difference 
approximations.

Classes:

matrix
	Description: This class provides the necessary framework
that finDiff uses.
	Constructor:
		matrix (n_rows, n_cols)
			defines a n_rows by n_cols matrix
	Public members:
		rref()
			reduces the matrix to reduced row-echelon form
		print()

finDiff
	Description: This class approximates the mesh points of 
a second-order boundary-value problem using the finite
difference method. It uses the following for notation purposes.
	d^2y/dx^2 + P(x) dy/dx + Q(x) y = f(x)
	y(a) = alpha
	y(b) = beta

	Constructor:
		finDiff (P, Q, f, a, b, alpha, beta)

	Public members:
		solve(n)
			n is the number of mesh points
			prints the values of y_i on the last column of 
			  a matrix


============================================================*/
class matrix
{
protected:
	int n_rows, n_cols;
	long double ** M;

public:
	friend class finDiff;
	// constructor: matrix (n_rows, n_cols)
	matrix (int tn_rows, int tn_cols) : n_rows(tn_rows), n_cols(tn_cols)
	{
		M = new long double * [tn_rows];
		for (int r=0; r<tn_rows; r++)
		{
			M[r] = new long double [tn_cols];
		}
	}
	// easier initialization for testing rref
	matrix& operator= (long double rhs [3][4])
	{
		for (int row=0; row<3;row++)
		{
			for (int col=0; col<4;col++)
			{
				M[row][col] = rhs[row][col];
			}
		}
		return *this;
	}
	// destructor
	~matrix ()
	{
		for (int r=0; r<n_rows; r++)
		{
			delete M[r];
		}
		delete M;
	}
	// printing
	void print ()
	{
		for (int row=0; row<n_rows; row++)
		{
			for (int col=0; col<n_cols; col++)
			{
				std::cout << M[row][col] << ", ";
			}
			std::cout << "\n";
		}
	}
	// using Gaussian elimination to turn *this into reduced row echelon form
	matrix& rref ()
	{
		for (int R=0; R<n_rows; R++)
		// R is the index of the working diagonal
		{
			// making the diagonal at R,R one by dividing its row by the value
			long double c = M[R][R];
			for (int col=R; col<n_cols; col++)
			{
				M[R][col] /= c;
			}
			// making col below the 1 equal 0
			for (int row=R+1; row<n_rows; row++)
			{
				long double cnst = M[row][R];
				for (int col=0; col<n_cols; col++)
				{
					M[row][col] -= M[R][col] * cnst;
				}
			}
			// making col above the 1 equal 0
			for (int row=0; row<R; row++)
			{
				long double cnst = M[row][R];
				for (int col=R; col<n_cols; col++)
				{
					M[row][col] -= M[R][col] * cnst;
				}
			}
		}
		return *this;
	}
};


class finDiff
{
protected:
	long double (*P)(long double);
	long double (*Q)(long double);
	long double (*f)(long double);
	long double a;
	long double b;
	long double alpha;
	long double beta;
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
