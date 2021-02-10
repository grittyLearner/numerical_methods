/*============================================================


Name: matrix.h

Author: Louis Horner

Date: February 9, 2021

Description:
	Here we define the class matrix to handle 2 dimentional
real-values matricies.


matrix
	Constructor:
		matrix (n_rows, n_cols)
			defines a n_rows by n_cols matrix
	Public members:
		rref()
			reduces the matrix to reduced row-echelon form
			using Gaussian elimination
		print()


============================================================*/


class matrix
{
protected:
	int n_rows, n_cols;
	long double ** M;

public:
	friend class finDiff;
	friend class laplace_eqn;
	friend class heat_eqn;
	friend class wave_eqn;
	matrix (int tn_rows, int tn_cols) : n_rows(tn_rows), n_cols(tn_cols)
	{
		M = new long double * [tn_rows];
		for (int r=0; r<tn_rows; r++)
		{
			M[r] = new long double [tn_cols];
		}
	}
	// initialization for testing
	matrix& operator= (long double rhs [2][3])
	{
		for (int row=0; row<2;row++)
		{
			for (int col=0; col<3;col++)
			{
				M[row][col] = rhs[row][col];
			}
		}
		return *this;
	}
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
	// matrix multiplication
	// C = A * B, A = *this, B = rhs
	long double ** operator* (matrix& rhs)
	{
		if (n_cols == rhs.n_rows)
		{
			long double ** C = new long double * [n_rows];
			for (int r=0; r<n_rows; r++)
			{
				C[r] = new long double [rhs.n_cols];
			}

			long double sum = 0.0;
			for (int row_a = 0; row_a < n_rows; row_a++)
			{
				for (int col_b = 0; col_b < rhs.n_cols; col_b++)
				{
					for (int col_a = 0; col_a < n_cols; col_a++)
					{
						sum += M[row_a][col_a]*rhs.M[col_a][col_b];
					}
					C[row_a][col_b] = sum;
					sum = 0.0;
				}
			}
			return C;
		}
		else
		{
			std::cout << "Error: you tried to multiply two matricies where the numer of rows in the second did not equal the numer of colums in the first\n";
			throw "Error";
		}
	}
	// copy assignment
	matrix& operator= (long double ** rhs)
	{
		// copying rhs.M to this->M
		for (int row_r=0; row_r<n_rows; row_r++)
		{
			for (int col_r=0; col_r<n_cols; col_r++)
			{
				M[row_r][col_r] = rhs[row_r][col_r];
			}
		}

		// deleting rhs
		for (int r=0; r<n_rows; r++)
		{
			delete rhs[r];
		}
		delete rhs;

		return *this;
	}
};
