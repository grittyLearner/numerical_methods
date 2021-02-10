/*============================================================


Name: laplace_eqn.h

Author: Louis Horner

Date: February 9, 2021

Description:
	Here we define the class laplace to approximate the 
interior mesh points of a second-order boundary-value problem
using the finite difference method. It solves problems of the
following form.

div^2 u = 0, where u(x, y)
0 < x < n
0 < y < n
u(x,0) = x0
u(0,y) = y0
u(x,n) = ux
u(n,y) = uy

laplace
	Constructor:
		laplace (u(x,n), u(n,y), x0, y0, n)

	Public members:
		solve(double h)
			prints a matrix with the mesh points on the right
			column.


============================================================*/
class laplace
{
protected:
	double n;
	double x0;
	double y0;
	double (*uy)(double);
	double (*ux)(double);
public:
	laplace (
		double (*uY)(double),
		double (*uX)(double),
		double X0,
		double Y0,
		double N
	) : uy(uY), ux(uX), x0(X0), y0(Y0), n(N)
	{}
	void solve (double h)
	{
		int r = (int) n/h - 1;
		int n_mesh_pts = r*r;

		// filling the matrix with respect to the square
		matrix M (n_mesh_pts, n_mesh_pts+1);
		// bottom left corner
		M.M[0][0] = -4.0;  // center
		M.M[0][1] = 1.0;  // right
		M.M[0][r] = 1.0;  // up
		M.M[0][n_mesh_pts] = -1*x0 - y0;  // rhs
		// bottom right corner
		M.M[r-1][r-1] = -4.0;  // center
		M.M[r-1][2*r-1] = 1.0;  // up
		M.M[r-1][r-2] = 1.0;  // left
		M.M[r-1][n_mesh_pts] = -1*x0 - (*uy)(h);  // rhs
		// top left corner
		M.M[n_mesh_pts-r][n_mesh_pts-r] = -4.0;  // center
		M.M[n_mesh_pts-r][n_mesh_pts-r+1] = 1.0;  // right
		M.M[n_mesh_pts-r][n_mesh_pts-2*r] = 1.0;  // down
		M.M[n_mesh_pts-r][n_mesh_pts] = -1*y0 - (*ux)(h);  // rhs
		// top right corner
		M.M[n_mesh_pts-1][n_mesh_pts-1] = -4.0;  // center
		M.M[n_mesh_pts-1][n_mesh_pts-2] = 1.0;  // left
		M.M[n_mesh_pts-1][n_mesh_pts-r-1] = 1.0;  // down
		M.M[n_mesh_pts-1][n_mesh_pts] = -1*(*ux)(n-h) - (*uy)(n-h);  // rhs

		// botton row:
		for (int P=1;P<r-1; P++)
		{
			M.M[P][P] = -4.0;  // center
			M.M[P][P+1] = 1.0;  // right
			M.M[P][P-1] = 1.0;  // left
			M.M[P][P+r] = 1.0;  // up
			M.M[P][n_mesh_pts] = -1*x0;  // rhs
		}
		// left column:
		for (int P=r;P<n_mesh_pts-r; P+=r)
		{
			M.M[P][P] = -4.0;  // center
			M.M[P][P+1] = 1.0;  // right
			M.M[P][P-r] = 1.0;  // down
			M.M[P][P+r] = 1.0;  // up
			M.M[P][n_mesh_pts] = -1*y0;  // rhs
		}
		// top row:
		for (int P=n_mesh_pts-r+1;P<n_mesh_pts-1; P++)
		{
			M.M[P][P] = -4.0;  // center
			M.M[P][P+1] = 1.0;  // right
			M.M[P][P-1] = 1.0;  // left
			M.M[P][P-r] = 1.0;  // down
			M.M[P][n_mesh_pts] = -1*(*ux)(P%r);  // rhs
		}
		// right column:
		for (int P=2*r-1;P<n_mesh_pts-r; P+=r)
		{
			M.M[P][P] = -4.0;  // center
			M.M[P][P-1] = 1.0;  // left
			M.M[P][P-r] = 1.0;  // down
			M.M[P][P+r] = 1.0;  // up
			M.M[P][n_mesh_pts] = -1*(*uy)(((P/r)+1)*h);  // rhs
		}
		// centered points
		for (int P=r+1;P<n_mesh_pts-r; P+=r)
		{
			M.M[P][P] = -4.0;  // center
			M.M[P][P-1] = 1.0;  // left
			M.M[P][P+1] = 1.0;  // right
			M.M[P][P-r] = 1.0;  // down
			M.M[P][P+r] = 1.0;  // up
		}

		// solving the matrix
		M.rref();
		// rref = solver
		M.print();
	}
};
