package linAlg;

/**
 * @author Steven Lawrence
 */

public class Matrix {
	protected double[][] A; // 2d array to hold the matrix information
	private int m, n; // the matrix has m rows and n columns

	/**
	 * Constructs a mxn matrix zero-d out. Note Matrix indexes start at 1 and not 0.
	 * 
	 * @param m
	 *            the rows of the matrix
	 * @param n
	 *            the columns of the matrix
	 */
	public Matrix(int m, int n) {
		this.m = m;
		this.n = n;
		A = new double[m][n];
	}

	/**
	 * Constructs a mxn matrix from an already established 2d array.
	 * 
	 * @param values
	 *            a 2d array list
	 */
	public Matrix(double[][] values) {
		A = values;
		m = values.length;
		n = values[0].length;
	}

	/**
	 * Construct a copy of a Matrix from one that is already made
	 * 
	 * @param C
	 */
	public Matrix(Matrix C) {
		m = C.m;
		n = C.n;
		A = new double[m][n];
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = C.A[i][j];
			}
		}
	}
	
	public Matrix(double[] vector) {
		m = 1;
		n = vector.length;
		A = new double[m][n];
		for (int i = 0; i < n; i++) {
			A[0][i] = vector[i];
		}
	}
	public static Matrix ones (Matrix C) {
		Matrix oneMat = new Matrix(C);
		for (int i = 0; i < oneMat.getRowSize(); i++) {
			for (int j = 0; j < oneMat.getColSize(); j++) {
				oneMat.A[i][j] = 1.0;
			}
		}
		return oneMat;
	}
	
	/**
	 * 
	 * @param poly
	 */
	public Matrix(Polynomial p) {
		m = p.deg;
		n = p.deg;
		A = new double[m][n];
		Polynomial poly = p;
		
		if (!poly.isMonic()) {
			poly = new Polynomial(p);
			poly.makeMonic();
		}
		for(int i = n; i > 0; i--) {
			double coeff = poly.coefficients[i];
			setEntry(1, i, -1 * coeff);
		}
		for (int k = 1; k < n; k++) {
			setEntry(k+1, k, 1);
		}
	}
	
	public int getRowCount() {
		return m;
	}
	
	public int getColCount() {
		return n;
	}

	/**
	 * Sets the entry in the ith row and jth column
	 */
	public void setEntry(int i, int j, double val) {
		A[i - 1][j - 1] = val;
	}

	/**
	 * get an entry
	 */
	public double getEntry(int i, int j) {
		return A[i - 1][j - 1];
	}

	/**
	 * Adds another matrix to itself
	 * 
	 * @param B
	 *            the matrix being added
	 */
	public Matrix add(Matrix B) {
		// make sure the dimensions are equal
		String s = "addition is undefined for matrices of different dimensions"+
				m + " " + n + "and" + B.m + " " + B.n;
		if (m != B.m || n != B.n) {
			throw new IllegalArgumentException(s);
		}
		else {
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < n; j++) {
					A[i][j] = A[i][j] + B.A[i][j];
				}
			}
		}
		return this;
	}

	/**
	 * Static function for matrix addition. Accepts two matrices and returns their
	 * sum
	 * 
	 * @param mat1
	 *            a matrix
	 * @param mat2
	 *            another matrix
	 * @return mat1 + mat2
	 */
	public static Matrix getSum(Matrix mat1, Matrix mat2) {
		Matrix newMat = new Matrix(mat1);
		newMat.add(mat2);
		return newMat;
	}

	/**
	 * subtract another matrix from itself
	 * 
	 * @param B
	 *            the matrix to subtract with
	 */
	public Matrix minus(Matrix B) {
		// make sure the dimensions are equal
		if (m == B.m && n == B.n) {
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < n; j++) {
					A[i][j] = A[i][j] - B.A[i][j];
				}
			}
		} else {
			System.out.println("Subtraction is undefined for matrices of different dimensions");
			System.out.println(m + " " + n);
			System.out.println(B.m + " " + B.n);
		}
		return this;
	}

	/**
	 * Static function for subtraction. Accepts two matrices and returns their
	 * difference
	 * 
	 * @param mat1
	 *            a matrix
	 * @param mat2
	 *            another matrix
	 * @return mat1 - mat2
	 */
	public static Matrix getDiff(Matrix mat1, Matrix mat2) {
		Matrix newMat = new Matrix(mat1);
		newMat.minus(mat2);
		return newMat;
	}

	/**
	 * Multiplies another matrix to itself
	 * 
	 * @param B
	 *            the matrix being multiplied with
	 */
	public Matrix times(Matrix B) {
		// Check to make sure the dimension of columns of the first matrix matches the
		// dimension of rows of the second
		if (n == B.m) {
			Matrix C = new Matrix(m, B.n); // the product matrix's dimensions are m x B.n
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < B.n; j++) {
					// dot the row vectors of the first matrix with the column vectors of the second
					// matrix
					double dotProd = 0;
					for (int k = 0; k < n; k++) {
						dotProd = dotProd + A[i][k] * B.A[k][j];
					}
					C.A[i][j] = dotProd;
				}
			}
			// set the product to the current matrix
			A = C.A;
			m = C.m;
			n = C.n;
		} else {
			System.out.println("Multiplication is undefined for the given matrices' dimensions");
			System.out.println(this);
			System.out.println(B);
		}
		return this;
	}

	/**
	 * Static function for multiplying matrices. Accepts two matrices and returns
	 * their product
	 * 
	 * @param mat1
	 *            a matrix
	 * @param mat2
	 *            another matrix
	 * @return mat1 x mat2
	 */
	public static Matrix getMult(Matrix mat1, Matrix mat2) {
		Matrix newMat = new Matrix(mat1);
		newMat.times(mat2);
		return newMat;
	}
	
	public Matrix hadamard(Matrix m2) {
		if (n == m2.n && m == m2.m) {
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < n; j++) {
					A[i][j] *= m2.A[i][j];
				}
			}
		}
		return this;
	}
	public static Matrix getHadamard(Matrix m1, Matrix m2) {
		Matrix newMat = new Matrix(m1);
		return newMat.hadamard(m2);
	}

	/**
	 * Scale up/down the entries of the matrix
	 * 
	 * @param s
	 *            scaling factor
	 */
	public Matrix scale(double s) {
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] *= s;
			}
		}
		return this;
	}

	/**
	 * Static method for scaling up a matrix. Scale up/down the entries of the
	 * matrix
	 * 
	 * @param mat
	 *            the matrix being scaled
	 * @param s
	 *            scaling factor
	 */
	public static Matrix scale(Matrix mat, double s) {
		Matrix C = new Matrix(mat);
		for (int i = 0; i < mat.m; i++) {
			for (int j = 0; j < mat.n; j++) {
				C.A[i][j] *= s;
			}
		}
		return C;
	}

	/**
	 * Switches the row vectors and the column vectors
	 */
	public void transpose() {
		Matrix Atp = new Matrix(n, m); // the transpose is nxm
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				Atp.A[j][i] = A[i][j];
			}
		}
		A = Atp.A;
		m = Atp.m;
		n = Atp.n;
	}
	
	public int getRowSize() {
		return m;
	}
	
	public int getColSize() {
		return n;
	}
	
	public double[] rowToDoubleArray(int row) {
		double[] d = new double[getColSize()];
		for (int i = 1; i <= getColSize(); i++) {
			double val = getEntry(row, i);
			d[i - 1] = val;
		}
		return d;
	}

	/**
	 * Static method for transposing a matrix. Switches the row vectors and the
	 * column vectors
	 * 
	 * @return the transpose matrix that is nxm
	 */
	public static Matrix getTranspose(Matrix At) {
		Matrix Atran = new Matrix(At.n, At.m);
		for (int i = 0; i < At.m; i++) {
			for (int j = 0; j < At.n; j++) {
				Atran.A[j][i] = At.A[i][j];
			}
		}
		return Atran;
	}

	/**
	 * switch rows in a matrix
	 * 
	 * @param row1
	 * @param row2
	 */
	public void swap(int row1, int row2) {
		if (row1 < 1 || row1 > m || row2 < 1 || row2 > m) {
			System.out.println("Index out of bounds");
		} else {
			double[] temp = A[row1 - 1].clone();
			A[row1 - 1] = A[row2 - 1];
			A[row2 - 1] = temp;
		}
	}

	/**
	 * Return the identity matrix of the matrix
	 * 
	 * @return identity matrix
	 */
	public Matrix getIdentity() {
		Matrix i = new Matrix(m, n);
		int d;
		d = (m <= n) ? m : n;
		for (int diag = 0; diag < d; diag++) {
			i.A[diag][diag] = 1;
		}
		return i;
	}

	/**
	 * Reduces the matrix into reduced row echelon form
	 * 
	 * @return integer flag on whether the matrix was sucessfully reduced
	 */
	public int reduce() {
		int E = 1; // flag

		// loop through the columns
		for (int j = 1; j <= m; j++) {
			double pivot = 0;
			int i, p; // row index i and pivot index p

			/*
			 * loop through the rows of the current column and find the largest magnitude to
			 * act as the pivot position. Record the value into pivot and the row index into
			 * p. The purpose of finding the largest magnitude is to reduce the error
			 * introduced by division.
			 */
			p = j;
			for (i = j; i <= n; i++) {
				if (Math.abs(A[i - 1][j - 1]) > pivot) {
					pivot = A[i - 1][j - 1];
					p = i;
				}
			}
			// interchange rows to match proper reduced row echelon form
			if (p > j) {
				swap(p, j);
			}

			// if a column vector is the 0 vector for nxn matrix then the system does not
			// have a unique solution
			if (pivot == 0) {
				E = 0;
				return E;
			}
			// divide the row by the pivot value so that the pivot position is equivalent to
			// 1
			for (int k = 0; k < A[0].length; k++) {
				A[j - 1][k] = A[j - 1][k] / pivot;
			}
			// perform Gauss-Jordan elimination on the other rows to zero out that column
			for (i = 1; i <= m; i++) { // i refers to the column and j refers to the row here
				if (i != j) {
					// since the pivot is one, all you have to do is subtract it scaled up from the
					// other row's corresponding column entry from the other none pivot rows
					double[] rowj = A[j - 1].clone();
					double scaler = A[i - 1][j - 1];
					for (int entry = 0; entry < A[0].length; entry++) {
						rowj[entry] *= scaler;
						A[i - 1][entry] -= rowj[entry];
					}
				}
			}
		}
		return E;
	}

	/**
	 * Augment a given matrix onto itself
	 * 
	 * @param b
	 *            the augmenting matrix
	 */
	private void augment(Matrix b) {
		if (m != b.m) {
			System.out.println("Can't augment matrices of different row dimemsions");
		} else {
			double[][] temp = A;
			A = new double[m][n + b.n]; // new matrix is mx(columns of first plus columns of second)
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < n + b.n; j++) {
					if (j < n) {
						A[i][j] = temp[i][j];
					} else {
						A[i][j] = b.A[i][j - n];
					}
				}
			}
		}
	}

	/**
	 * De-augment the matrix
	 * 
	 * @return the matrix that was de-augmented from the augmented matrix
	 */
	public Matrix deaug() {
		if (A[0].length == n) {
			return this; // return original if matrix is not augmented
		} else {
			Matrix b = new Matrix(m, A[0].length - n);
			double[][] temp = new double[m][n];
			for (int i = 0; i < m; i++) {
				for (int j = 0; j < A[0].length; j++) {
					if (j < n) {
						temp[i][j] = A[i][j];
					} else {
						b.A[i][j - n] = A[i][j];
					}
				}
			}
			A = temp;
			return b;
		}
	}

	/**
	 * Solves the matrix system and gives the matrix. If matrix does not have a
	 * unique solution then it returns the failed to reduced matrix.
	 * 
	 * @param b
	 *            solution vector
	 * @return flag value on whether or not the matrix has a unique solution or not
	 */
	public int solve(Matrix b) {
		augment(b);
		int E = reduce();
		System.out.println(this);
		Matrix soln = deaug();
		if (E == 1) {
			A = soln.A;
			m = soln.m;
			n = soln.n;
		}
		return E;
	}

	/**
	 * Solves matrix system and allows the user to reduce a copy of the matrix so
	 * that the original can be preserved.
	 * 
	 * @param b
	 *            the matrix trying to solve
	 * @param reduceCopy
	 *            flag on whether to create a copy of the solution or to just
	 *            convert the original into the solution
	 * @return the answer to the matrix equation
	 */
	public Matrix solve(Matrix b, boolean reduceCopy) {
		if (reduceCopy) {
			Matrix copyM = new Matrix(this);
			copyM.augment(b);
			copyM.reduce();
			System.out.println(copyM);
			b = copyM.deaug();
			return b;
		} else {
			this.solve(b);
			return this;
		}
	}

	/**
	 * Gaussian elimination solve a matrix by finding the upper triangular matrix
	 * and then using substitution to solve the system of equation
	 * 
	 * @param b
	 *            the matrix trying to be solved
	 * @return the solution matrix to the system of equation
	 */
	public int gaussianElim(Matrix b) {
		int E = 1; // flag
		augment(b);

		// loop through the columns
		for (int j = 1; j <= m; j++) {
			double pivot = 0;
			int i, p; // row index i and pivot index p

			/*
			 * loop through the rows of the current column and find the largest magnitude to
			 * act as the pivot position. Record the value into pivot and the row index into
			 * p. The purpose of finding the largest magnitude is to reduce the error
			 * introduced by division.
			 */
			p = j;
			for (i = j; i <= n; i++) {
				if (Math.abs(A[i - 1][j - 1]) > pivot) {
					pivot = A[i - 1][j - 1];
				}
				// interchange rows in order to get upper triangular matrix
				if (p > j) {
					swap(p, j);
				}
			}

			// if a column vector is the 0 vector for nxn matrix then the system does not
			// have unique solution
			if (pivot == 0) {
				E = 0;
				return E;
			}
			// Note: don't have to divide row by pivot first

			// Gaussian elimnation process: part 1 involves constructing the triangular
			// upper matrix
			for (i = j; i <= m; i++) { // start from j and not 1
				if (i != j) {
					// Subtract the 1/pivot scaler times rowj from the rows below the jth column
					// entry
					double[] rowj = A[j - 1].clone();
					double scaler = A[i - 1][j - 1] / A[j - 1][j - 1]; // Cij / Cjj
					for (int entry = 0; entry < A[0].length; entry++) {
						rowj[entry] *= scaler;
						if (i != j) {
							A[i - 1][entry] -= rowj[entry];
						}
					}
				}
				if (j == m) {
					// solve for the last rows equation
					double c = A[m - 1][n - 1];
					A[m - 1][n - 1] /= c;
					A[m - 1][n] /= c;
				}
			}
		}
		// Part 2: back substitute the unknown from the last equation into the equations
		// to get the other unknowns
		Matrix x = deaug();
		for (int j = m - 1; j >= 1; j--) {
			double D = 1 / A[j - 1][j - 1];
			for (int i = m; i > j; i--) {
				x.A[j - 1][0] = D * (x.A[j - 1][0] - (A[j - 1][i - 1] * x.A[i - 1][0]));
			}

		}
		augment(x);
		return E;
	}

	/**
	 * Convert the inverse of the matrix
	 * 
	 * @return flag on success of inverting
	 */
	public int invert() {
		int E = 1; // flag
		Matrix identity = getIdentity();
		augment(identity);
		// loop through the columns
		for (int j = 1; j <= m; j++) {
			double pivot = 0;
			int i, p; // row index i and pivot index p

			/*
			 * loop through the rows of the current column and find the largest magnitude to
			 * act as the pivot position. Record the value into pivot and the row index into
			 * p. The purpose of finding the largest magnitude is to reduce the error
			 * introduced by division.
			 */
			p = j;
			for (i = j; i <= n; i++) {
				if (Math.abs(A[i - 1][j - 1]) > pivot) {
					pivot = A[i - 1][j - 1];
					p = i;
				}
			}
			// interchange rows to match proper reduced row echelon form
			if (p > j) {
				swap(p, j);
			}

			// if a column vector is the 0 vector for nxn matrix then it is invertible
			// because the determinant = 0
			if (pivot == 0) {
				E = 0;
				return E;
			}
			// divide the row by the pivot value so that the pivot position is equivalent to
			// 1
			// Note: A[0].length is used instead of n since it is an augmented matrix
			for (int k = 0; k < A[0].length; k++) {
				A[j - 1][k] = A[j - 1][k] / pivot;
			}
			// perform Gauss-Jordan elimination on the other rows to zero out that column
			for (i = 1; i <= m; i++) { // i refers to the column and j refers to the row here
				if (i != j) {
					// since the pivot is one, all you have to do is subtract it scaled up from the
					// other rows corresponding column entry from the other none pivot rows
					double[] rowj = A[j - 1].clone();
					double scaler = A[i - 1][j - 1];
					for (int entry = 0; entry < A[0].length; entry++) {
						rowj[entry] *= scaler;
						A[i - 1][entry] -= rowj[entry];
					}
				}
			}
		}
		Matrix Aneg = deaug(); // get the inverse
		A = Aneg.A;
		return E;
	}

	/**
	 * Static function for inverse. Retrieves the inverse of the matrix if it exists
	 * 
	 * @param mat1
	 *            a matrix to get the inverse of
	 * @return the inverse matrix
	 */
	public static Matrix getInverse(Matrix mat1) {
		if (mat1.determinant() != 0) {
			Matrix copyM = new Matrix(mat1);
			copyM.invert();
			return copyM;
		} else {
			System.out.println("Matrix is invertible");
			return null;
		}
	}

	/**
	 * calculates the upper triangular determinant
	 * 
	 * @return a double value of the determinant of the matrix
	 */
	public double determinant() {
		int r = 0; // number of row swaps
		double det;
		Matrix copyM = new Matrix(this);

		// loop through the columns
		for (int j = 1; j <= n; j++) {
			double pivot = 0;
			int i, p; // row index i and pivot index p

			/*
			 * loop through the rows of the current column and find the largest magnitude to
			 * act as the pivot position. Record the value into pivot and the row index into
			 * p. The purpose of finding the largest magnitude is to reduce the error
			 * introduced by division.
			 */
			p = j;
			for (i = j; i <= n; i++) {
				if (Math.abs(copyM.A[i - 1][j - 1]) > pivot) {
					pivot = copyM.A[i - 1][j - 1];
					p = i;
				}
			}
			// interchange rows to form triangular matrix
			if (p > j) {
				copyM.swap(p, j);
				r++; // keep track of permutations
			}

			// if the pivot value is zero then the determinant is zero
			if (pivot == 0) {
				det = 0;
				return det;
			}

			// zero out the lower triangle of the matrix
			for (i = j; i <= m; i++) { // start from j and not 1
				if (i > j) {
					// Subtract the 1/pivot scaler times rowj from the rows below the jth column
					// entry
					double[] rowj = copyM.A[j - 1].clone();
					double scaler = copyM.A[i - 1][j - 1] / copyM.A[j - 1][j - 1]; // Cij / Cjj
					for (int entry = 0; entry < A[0].length; entry++) {
						rowj[entry] *= scaler;
						if (i != j) {
							copyM.A[i - 1][entry] -= rowj[entry];
						}
					}
				}
			}
		}
		det = Math.pow(-1, r); // handle the swaps that occurred
		for (int k = 0; k < n; k++) {
			det *= copyM.A[k][k];
		}
		return det;
	}

	/**
	 * Find the absolute sum of every row. Take the maximum of those values.
	 * 
	 * @return the max row sum
	 */
	public double infNorm() {
		double maxRowSum = 0;
		for (int i = 0; i < m; i++) {
			double rowSum = 0;
			for (int j = 0; j < n; j++) {
				rowSum += Math.abs(A[i][j]);
			}
			if (maxRowSum < rowSum) {
				maxRowSum = rowSum;
			}
		}
		return maxRowSum;
	}

	/**
	 * the infinity condition number is a measurment for checking the sensitivity of
	 * the system to errors
	 * 
	 * @return the conditioning number, or -1 if the matrix is not invertible
	 */
	public double infCond() {
		if (determinant() != 0) {
			Matrix invert = Matrix.getInverse(this);
			double aInff = infNorm();
			double invInff = invert.infNorm();
			// System.out.println(aInff); for testing purposes
			// System.out.println(invInff);
			double ans = aInff * invInff;
			return ans;
		} else {
			return -1;
		}
	}

	/**
	 * String representation of the matrix for testing purposes
	 */
	public String toString() {
		String s = "";
		for (int i = 0; i < m; i++) {
			s += "| ";
			for (int j = 0; j < n; j++) {
				s += String.format("%10f ", A[i][j]);
			}
			s += "|\n";
		}
		return s;
	}
	
	/**
	 * Crout's Method with partial pivoting Takes a matrix and decomposes it into a
	 * lower and upper triangular matrices
	 * 
	 * @return the determinant
	 */
	public double LU() {
		double accumalator, Dkj, pivot, Djk;
		double det = 1; // determinant

		augment(getIdentity()); // augment the identity matrix to keep track of swaps

		// loop through the columns
		for (int j = 1; j <= n; j++) {
			pivot = 0;
			int p = j;

			// above diagonal computations, find Dkj = Dkj - Summation i = 1 to j-1: (Dki *
			// Dij) for augmented matrix D
			for (int k = j; k <= n; k++) {
				Dkj = A[k - 1][j - 1];
				accumalator = 0;
				for (int i = 1; i <= j - 1; i++) {
					accumalator += A[k - 1][i - 1] * A[i - 1][j - 1];
				}
				Dkj -= accumalator;
				A[k - 1][j - 1] = Dkj;
				if (Math.abs(A[k - 1][j - 1]) > pivot) {
					pivot = A[k - 1][j - 1];
					p = k;
				}

			}
			// swap to bring max value to pivot
			if (p > j) {
				swap(p, j);
				det = -1 * det; // when swapping multiply determinant by negative 1
			}
			// if a column is the zero vector then exit
			if (pivot == 0) {
				det = 0;
				return det;
			}
			// below diagonal computations, find Djk = 1/Djj * (Djk - Summation i = 1 to
			// j-1: (Dji * Dik)) for augmented matrix D
			for (int k = j + 1; k <= n; k++) {
				accumalator = 0;
				Djk = A[j - 1][k - 1];
				for (int i = 1; i <= j - 1; i++) {
					accumalator += A[j - 1][i - 1] * A[i - 1][k - 1];
				}
				Djk -= accumalator;
				Djk *= 1.0 / A[j - 1][j - 1];
				A[j - 1][k - 1] = Djk;

				det *= A[j - 1][j - 1]; // compute determinant
			}
		}

		return det;
	}
	
	/**
	 * Decomposes a Matrix into its LU factors
	 * 
	 * @return An array of Matrices with the first element being the swap matrix,
	 *         the second being L, and the third being U
	 */
	public Matrix[] LUDecompose() {
		Matrix L, U, P;

		LU();

		// get swap matrix
		P = deaug();

		L = new Matrix(m, n);
		U = new Matrix(m, n);

		// get lower matrix
		for (int j = 0; j < n; j++) {
			for (int i = j; i < m; i++) {
				L.A[i][j] = A[i][j];
			}
		}

		// get upper matrix
		for (int j = 1; j < n; j++) {
			for (int i = j - 1; i >= 0; i--) {
				U.A[i][j] = A[i][j];
			}
		}
		// diagnol for u is all 1's
		for (int d = 0; d < n; d++) {
			U.A[d][d] = 1;
		}
		Matrix[] PLU = { P, L, U };
		return PLU;
	}

	/**
	 * Static Method for breaking up a matrix into its LU factorization
	 * 
	 * @param mat
	 *            the matrix to decompose
	 * @return An array of Matrices with the first element being the swap matrix,
	 *         the second being L, and the third being U
	 */
	public static Matrix[] LUDecompose(Matrix mat) {
		Matrix copyM = new Matrix(mat);
		Matrix[] PLU = copyM.LUDecompose();
		return PLU;
	}
	
	public Matrix LUsolve(Matrix b) {
		double det, summation;
		
		// compute d = Pb
		Matrix Q = new Matrix(this);
		det = Q.LU();
		System.out.println("determinant = " + det);
		Matrix P = Q.deaug();
		Matrix d = Matrix.getMult(P, b);
		// calculate intermediate
		// for k = 1 to n find dk = 1 / Qkk (dk - summation i = 1 to k - 1: Qki * di) i is row, k is column (remember d is a vector)
		for (int k = 1; k <= n; k++) {
			summation = 0;
			for (int i = 1; i <= k - 1; i++) {
				summation += Q.A[k - 1][i - 1] * d.A[i - 1][0];
			}
			d.A[k - 1][0] -= summation;
			d.A[k - 1][0] /= Q.A[k - 1][k - 1];
		}
		// for k = n down to 1 find dk = dk - summation i = k + 1 to n: Qki * di
		for (int k = n; k >= 1; k--) {
			summation = 0;
			for (int i = k + 1; i <= n; i++) {
				summation += Q.A[k - 1][i - 1] * d.A[i - 1][0];
			}
			d.A[k - 1][0] -= summation;
		}
		return d;
	}

	// ******************************Project
	// 2******************************************

	/**
	 * The trace of a matrix is the sum of its diagonal entries
	 * 
	 * @return the sum of the diagonal entries
	 */
	public double trace() {
		if (m == n) {
			double sum = 0;
			for (int i = 0; i < n; i++) {
				sum += A[i][i];
			}
			return sum;
		} else {
			System.out.println("Not a nxn matrix");
			return -1;
		}
	}

	/**
	 * Static method for obtaining the trace of a matrix.
	 * 
	 * @param mat
	 *            the matrix to calculate the trace of
	 * @return the trace of the matrix
	 */
	public static double trace(Matrix mat) {
		if (mat.m == mat.n) {
			double sum = 0;
			for (int i = 0; i < mat.n; i++) {
				sum += mat.A[i][i];
			}
			return sum;
		} else {
			System.out.println("Not a nxn matrix");
			return -1;
		}
	}

	/*
	 * the determinant of a matrix A is equal to product of eigen values 1.) get
	 * eigenvalues 2.) get determinant 3.) check
	 */

	public Polynomial LeVerrier() {
		double[] coefficients;
		int indexStart;
		// remeber the first coefficent for monic polynomials is always one
		coefficients = new double[n + 1];
		coefficients[0] = 1;
		indexStart = 1;
		

		Matrix Bn = new Matrix(this);

		// initialize entry value to be negative the trace
		double an = -1 * trace();
		//System.out.println("first an " + an);

		// add it to the solution set because it will be a coefficient
		coefficients[indexStart] = an;

		Matrix I = getIdentity();

		for (int k = n - 1; k >= 1; k--) {

			// Bk = A(Bk+1 + an*I) or multiply the original matrix by the current B one plus
			// the identity matrix scaled by an
			//System.out.println("For k " + k);
			Matrix anI = Matrix.scale(I, an);
			Bn.add(anI);
			Bn = Matrix.getMult(this, Bn);
			//System.out.println("new Bk");
			//System.out.println(Bn);

			// the new an value is: an = -trBk / (n - k + 1)
			an = (-1.0 * trace(Bn)) / (n - k + 1);
			//System.out.println("new an " + an);
			coefficients[n - k + indexStart] = an;
		}
		Polynomial p = new Polynomial(coefficients);
		return p;
	}
	
	
	public void makeUnitVector() {
		if (n != 1 && m != 1) {
			System.out.println("Not a vector");
		}else {
		double magnitude = vectorMag();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = A[i][j] / magnitude;
			}
		}
		}
	}

	public double vectorMag() {
		if (m == 1) {
			Matrix result = Matrix.getMult(this, Matrix.getTranspose(this));
			double mag = result.getEntry(1, 1);
			return Math.sqrt(mag);
		} else if (n == 1) {
			Matrix result = Matrix.getMult(Matrix.getTranspose(this), this);
			double mag = result.getEntry(1, 1);
			return Math.sqrt(mag);
		} else {
			System.out.println("Not a vector but a matrix");
			return -1;
		}
	}

	public double powerMethod(Matrix initEigenVector, double epsi, int iterationLimit) {
		if (initEigenVector.equals(new Matrix(initEigenVector.m, initEigenVector.n))) {
			System.out.println("zero vector not allowed");
			return -1;
		} else if (n != initEigenVector.m) {
			System.out.println("invalid initial vector");
			return -1;
		} else if (m != n) {
			System.out.println("Not an nxn matrix cant perform power method");
			return -1;
		} else {
			double epsilon = epsi; // acceptable error
			double m = iterationLimit; // acceptable iterations to do
			Matrix y = initEigenVector;
			int i = 0;
			double rMag, newEigen, xInvMagnitude;
			Matrix x = Matrix.getMult(this, y); // calculate the new eigenvector for next iteration

			do {
				xInvMagnitude = 1 / x.vectorMag();
				y = Matrix.scale(x, xInvMagnitude); // get unit vector of the new eigenvector (normalize it)
				x = Matrix.getMult(this, y); // calculate the new eigenvector for next iteration
				Matrix yTran = Matrix.getTranspose(y);
				// solve mu = y^T * x / y^T * y mu here aproximates the dominant eigenvalue
				Matrix numerator = Matrix.getMult(yTran, x);
				Matrix denom = Matrix.getMult(yTran, y);
				/*
				 * denom.invert(); Matrix u = Matrix.getMult(numerator, denom); // mu
				 * aproximates the dominant eigenvalue newEigen = u.getEntry(1, 1); or
				 */
				newEigen = numerator.getEntry(1, 1) / denom.getEntry(1, 1);
				Matrix r = Matrix.getDiff(Matrix.scale(y, newEigen), x); // the residual error vector
				rMag = r.vectorMag();
				i++;
			} while (rMag > epsilon && i < m);
			System.out.println("x");
			System.out.println(x);
			System.out.println("y");
			System.out.println(y);
			return newEigen;
		}
	}

	public double inversePowerMethod(Matrix initEigenVector, double epsi, int iterationLimit) {
		if (initEigenVector.equals(new Matrix(initEigenVector.m, initEigenVector.n))) {
			System.out.println("zero vector not allowed");
			return -1;
		} else if (n != initEigenVector.m) {
			System.out.println("invalid initial vector");
			return -1;
		} else if (m != n) {
			System.out.println("Not an nxn matrix cant perform power method");
			return -1;
		} else {
			double epsilon = epsi; // acceptable error
			double m = iterationLimit; // acceptable iterations to do
			Matrix y = initEigenVector;
			int i = 0;
			double rMag, newEigen, xInvMagnitude;
			Matrix x = LUsolve(y); // calculate the new eigenvector for next iteration

			do {
				xInvMagnitude = 1 / x.vectorMag();
				y = Matrix.scale(x, xInvMagnitude); // get unit vector of the new eigenvector (normalize it)
				x = LUsolve(y); // calculate the new eigenvector for next iteration
				Matrix yTran = Matrix.getTranspose(y);
				// solve mu = y^T * x / y^T * y mu here aproximates the dominant eigenvalue
				Matrix numerator = Matrix.getMult(yTran, x);
				Matrix denom = Matrix.getMult(yTran, y);
				/*
				 * denom.invert(); Matrix u = Matrix.getMult(numerator, denom); // mu
				 * aproximates the dominant eigenvalue newEigen = u.getEntry(1, 1); or
				 */
				newEigen = numerator.getEntry(1, 1) / denom.getEntry(1, 1);
				Matrix r = Matrix.getDiff(Matrix.scale(y, newEigen), x); // the residual error vector
				rMag = r.vectorMag();
				i++;
			} while (rMag > epsilon && i < m);
			System.out.println("x");
			System.out.println(x);
			System.out.println("y");
			System.out.println(y);
			return 1 / newEigen;
		}
	}
	
	
	public double inversePowerMethod2(Matrix initEigenVector, double epsi, int iterationLimit) {
	 Matrix inv = new Matrix(this);
	 inv.invert();
	 double smallEig = inv.powerMethod(initEigenVector, epsi, iterationLimit);
	 return 1 / smallEig;
	 
	}

	private double sgn(double x) {
		if (x >= 0)
			return 1;
		else
			return -1;
	}
	
	/**
	 * Method for tridiagonalizing a matrix, or a matrix in the form with its diagonal and both the first upper and lower sub-diagonals.
	 * All the other diagonals are zerod out. The matrix works by finding the the transformation matrix P who will linearly diagonalize 
	 * a matrix A through several iterations of multiplication. A(i+1) = P * A(i) * P and the matrix P is found by the equation:
	 * P = I - 2(u*u^T / u^T * u)
	 * example:
	 *        |val val val val |            |val val 0   0   |             |val val 0   0   |
	 * A(1) = |val val val val |  -> A(2) = |val val val val |   -> A(3) = |val val val 0   | finished  
	 *        |val val val val |     i = 1  |0   val val val |      i = 2  |0   val val val | i to  n - 2 iterations
	 *        |val val val val |            |0   val val val |             |0   0   val val | 
	 *                                  use minor to find next P                                             
	 */
	public Matrix houseHolder() {
		/*
		 * since we know what we want to get as a result of multiplying P we can solve for the column transformation matrix: 
		 * P(minor) * x = | alpha | = (I - 2uuT/uTu)x = x - 2(uuT/uTu)x -> solve for all of the entries (a21, a31, ...)
		 *                |   0   |
		 *                |   0   |
		 *                |   0   |
		 * Solving for alpha:
		 * alpha = (-1)(sign of entry in the alpha eq.)(summation k=2 to n (aj1)^2) ^ (1/2)
		 * Then solving for the entries of vector x one can figure out the vector u.
		 * all of the entries corresponding to the major of P(minor) are zeros and would be zeros in u.
		 * Finally all of the other entries in x are homogeneous except for the entry with alpha. From that we derive the equation
		 * Thus u^T = [0,0,...,Ak+1 + alpha, Ak+2,...,An]
		 */

		double alpha;
		Matrix u = new Matrix(m, 1);
		Matrix Q = getIdentity();
		Matrix B = new Matrix(this);
		for (int k = 1; k <= m - 2; k++) {
			// solve for alpha
			double pythag = 0;
			for (int j = k; j < m; j++) {
				double square = B.A[j][k - 1] * B.A[j][k - 1]; // same as Ak+1,k but array indexing starts at 0
				pythag += square;
			}
			pythag = Math.sqrt(pythag);
			alpha = sgn(B.A[k][k - 1]) * pythag; // remeber index starts at zero for arrays

			// now that we know alpha we can construct the u vector
			for (int i = 0; i < u.m; i++) {
				// zero out the corresponding positions to the major tridiagonal entries
				if (i < k) {
					u.setEntry(i + 1, 1, 0);
				} else if (i == k) {
					// case where its the alpha entry
					u.setEntry(i + 1, 1, B.A[i][k - 1] + alpha);
				} else {
					// other entries
					u.setEntry(i + 1, 1, B.A[i][k - 1]);
				}
			}
			// calculate P
			Matrix uTran = Matrix.getTranspose(u);
			Matrix w = Matrix.getMult(u, uTran);
			double dotProd = Matrix.getMult(uTran, u).getEntry(1, 1);
			w.scale(2.0 * (1.0 / dotProd));
			Matrix P = getIdentity();
			P.minus(w);
			Q.times(P);
			
			// calculate new A
			B = Matrix.getMult(P, B);
			B = Matrix.getMult(B, P);
		}
		return B;
	}
	
	/**
	 * The QR method involves repeated QR factorization of a Matrix A, A = QR, where Q is an orthogonal matrix and R is an upper
	 * triangular matrix.  Then the factors are multiplied in reverse order. This process is repeated so that the lower triangular
	 * entries converge to zero and the diagonal converges to the eigenvalues. 
	 */
	public double[] qr(double epsilon, int iterationLimit) {
		/*
		 * Finding the upper hesenberg form of matrix A is ideal because it makes the computation for A's eigenvalues much easier.
		 * Suppose B is the hesenberg matrix of A. Then B = Q^T*A*Q hence Q*B = Q*QT*A*Q = A*Q -> Q*B = A*Q
		 * If ek is an eigenvector of B with eigenvalue Lk, then B*ek = Lk*ek by definition of eigenvectors
		 * A*Q = Q*B -> Q*B*ek = Q*Lk*ek by substitution. Since Lk is a scalar Q*Lk*ek = Lk*Q*ek 
		 * Lk*Q*ek = Q*Lk*ek = Q*B*ek = A*Q*ek -> Lk*Q*ek = A*Q*ek
		 * This means that Lk is an eigenvalue of A and Q*ek is its corresponding eigenvector
		 */
		double c, s;
		Matrix B = houseHolder();
		int i = 0;
		boolean biggerThanEpsi; // flag for epsilon
		do {
			Matrix Qtran = getIdentity();
			for (int k = 1; k <= B.n - 1; k++) { // 1 to n - 1
				c = B.getEntry(k, k) / Math.hypot(B.getEntry(k, k), B.getEntry(k + 1, k)); // c = Bk,k / sqrt(Bk,k^2 + Bk+1,k^2)
				s = B.getEntry(k + 1, k) / Math.hypot(B.getEntry(k, k), B.getEntry(k + 1, k)); // s = Bk+1,k / sqrt(Bk,k^2 + Bk+1,k^2)

				Matrix P = getIdentity();
				// diagonal values
				P.setEntry(k, k, c);
				P.setEntry(k + 1, k + 1, c);
				// left diagonal values
				P.setEntry(k + 1, k, -1 * s);
				P.setEntry(k, k + 1, s);

				B = Matrix.getMult(P, B);
				Qtran = Matrix.getMult(P, Qtran);
			}
			Matrix Q = Matrix.getTranspose(Qtran);
			B = Matrix.getMult(B, Q);
			i++;

			// check the lower triangular values
			
			biggerThanEpsi = false;
			for (int col = 1; col <= B.n; col++) {
				for (int row = col; row <= B.m; row++) {
					if (B.getEntry(row, col) > epsilon) {
						biggerThanEpsi = true;
					}
				}
			}
		} while (biggerThanEpsi && i < iterationLimit);

		// return the eigenvalues
		double[] eigenvalues = new double[m];
		for (int h = 1; h <= B.m; h++) {
			eigenvalues[h-1] = B.getEntry(h, h);
		}
		System.out.println(B);
		return eigenvalues;
	}
	
}
