
public class BFGSUpdate {
	 public static double[][] updateHessian(double[][] H_old, double[] y0, double[] deltaX) {
	        int n = deltaX.length;
	        double[][] H_new = new double[n][n];

	        // Calculate the denominators
	        double denom1 = dotProduct(y0, deltaX);
	        double denom2 = dotProduct(deltaX, matrixVectorMultiply(H_old, deltaX));

	        // Calculate the first correction term (y0 * y0^T) / (y0^T * deltaX)
	        double[][] firstCorrection = scalarMultiply(outerProduct(y0, y0), 1 / denom1);

	        // Calculate the second correction term (H_old * deltaX * deltaX^T * H_old) / (deltaX^T * H_old * deltaX)
	        double[][] deltaX_outer = outerProduct(deltaX, deltaX);
	        double[][] temp = matrixMultiply(H_old, deltaX_outer);
	        double[][] secondCorrection = matrixMultiply(temp, H_old);
	        secondCorrection = scalarMultiply(secondCorrection, 1 / denom2);

	        // Update Hessian matrix
	        H_new = matrixAdd(H_old, firstCorrection);
	        H_new = matrixSubtract(H_new, secondCorrection);

	        return H_new;
	    }

	    // Method to multiply two matrices
	    private static double[][] matrixMultiply(double[][] A, double[][] B) {
	        int rowsA = A.length;
	        int colsA = A[0].length;
	        int colsB = B[0].length;

	        double[][] result = new double[rowsA][colsB];
	        for (int i = 0; i < rowsA; i++) {
	            for (int j = 0; j < colsB; j++) {
	                for (int k = 0; k < colsA; k++) {
	                    result[i][j] += A[i][k] * B[k][j];
	                }
	            }
	        }
	        return result;
	    }

// Utility method to perform matrix-vector multiplication
private static double[] matrixVectorMultiply(double[][] matrix, double[] vector) {
   double[] result = new double[vector.length];
   for (int i = 0; i < matrix.length; i++) {
       result[i] = 0;
       for (int j = 0; j < vector.length; j++) {
           result[i] += matrix[i][j] * vector[j];
       }
   }
   return result;
}

// Utility method to compute the dot product of two vectors
private static double dotProduct(double[] a, double[] b) {
   double result = 0;
   for (int i = 0; i < a.length; i++) {
       result += a[i] * b[i];
   }
   return result;
}

// Utility method to perform outer product of two vectors
private static double[][] outerProduct(double[] a, double[] b) {
   double[][] result = new double[a.length][b.length];
   for (int i = 0; i < a.length; i++) {
       for (int j = 0; j < b.length; j++) {
           result[i][j] = a[i] * b[j];
       }
   }
   return result;
}

// Utility method to multiply a matrix by a scalar
private static double[][] scalarMultiply(double[][] matrix, double scalar) {
   double[][] result = new double[matrix.length][matrix[0].length];
   for (int i = 0; i < matrix.length; i++) {
       for (int j = 0; j < matrix[i].length; j++) {
           result[i][j] = matrix[i][j] * scalar;
       }
   }
   return result;
}

// Utility method to add two matrices
private static double[][] matrixAdd(double[][] A, double[][] B) {
   double[][] result = new double[A.length][A[0].length];
   for (int i = 0; i < A.length; i++) {
       for (int j = 0; j < A[i].length; j++) {
           result[i][j] = A[i][j] + B[i][j];
       }
   }
   return result;
}

// Utility method to subtract one matrix from another
private static double[][] matrixSubtract(double[][] A, double[][] B) {
   double[][] result = new double[A.length][A[0].length];
   for (int i = 0; i < A.length; i++) {
       for (int j = 0; j < A[i].length; j++) {
           result[i][j] = A[i][j] - B[i][j];
       }
   }
   return result;
}

// Example usage of the BFGS Hessian update function
public static void main(String[] args) {
   double[][] H = {{1, 0}, {0, 1}}; // Initial Hessian approximation (identity matrix)
   double[] s = {-0.5, -2.25}; // Change in variables (deltaX)
   double[] y = {-21, -7}; // Change in gradient of Lagrangian (y0)

   double[][] H_new = updateHessian(H, y, s);
   // Print the new Hessian matrix
   for (double[] row : H_new) {
       for (double value : row) {
           System.out.print(value + " ");
       }
       System.out.println();
   }
}
}
