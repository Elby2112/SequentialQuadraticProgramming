import java.util.Arrays;

import org.apache.commons.math4.legacy.linear.DecompositionSolver;
import org.apache.commons.math4.legacy.linear.LUDecomposition;
import org.apache.commons.math4.legacy.linear.MatrixUtils;
import org.apache.commons.math4.legacy.linear.RealMatrix;
import org.apache.commons.math4.legacy.linear.RealVector;

public class SQP_Helpers {
	//1- Function to print a Matrix
		 public static void printMatrix(RealMatrix matrix) {
		        System.out.println("Matrix (" + matrix.getRowDimension() + "x" + matrix.getColumnDimension() + "):");
		        for (int i = 0; i < matrix.getRowDimension(); i++) {
		            for (int j = 0; j < matrix.getColumnDimension(); j++) {
		                System.out.print(matrix.getEntry(i, j) + " ");
		            }
		            System.out.println();  
		        }
		    }
		 
		 //2- Calculation of the Gradient of Lagrangian Function where GradL = GradF + Sum(lambda(i)*GradGi)_i=1......N, where N= number of constraints
		 public static double[] calculateGradientLagrangian(double[] gradF, double[][] gradGs, double[] lambdas) {
		        double[] gradL = new double[gradF.length]; // Assumes gradF and each gradG are the same size

		        // Start with the gradient of the objective function
		        for (int i = 0; i < gradF.length; i++) {
		            gradL[i] = gradF[i];
		        }

		        // Subtract each lambda * corresponding gradient of the constraints
		        for (int j = 0; j < lambdas.length; j++) {
		            double lambda = lambdas[j];
		            double[] gradG = gradGs[j];
		            for (int i = 0; i < gradL.length; i++) {
		                gradL[i] -= lambda * gradG[i];
		            }
		        }

		        return gradL;
		    }
		 
		 //3-KKT Solver
		 public static RealVector solveKKT(RealMatrix H, RealMatrix J, RealVector gradientF, RealVector g) {
		        // Construct KKT matrix A and vector B
		        RealMatrix A = constructKKTMatrix(H, J);
		        RealVector B = constructKKTVector(gradientF, g, J);

		        // Solve the KKT system A * x = B
		        DecompositionSolver solver = new LUDecomposition(A).getSolver();
		        return solver.solve(B);
		    }

		public static RealMatrix constructKKTMatrix(RealMatrix H, RealMatrix J) {
		        int dimH = H.getRowDimension();
		        int dimJ = J.getRowDimension();
		        int dimJT = J.getColumnDimension();
		        
		        // Create a large enough matrix to hold both H and J^T, J and zero matrix
		        RealMatrix KKTMatrix = MatrixUtils.createRealMatrix(dimH + dimJ, dimH + dimJ);
		        
		        // Set H in the top left
		        KKTMatrix.setSubMatrix(H.getData(), 0, 0);
		        
		        // Set J^T in the top right
		        RealMatrix JT = J.transpose();
		        KKTMatrix.setSubMatrix(JT.getData(), 0, dimH);
		        
		        // Set J in the bottom left
		        KKTMatrix.setSubMatrix(J.getData(), dimH, 0);
		        
		        // The bottom right remains a zero matrix (already initialized as zero)
		        
		        return KKTMatrix;
		    }

		    public static RealVector constructKKTVector(RealVector gradientF, RealVector g, RealMatrix J) {
		        // Construct the KKT vector using gradients and constraint evaluations
		      //  RealVector lambdaVector = new ArrayRealVector(lambda);
		        RealVector BUpper = gradientF.mapMultiply(-1);
		        RealVector BLower = g.mapMultiply(-1);
		        return BUpper.append(BLower);
		    }
		    
		    //4- Dot Product Calculations
		    
		    private static double dotProduct(double[] a, double[] b) {
		        double result = 0;
		        for (int i = 0; i < a.length; i++) {
		            result += a[i] * b[i];
		        }
		        return result;
		    }
		    
		    //5- Another helper calculations

		    private static double[] matrixVectorMultiply(double[][] matrix, double[] vector) {
		        double[] result = new double[vector.length];
		        for (int i = 0; i < matrix.length; i++) {
		            for (int j = 0; j < vector.length; j++) {
		                result[i] += matrix[i][j] * vector[j];
		            }
		        }
		        return result;
		    }
		    
		    //6- Penality Functions
		    public static double penaltyFunction(double[] x, double[] lambda, double[] gradF, double[] gValues, double fvalue) {
				double penalty = fvalue;
			    for (int i = 0; i < gValues.length; i++) {
			        // Adding a quadratic penalty for constraint violations
			        penalty += lambda[i] * gValues[i] * gValues[i];  // Note: using gValues[i] squared to penalize the violation quadratically
			    }
			    return penalty;
			}

			public static double norm(double[] vec) {
			    double sum = 0;
			    for (double v : vec) {
			        sum += v * v;
			    }
			    return Math.sqrt(sum);
			}
			
			//7- Update Lambda and DeltaX
			
			public static double[] updateX(double[] x, double[] deltaX, double alpha) {
			    double[] newX = new double[x.length];
			    for (int i = 0; i < x.length; i++) {
			        newX[i] = x[i] + alpha * deltaX[i];
			    }
			    return newX;
			}
			
			public static double[] updateLambda(double[] lambda, double[] deltaLambda) {
			    if (lambda.length != deltaLambda.length) {
			        throw new IllegalArgumentException("The length of lambda and deltaLambda must be the same.");
			    }
			    
			    for (int i = 0; i < lambda.length; i++) {
			        lambda[i] = -(deltaLambda[i]);
			    }
			     return lambda;
			}
			
			//8- PreformCalculation are th KKT solver when all the constraints are dropped off if their lagrange multiplier is Negative 
			public static RealVector performCalculation(RealMatrix hessian, RealVector gradient) {
				// Negate the gradient
		        RealVector negatedGradient = gradient.mapMultiply(-1.0);

		        // Solve Hessian * x = -gradient
		        DecompositionSolver solver = new LUDecomposition(hessian).getSolver();
		        RealVector solution = solver.solve(negatedGradient);

		        return solution;
		    }
			
			//9- Rounding to Nearest commericial diamater set available 
			public static double roundToNearest(double value, double[] discreteSet) {
				// Sort the discrete set to ensure proper ordering
		        Arrays.sort(discreteSet);
		        
		        for (double discreteValue : discreteSet) {
		            if (discreteValue >= value) {
		                return discreteValue;
		            }
		        }
		        
		        // If no value is greater than or equal to the given value, return the greatest value
		        return discreteSet[discreteSet.length - 1];
			}
			
			public static double[] roundVariables(double[] values, double[] discreteSet) {
			    double[] roundedValues = new double[values.length];
			    for (int i = 0; i < values.length; i++) {
			        roundedValues[i] = roundToNearest(values[i], discreteSet);
			    }
			    return roundedValues;
			}


}
