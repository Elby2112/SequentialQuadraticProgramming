import java.util.Arrays;

public class Main {
	public static void main(String[] args) {
		
		        // Define the objective function
		       // Objective: f(x1, x2) = x1^4 - 2 * x2 * x1^2 + x2^2 + x1^2 - 2 * x1 + 5
		        SQP.FunctionEvaluator objectiveEvaluator = new SQP.FunctionEvaluator() {
		            @Override
		            public double evaluate(double[] x) {
		                double x1 = x[0];
		                double x2 = x[1];
		                return Math.pow(x1, 4) - 2 * x2 * Math.pow(x1, 2) + Math.pow(x2, 2) + Math.pow(x1, 2) - 2 * x1 + 5;
		            }
		        };

		        // Define the constraint function
		        // Constraint: g(x1, x2) = -(x1 + 0.25)^2 + 0.75 * x2 ≤ 0
		        SQP.FunctionEvaluator constraintEvaluator = new SQP.FunctionEvaluator() {
		            @Override
		            public double evaluate(double[] x) {
		                double x1 = x[0];
		                double x2 = x[1];
		                return -Math.pow(x1 + 0.25, 2) + 0.75 * x2;
		            }
		        };

		        // Create an array of constraints (just one in this case)
		        SQP.FunctionEvaluator[] constraintEvaluators = new SQP.FunctionEvaluator[]{constraintEvaluator};

		        // Set initial guess for variables
		        double[] initialX = {-1, 4};  // Starting point for the optimization
		        double[] initialLambda = {0.0};  // Initial guess for the Lagrange multiplier (each lagrange multiplier refers to one constraint)

		        // Define the discrete set for rounding (if needed)
		        double[] discreteSet = {0.1, 0.5, 1.0};  // Example discrete set, modify as needed

		        // Create an instance of the SQP solver
		        SQP sqp = new SQP(objectiveEvaluator, constraintEvaluators, initialX, initialLambda);

		        // Run the optimization process
		        SQP.OptimizationResult result = sqp.optimize(discreteSet);

		        // Print the optimization result
		        System.out.println("Optimal solution: " + Arrays.toString(result.solution));
		        System.out.println("Lagrange multipliers: " + Arrays.toString(result.multipliers));
		        System.out.println("Objective value at optimal solution: " + result.objectiveValue);
		    }
		}


