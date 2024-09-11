import java.util.Arrays;

import org.apache.commons.math4.legacy.linear.Array2DRowRealMatrix;
import org.apache.commons.math4.legacy.linear.ArrayRealVector;
import org.apache.commons.math4.legacy.linear.MatrixUtils;
import org.apache.commons.math4.legacy.linear.RealMatrix;
import org.apache.commons.math4.legacy.linear.RealVector;



public class SQP {
	// Optimization result structure
    public static class OptimizationResult {
        public double[] solution;
        public double[] multipliers;
        public double objectiveValue;

        public OptimizationResult(double[] solution, double[] multipliers, double objectiveValue) {
            this.solution = solution;
            this.multipliers = multipliers;
            this.objectiveValue = objectiveValue;
        }
    }

    // Interface for function evaluation
    public interface FunctionEvaluator {
        double evaluate(double[] x);
    }

    private final FunctionEvaluator objectiveEvaluator;
    private final FunctionEvaluator[] constraintEvaluators;
    private double[] x;
    private double[] lambda;
    private RealMatrix hessianLagrangian;
    private double tol = 1e-6;
    private double alpha = 1.0;
    private double beta = 0.5;
    private int maxIter = 100;
    private double h = 0.001;

    // Constructor to initialize the SQP problem
    public SQP(FunctionEvaluator objectiveEvaluator, FunctionEvaluator[] constraintEvaluators, double[] initialX, double[] initialLambda) {
        this.objectiveEvaluator = objectiveEvaluator;
        this.constraintEvaluators = constraintEvaluators;
        this.x = Arrays.copyOf(initialX, initialX.length);
        this.lambda = Arrays.copyOf(initialLambda, initialLambda.length);
        this.hessianLagrangian = MatrixUtils.createRealIdentityMatrix(initialX.length);
    }

 // Main optimization method
    public OptimizationResult optimize( double[] discreteSet) {
        System.out.println("Starting optimization process...");
        int dimension = x.length;
        int optimalIteration = -1;
        double objectiveValue = 0.0;

        for (int iter = 0; iter < maxIter; iter++) {
            System.out.println("Iteration: " + iter);
            
            // **** Make Necessary Calculations **** //

            // 1. Gradient of the objective function
            double[] gradF = GradientCalculator.calculateGradient(objectiveEvaluator, x, h);
            RealVector gradientF = new ArrayRealVector(gradF);
            System.out.println("GradF: " + Arrays.toString(gradF));
            
            // 2. Gradient of constraints + value of constraints at the initial guess
            double[][] jacobianData = new double[constraintEvaluators.length][];
            double[] gValues = new double[constraintEvaluators.length];
            for (int i = 0; i < constraintEvaluators.length; i++) {
                jacobianData[i] = GradientCalculator.calculateGradient(constraintEvaluators[i], x, h);
                gValues[i] = constraintEvaluators[i].evaluate(x);
            }
            RealMatrix J = new Array2DRowRealMatrix(jacobianData);
            System.out.println("Jacobian data: " + J.toString());
            RealVector g = new ArrayRealVector(gValues);
            System.out.println("gValues: " + g.toString());
            
            // 3. Value of the objective function at the current point
            double fValue = objectiveEvaluator.evaluate(x);
            System.out.println("Objective function value  " + fValue);

            // 4. Penalty Calculation
            ViolsResult result = Violations.viols(gValues, lambda);
            double P = PenaltyCalculation.penalty(fValue, result.gViol, result.lViol);

            // **** Solve KKT Equations **** //
            RealVector results = SQP_Helpers.solveKKT(hessianLagrangian, J, gradientF, g);

            // Extract deltaX and deltaLambda from the results
            double[] deltaX = Arrays.copyOfRange(results.toArray(), 0, dimension);
            System.out.println("deltaX: " + Arrays.toString(deltaX));
            double[] deltaLambda = Arrays.copyOfRange(results.toArray(), dimension, dimension + lambda.length);

            // Update lambda using deltaLambda
            lambda = SQP_Helpers.updateLambda(lambda, deltaLambda);
            System.out.println("Lambda: " + Arrays.toString(lambda));
            
            // Check for active constraints and update deltaX appropriately
            boolean[] isActive = new boolean[lambda.length];
            Arrays.fill(isActive, true);
            for (int i = 0; i < lambda.length; i++) {
                if (lambda[i] < 0) {
                    isActive[i] = false;
                    lambda[i] = 0;
                }
            }
            int activeCount = 0;
            for (boolean active : isActive) {
                if (active) {
                    activeCount++;
                }
            }

            if (activeCount == 0) {
                results = SQP_Helpers.performCalculation(hessianLagrangian, gradientF);
            } else {
                double[][] activeJacobianData = new double[activeCount][];
                double[] activeG = new double[activeCount];
                int index = 0;
                for (int i = 0; i < lambda.length; i++) {
                    if (isActive[i]) {
                        activeJacobianData[index] = GradientCalculator.calculateGradient(constraintEvaluators[i], x, h);
                        activeG[index] = gValues[i];
                        index++;
                    }
                }
                RealMatrix activeJ = new Array2DRowRealMatrix(activeJacobianData);
                RealVector activeGVector = new ArrayRealVector(activeG);
                results = SQP_Helpers.solveKKT(hessianLagrangian, activeJ, gradientF, activeGVector);
            }

            deltaX = Arrays.copyOfRange(results.toArray(), 0, dimension);

            // Update x
            double[] xNew = SQP_Helpers.updateX(x, deltaX, alpha);
            System.out.println("xNew: " + Arrays.toString(xNew));
            
            // Rounding
             xNew = SQP_Helpers.roundVariables(xNew, discreteSet);
             System.out.println("xNew after Rounding : " + Arrays.toString(xNew));

            // Evaluate new function values
            double fEvalNew = objectiveEvaluator.evaluate(xNew);
            double[] gEvalNew = new double[constraintEvaluators.length];
            for (int i = 0; i < constraintEvaluators.length; i++) {
                gEvalNew[i] = constraintEvaluators[i].evaluate(xNew);
            }

            // Calculate penalty for the new values
            ViolsResult resultNew = Violations.viols(gEvalNew, lambda);
            double PNew = PenaltyCalculation.penalty(fEvalNew, resultNew.gViol, resultNew.lViol);

            // Reduce step size if the penalty increased
            while (PNew - P > 1e-4) {
                for (int i = 0; i < deltaX.length; i++) {
                    deltaX[i] *= beta;
                }
                xNew = SQP_Helpers.updateX(x, deltaX, alpha);
                fEvalNew = objectiveEvaluator.evaluate(xNew);
                for (int i = 0; i < constraintEvaluators.length; i++) {
                    gEvalNew[i] = constraintEvaluators[i].evaluate(xNew);
                }
                resultNew = Violations.viols(gEvalNew, lambda);
                PNew = PenaltyCalculation.penalty(fEvalNew, resultNew.gViol, resultNew.lViol);
            }

            // **** Stopping Criterion **** //
            if (StoppingCriterion.converged(xNew, x, tol)) {
                optimalIteration = iter;
                objectiveValue = fEvalNew;
                System.out.println("Stopping criterion met at iteration: " + iter);
                break;
            }

            // **** Hessian Update **** //
            double[] gradfEvalNew = GradientCalculator.calculateGradient(objectiveEvaluator, xNew, h);
            double[][] gradgEvalNew = new double[constraintEvaluators.length][];
            for (int i = 0; i < constraintEvaluators.length; i++) {
                gradgEvalNew[i] = GradientCalculator.calculateGradient(constraintEvaluators[i], xNew, h);
            }
            double[] gradLEvalNew = SQP_Helpers.calculateGradientLagrangian(gradfEvalNew, gradgEvalNew, lambda);
            double[] gradLEval = SQP_Helpers.calculateGradientLagrangian(gradF, jacobianData, lambda);

            double[] y0 = new double[gradLEval.length];
            for (int i = 0; i < gradLEval.length; i++) {
                y0[i] = gradLEvalNew[i] - gradLEval[i];
            }

            double[] hessianDeltaX = new double[x.length];
            for (int i = 0; i < x.length; i++) {
                hessianDeltaX[i] = xNew[i] - x[i];
            }

            double[][] hessianNew = BFGSUpdate.updateHessian(hessianLagrangian.getData(), y0, hessianDeltaX);
            System.out.println("hessianNew: " + Arrays.deepToString(hessianNew));
            hessianLagrangian = MatrixUtils.createRealMatrix(hessianNew);

            // **** Update all Values for Next Iteration **** //
            fValue = fEvalNew;
            gValues = gEvalNew;
            gradF = gradfEvalNew;
            jacobianData = gradgEvalNew;
            P = PNew;
            x = Arrays.copyOf(xNew, xNew.length);
            System.out.println("x: " + Arrays.toString(x));
        }
        System.out.println("Optimum found at iteration: " + optimalIteration);
        
        // Return optimization result (Optional)
        //  x = SQP_Helpers.roundVariables(x, discreteSet);
         
        // Return optimization result
        return new OptimizationResult(x, lambda, objectiveValue);
    }
    

}
