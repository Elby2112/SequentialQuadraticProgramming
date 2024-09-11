

public class GradientCalculator {
	 public static double[] calculateGradient(SQP.FunctionEvaluator evaluator, double[] x, double h) {
	        double[] gradient = new double[x.length];
	        double[] xPlusH = new double[x.length];
	        double[] xMinusH = new double[x.length];
	        System.arraycopy(x, 0, xPlusH, 0, x.length);
	        System.arraycopy(x, 0, xMinusH, 0, x.length);

	        for (int i = 0; i < x.length; i++) {
	            xPlusH[i] += h;
	            xMinusH[i] -= h;

	            double fxPlusH = evaluator.evaluate(xPlusH);
	            double fxMinusH = evaluator.evaluate(xMinusH);
	            gradient[i] = (fxPlusH - fxMinusH) / (2 * h);

	            xPlusH[i] = x[i];
	            xMinusH[i] = x[i];
	        }

	        return gradient;
	    }
}
