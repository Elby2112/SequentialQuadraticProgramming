
public class PenaltyCalculation {
	// Method to calculate the penalty
    public static double penalty(double fEval, double[] gViol, double[] lViol) {
        double sum = 0;
        for (int i = 0; i < gViol.length; i++) {
            sum += lViol[i] * Math.abs(gViol[i]);
        }
        return fEval + sum;
    }

    // Main method to test the functionality
    public static void main(String[] args) {
        double fEval = 100;  // Example value for fEval
        double[] gViol = {-5.0, 3.0, -2.0}; // Example gViol values
        double[] lViol = {1.0, 2.0, 3.0};   // Example lViol values

        double result = penalty(fEval, gViol, lViol);
        System.out.println("Penalty result: " + result);
    }
}
