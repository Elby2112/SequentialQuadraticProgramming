
public class StoppingCriterion {
    public static void main(String[] args) {
        // Example initialization for demonstration
        double[] xNew = {1.0, 2.0, 3.0}; // New values
        double[] x = {0.98, 2.01, 2.99}; // Previous values

        // Looping to check the condition
        while (true) {
            if (converged(xNew, x, 1e-2)) {
                break;
            }
        }
    }

 
    
    public static boolean converged(double[] xNew, double[] x, double tolerance) {
        double sum = 0;
        for (int i = 0; i < x.length; i++) {
            double diff = xNew[i] - x[i];
            sum += diff * diff;
        }
        double norm = Math.sqrt(sum);
        System.out.println("Convergence norm: " + norm);
        return norm <= tolerance;
    }
}
