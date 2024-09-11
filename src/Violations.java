

public class Violations {
	// Method to evaluate violations
		public static ViolsResult viols(double[] gEval ,double[] l) {
			 ViolsResult result = new ViolsResult(gEval.length);

		     for (int i = 0; i < gEval.length; i++) {
		         if (gEval[i] < 0) {
		             if (result.count < gEval.length) { // Check to prevent array overflow
		                 result.gViol[result.count] = gEval[i];
		                 result.lViol[result.count] = l[i];
		                 result.count++;
		             }
		         }
		     }

		     return result;
		}
	 // Main method to test the functionality
	 public static void main(String[] args) {
	     double[] gEval = {-5.0, 3.0, -2.0};
	     double[] l = {1.0, 2.0, 3.0};

	     ViolsResult result = viols(gEval, l);

	     System.out.println("Violations found: " + result.count);
	     for (int i = 0; i < result.count; i++) {
	         System.out.println("gViol[" + i + "]: " + result.gViol[i]);
	         System.out.println("lViol[" + i + "]: " + result.lViol[i]);
	     }
	 }

}
