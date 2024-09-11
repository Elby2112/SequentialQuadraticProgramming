
public class ViolsResult {
	 public double[] gViol;
	 public double[] lViol;
	 public int count; // To keep track of the number of violations


	public ViolsResult(int size) {
		 gViol = new double[size];
	     lViol = new double[size];
	     count = 0;
	}
}
