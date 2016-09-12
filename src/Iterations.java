import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;

/**
 * Created by Kris on 4/14/2016.
 */
public class Iterations {
    private int iterations;
    private boolean didSucceed;
    private Array2DRowRealMatrix A, b, eigenvector;

    public Array2DRowRealMatrix jacobi_iter(Array2DRowRealMatrix Ab, double tolerance, int maxIterations, Array2DRowRealMatrix x0) {
        iterations = 1;
        extract(Ab);
        int n = A.getRowDimension();

        Array2DRowRealMatrix X = new Array2DRowRealMatrix(n, 1);
        Array2DRowRealMatrix previousX = x0;

        boolean stop;
        do {
            for (int i = 0; i < n; i++) {
                double sum = b.getEntry(i, 0);

                for (int j = 0; j < n; j++) {
                    if (j != i)
                        sum -= A.getEntry(i, j) * previousX.getEntry(j, 0);
                }
                X.setEntry(i, 0, 1 / A.getEntry(i, i) * sum);
            }

            stop = true;
            for (int i = 0; i < n && stop; i++)
                if (Math.abs(X.getEntry(i, 0) - previousX.getEntry(i, 0)) > tolerance) {
                    stop = false;
                }

            previousX = (Array2DRowRealMatrix) X.copy();
            iterations++;
        } while (!stop && iterations < maxIterations);

        didSucceed = iterations < maxIterations; //TODO: Should this be <= instead?
        return previousX;
    }

    public Array2DRowRealMatrix gs_iter(Array2DRowRealMatrix Ab, double tolerance, int maxIterations, Array2DRowRealMatrix x0) {
        iterations = 1;
        extract(Ab);
        int n = A.getRowDimension();

        Array2DRowRealMatrix X = new Array2DRowRealMatrix(n, 1);
        Array2DRowRealMatrix previousX = x0;

        boolean stop;
        do {
            for (int i = 0; i < n; i++) {
                double sum = b.getEntry(i, 0);

                for (int j = 0; j < n; j++) {
                    if (j != i)
                        sum -= A.getEntry(i, j) * X.getEntry(j, 0);
                }
                X.setEntry(i, 0, 1 / A.getEntry(i, i) * sum);
            }


            stop = true;
            for (int i = 0; i < n && stop; i++)
                if (Math.abs(X.getEntry(i, 0) - previousX.getEntry(i, 0)) > tolerance) {
                    stop = false;
                }

            previousX = (Array2DRowRealMatrix) X.copy();
            iterations++;
        } while (!stop && iterations < maxIterations);

        didSucceed = iterations < maxIterations;
        return previousX;
    }

    public double power_method(Array2DRowRealMatrix A, ArrayRealVector u0,
                                                    ArrayRealVector w, double tolerance, int maxIterations) {

        int n = A.getRowDimension();
        Array2DRowRealMatrix u = new Array2DRowRealMatrix(n, 1);
        Array2DRowRealMatrix b = new Array2DRowRealMatrix(u0.toArray());
        Array2DRowRealMatrix uPrev = (Array2DRowRealMatrix) b.copy();
        Array2DRowRealMatrix wTranspose = new Array2DRowRealMatrix(w.toArray());
        wTranspose = (Array2DRowRealMatrix) wTranspose.transpose(); // TODO: transpose()

        boolean stop;
        iterations = 1;
        double maxEigenvalue;
        double prevMaxEigenvalue = 0;

        do {
            for(int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    u.setEntry(i, 0, u.getEntry(i, 0) + A.getEntry(i, j) * uPrev.getEntry(j, 0));
                }
            }

            double numerator = wTranspose.multiply(u).getEntry(0, 0);
            double denominator = wTranspose.multiply(uPrev).getEntry(0, 0);
            maxEigenvalue = numerator / denominator;

            stop = true;
            if (Math.abs(maxEigenvalue - prevMaxEigenvalue) > tolerance) {
                stop = false;
            }

            prevMaxEigenvalue = maxEigenvalue;
            uPrev = (Array2DRowRealMatrix) u.copy();

            double norm_sq = 0;
            for (int k = 0; k < n; k++)
                norm_sq += u.getEntry(k, 0) * u.getEntry(k, 0);

            double norm = Math.sqrt(norm_sq);
            for (int i = 0; i < n; i++) {
                u.setEntry(i, 0, u.getEntry(i, 0) / norm);
            }
            eigenvector = (Array2DRowRealMatrix) u.copy();
            iterations++;
        } while(!stop && iterations < maxIterations);

        didSucceed = iterations < maxIterations;
        return Math.abs(maxEigenvalue);
    }

    public Array2DRowRealMatrix getEigenvector() {
        return eigenvector;
    }

    public double computeError(Array2DRowRealMatrix Ab, Array2DRowRealMatrix x) {
        extract(Ab);
        return A.multiply(x).subtract(b).getNorm(); // TODO: norm()
    }

    public void extract(Array2DRowRealMatrix Ab) {
        int rows = Ab.getRowDimension();
        int cols = Ab.getColumnDimension();
        A = (Array2DRowRealMatrix) Ab.getSubMatrix(0, rows - 1, 0, cols - 2);
        b = (Array2DRowRealMatrix) Ab.getSubMatrix(0, rows - 1, cols - 1, cols - 1);
    }

    public boolean didSucceed() {
        return didSucceed;
    }

    public int getIterations() {
        return iterations;
    }

    public static void main(String[] args) {
        // power method
        double[][] APrep = {{3, 4}, {3, 1}};
        double[] u0prep = {1, 0};
        double[] wPrep = {1, 0};
        Array2DRowRealMatrix A = new Array2DRowRealMatrix(APrep);
        ArrayRealVector u0 = new ArrayRealVector(u0prep);
        ArrayRealVector w = new ArrayRealVector(wPrep);

        // jacobi and GS
        double[][] AbPrep = {{2, 4, 7}, {1, 3, 5}};
        double[][] AbPrep1 = {{5, -2, 3, -1}, {-3, 9, 1, 2}, {2, -1, -7, 3}};
        double[][] x0Prep = {{1}, {0}, {0}};

        Array2DRowRealMatrix Ab = new Array2DRowRealMatrix(AbPrep);
        Array2DRowRealMatrix Ab1 = new Array2DRowRealMatrix(AbPrep1);
        Array2DRowRealMatrix x0 = new Array2DRowRealMatrix(x0Prep);

        Iterations it = new Iterations();

        System.out.println("Power: " + it.power_method(A, u0, w, 1e-15, 100));
        System.out.println("Power eigenvector: " + it.getEigenvector());
        System.out.println("Power succeeded? " + it.didSucceed());

        Array2DRowRealMatrix x = it.jacobi_iter(Ab, 1e-10, 100, x0);

        //debug
        System.out.println();
        FileHandler.printMatrix(Ab);

        System.out.println("\nJacobi X: " + x); // .5, 1.5
        System.out.println("Jacobi error: " + it.computeError(Ab, x)); // .5, 1.5
        System.out.println("Jacobi iterations: " + it.iterations);
        System.out.println("Jacobi succeeded? " + it.didSucceed()+ "\n");

//        Array2DRowRealMatrix x1 = it.jacobi_iter(Ab1, 1e-10, 100);
//        System.out.println("Jacobi1 X: " + x1);
//        System.out.println("Jacobi1 error: " + it.computeError(Ab1, x1));
//        System.out.println("Jacobi1 iterations: " + it.iterations);
//        System.out.println("Jacobi1 succeeded? " + it.didSucceed()+ "\n");

        Array2DRowRealMatrix gsX = it.gs_iter(Ab, 1e-10, 100, x0);
        System.out.println("\nGS X: " + gsX); // .5, 1.5
        System.out.println("GS error: " + it.computeError(Ab, gsX));
        System.out.println("GS iterations: " + it.iterations);
        System.out.println("GS succeeded? " + it.didSucceed()+ "\n");

//        Array2DRowRealMatrix gsX1 = it.jacobi_iter(Ab1, 1e-10, 100);
//        System.out.println("GS1 X: " + gsX1);
//        System.out.println("GS1 error: " + it.computeError(Ab1, gsX1));
//        System.out.println("GS1 iterations: " + it.iterations);
//        System.out.println("GS1 succeeded? " + it.didSucceed()+ "\n");
    }

}
