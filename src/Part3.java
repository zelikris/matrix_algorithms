import com.sun.org.apache.xpath.internal.SourceTree;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Random;

/**
 * Created by Kris on 4/19/2016.
 */
public class Part3 {
    private static Array2DRowRealMatrix[] matrices = new Array2DRowRealMatrix[200];
    private static Array2DRowRealMatrix[] inverses = new Array2DRowRealMatrix[200];
    private static double[] determinants = new double[200];
    private static double[] traces = new double[200];
    private static double[] Ns = new double[200];
    private static double[] inverseNs = new double[200];

    private static void generateMatrices() {
        Array2DRowRealMatrix matrix;
        Random rand = new Random();
        for (int i = 0; i < 200; i++) {
            matrix = new Array2DRowRealMatrix(2, 2);
            for (int j = 0; j < 2; j++) {
                for (int k = 0; k < 2; k++) {
                    matrix.setEntry(j, k, rand.nextDouble() * 4 - 2);
                }
            }
            matrices[i] = matrix;
        }
    }

    public static Array2DRowRealMatrix inverse(Array2DRowRealMatrix A) {
        Array2DRowRealMatrix result = (Array2DRowRealMatrix) A.copy();
        double denominator = result.getEntry(0, 0) * result.getEntry(1, 1)
                                - result.getEntry(0, 1) * result.getEntry(1, 0);
        double multiplier = 1 / denominator;
        double a = result.getEntry(0, 0);
        result.setEntry(0, 0, result.getEntry(1, 1));
        result.setEntry(1, 1, a);
        result.setEntry(0, 1, result.getEntry(0, 1) * -1);
        result.setEntry(1, 0, result.getEntry(1, 0) * -1);
        result = (Array2DRowRealMatrix) result.scalarMultiply(multiplier);
        return result;
    }

    private static void part1ToCSV() throws FileNotFoundException {
        PrintWriter pw = new PrintWriter(new File("power.csv"));
        StringBuilder sb = new StringBuilder();
        sb.append("Determinant,");
        for (int i = 0; i < matrices.length; i++) {
            sb.append(determinants[i] + ",");
        }
        sb.append('\n');

        sb.append("Trace,");
        for (int i = 0; i < matrices.length; i++) {
            sb.append(traces[i] + ",");
        }
        sb.append('\n');

        sb.append("# of Iterations for A,");
        for (int i = 0; i < matrices.length; i++) {
            sb.append(Ns[i] + ",");
        }
        sb.append('\n');
        sb.append("# of Iterations for A^-1,");
        for (int i = 0; i < matrices.length; i++) {
            sb.append(inverseNs[i] + ",");
        }

        pw.write(sb.toString());
        pw.close();
    }

    public static void main(String[] args) {
        // part A
        generateMatrices();
        for (int i = 0; i < matrices.length; i++) {
            inverses[i] = inverse(matrices[i]);
            System.out.println("\noriginal: " + matrices[i]);
            System.out.println("inverse: " + inverses[i]);
        }
        double[] u0Prep = {1, 0};
        double[] u0Prep1 = {0, 1};
        double[] u0Prep2 = {1, 1};
        ArrayRealVector u0 = new ArrayRealVector(u0Prep);
        ArrayRealVector u01 = new ArrayRealVector(u0Prep1);
        ArrayRealVector u02 = new ArrayRealVector(u0Prep2);

        ArrayRealVector w0 = u0.copy();
        ArrayRealVector w01 = u01.copy();
        ArrayRealVector w02 = u02.copy();

        Iterations it = new Iterations();

        for (int i = 0; i < matrices.length; i++) {
            double maxEigen = it.power_method(matrices[i],u0, w0, 0.00005, 100);
            System.out.println("\nMax Result: " + maxEigen);
            System.out.println("Max M: " + it.getIterations());
            System.out.println("Succeeded? " + it.didSucceed());
            if (!it.didSucceed()) {
                System.out.println("FAILED 1");
                maxEigen = it.power_method(matrices[i], u01, w01, 0.00005, 100);
                System.out.println("\nResult: " + maxEigen);
                System.out.println("M: " + it.getIterations());
                System.out.println("Succeeded? " + it.didSucceed());

                if (!it.didSucceed()) {
                    System.out.println("FAILED 2");
                    maxEigen = it.power_method(matrices[i], u02, w02, 0.00005, 100);
                    System.out.println("\nResult: " + maxEigen);
                    System.out.println("M: " + it.getIterations());
                    System.out.println("Succeeded? " + it.didSucceed());
                }
            }
            Ns[i] = it.getIterations();

            double minEigen = it.power_method(inverses[i], u0, w0, 0.00005, 100);
            System.out.println("\nMin Result: " + minEigen);
            System.out.println("Min M: " + it.getIterations());
            System.out.println("Succeeded? " + it.didSucceed());
            if (!it.didSucceed()) {
                System.out.println("FAILED 1");
                minEigen = it.power_method(inverses[i], u01, w01, 0.00005, 100);
                System.out.println("\nResult: " + minEigen);
                System.out.println("M: " + it.getIterations());
                System.out.println("Succeeded? " + it.didSucceed());

                if (!it.didSucceed()) {
                    System.out.println("FAILED 2");
                    minEigen = it.power_method(inverses[i], u02, w02, 0.00005, 100);
                    System.out.println("\nResult: " + minEigen);
                    System.out.println("M: " + it.getIterations());
                    System.out.println("Succeeded? " + it.didSucceed());
                }
            }
            inverseNs[i] = it.getIterations();

            minEigen = 1 / minEigen;
            double trace = maxEigen + minEigen;
            double determinant = maxEigen * minEigen;

            traces[i] = trace;
            determinants[i] = determinant;
//            System.out.println("A = " + matrices[i]);
            System.out.println("Max: " + maxEigen);
            System.out.println("Min: " + minEigen);
            System.out.println("Trace: " + trace);
            System.out.println("Determinant: " + determinant);
            System.out.println("===================");
        }
        try {
            part1ToCSV();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }
}
