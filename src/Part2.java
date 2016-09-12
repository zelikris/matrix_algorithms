import org.apache.commons.math3.linear.Array2DRowRealMatrix;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Random;

/**
 * Created by Kris on 4/18/2016.
 */
public class Part2 {
    private static Array2DRowRealMatrix[] x0 = new Array2DRowRealMatrix[100];
    private static int[] numIterationsJacobi, numIterationsGS;
    private static double[] initialError;
    private static double[] jacobiSuccessiveErrors, gsSuccessiveErrors;

    public static void generateX0() {
        Random rand = new Random();
        for (int i = 0; i < 100; i++) {
            x0[i] = new Array2DRowRealMatrix(3, 1);
            for (int j = 0; j < 3; j++) {
                x0[i].setEntry(j, 0, rand.nextDouble() * 20 - 10);
            }
        }
    }

    private static void Part2ToCSV() throws FileNotFoundException {

        PrintWriter pw = new PrintWriter(new File("iterations.csv"));
        StringBuilder sb = new StringBuilder();
        sb.append("|x_0 - x_exact|,");
        for (int i = 0; i < 100; i++) {
            sb.append(initialError[i] + ",");
        }
        sb.append('\n');

        sb.append("Jacobi # of Iterations,");
        for (int i = 0; i < 100; i++) {
            sb.append(numIterationsJacobi[i] + ",");
        }
        sb.append('\n');

        sb.append("GS # of Iterations,");
        for (int i = 0; i < 100; i++) {
            sb.append(numIterationsGS[i] + ",");
        }

        pw.write(sb.toString());
        pw.close();
    }

    private static void part2ToCSV() throws FileNotFoundException {
        PrintWriter pw = new PrintWriter(new File("iterations1.csv"));
        StringBuilder sb = new StringBuilder();
        sb.append("Iteration #,");
        for (int i = 2; i <= 10; i++) {
            sb.append(i + ",");
        }
        sb.append('\n');

        sb.append("Jacobi Successive Errors,");
        for (int i = 0; i < 9; i++) {
            sb.append(jacobiSuccessiveErrors[i] + ",");
        }
        sb.append('\n');

        sb.append("GS Successive Errors,");
        for (int i = 0; i < 9; i++) {
            sb.append(gsSuccessiveErrors[i] + ",");
        }

        pw.write(sb.toString());
        pw.close();
    }

    public static void main(String[] args) {
        // part a
        generateX0();
        double[][] AbPrep = {{1.0, 1.0 / 3.0, 1.0 / 9.0, 0.9}, {1.0 / 3.0, 1.0, 1.0 / 3.0, 0.1}, {1.0 / 9.0, 1.0 / 3.0, 1.0, 0.3}};
        Array2DRowRealMatrix Ab = new Array2DRowRealMatrix(AbPrep);
        Iterations it = new Iterations();

        Array2DRowRealMatrix[] jacobiResults = new Array2DRowRealMatrix[100];
        Array2DRowRealMatrix[] gsResults = new Array2DRowRealMatrix[100];
        numIterationsJacobi = new int[100];
        numIterationsGS = new int[100];
        for (int i = 0; i < 100; i++) {
            Array2DRowRealMatrix x = it.jacobi_iter(Ab, 0.00005, 100, x0[i]);
            jacobiResults[i] = x;
            System.out.println(x);
            System.out.println(it.didSucceed());
            System.out.println(it.getIterations());
            numIterationsJacobi[i] = it.getIterations();
        }
        for (int i = 0; i < 100; i++) {
            Array2DRowRealMatrix x = it.gs_iter(Ab, 0.00005, 100, x0[i]);
            gsResults[i] = x;
            System.out.println(x);
            System.out.println(it.didSucceed());
            System.out.println(it.getIterations());
            numIterationsGS[i] = it.getIterations();
        }

        // part b
        Array2DRowRealMatrix jacobiAvg = new Array2DRowRealMatrix(3, 1);
        Array2DRowRealMatrix gsAvg = new Array2DRowRealMatrix(3, 1);

        for(int i = 0; i < 100; i++) {
            jacobiAvg = jacobiAvg.add(jacobiResults[i]);
        }
        jacobiAvg = (Array2DRowRealMatrix)jacobiAvg.scalarMultiply(1.0 / 100.0);
        double[][] exactPrep = {{39.0 / 40.0}, {-13.0 / 40.0}, {12.0 / 40.0}};
        Array2DRowRealMatrix exact = new Array2DRowRealMatrix((exactPrep));
        double jacobiError = jacobiAvg.subtract(exact).getNorm();
        System.out.println("Jacobi Error: " + jacobiError);

        for(int i = 0; i < 100; i++) {
            gsAvg = gsAvg.add(gsResults[i]);
        }
        gsAvg = (Array2DRowRealMatrix)gsAvg.scalarMultiply(1.0 / 100.0);
        double gsError = gsAvg.subtract(exact).getNorm();
        System.out.println("GS Error: " + gsError);

        // part c
        double[] ratios = new double[100];
        for (int i = 0; i < 100; i++) {
            ratios[i] = numIterationsJacobi[i] / numIterationsGS[i];
        }
        double ratioAvg = 0;
        for (int i = 0; i < 100; i++) {
            ratioAvg += ratios[i];
        }
        ratioAvg = ratioAvg / 100;
        System.out.println("ratio avg: " + ratioAvg);

        // part d
        initialError = new double[100];
        for (int i = 0; i < 100; i++) {
            initialError[i] = x0[i].subtract(exact).getNorm();
        }

        try {
            Part2ToCSV();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }

        // part f
        double[][] badAbPrep = {{1, 2, 3}, {2, 1, 3}};
        double[][] badx0Prep = {{0}, {0}};
        Array2DRowRealMatrix badAb = new Array2DRowRealMatrix(badAbPrep);
        Array2DRowRealMatrix badX0 = new Array2DRowRealMatrix(badx0Prep);

        // jacobi
        Array2DRowRealMatrix jacobiResult;
        jacobiSuccessiveErrors = new double[9];
        for (int i = 2; i <= 10; i++) {
            jacobiResult = it.jacobi_iter(badAb, 0.00005, i, badX0);
            System.out.println("\nM: " + i);
            double error = it.computeError(badAb, jacobiResult);
            System.out.println("Jacobi Bad Error: " + error);
            jacobiSuccessiveErrors[i - 2] = error;
        }

        // GS
        gsSuccessiveErrors = new double[9];
        badX0 = new Array2DRowRealMatrix(badx0Prep);
        Array2DRowRealMatrix gsResult;
        for (int i = 2; i <= 10; i++) {
            gsResult = it.gs_iter(badAb, 0.00005, i, badX0);
            System.out.println("\nM: " + i);
            double error = it.computeError(badAb, gsResult);
            System.out.println("GS Bad Error: " + error);
            gsSuccessiveErrors[i - 2] = error;
        }
        try {
            part2ToCSV();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }
}
