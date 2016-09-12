import org.apache.commons.math3.linear.Array2DRowRealMatrix;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;

/**
 * Created by Kris on 4/18/2016.
 */
public class Part1 {
    public static Array2DRowRealMatrix generateHilbert(int n) {
        Array2DRowRealMatrix hilbert = new Array2DRowRealMatrix(n, n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                hilbert.setEntry(i, j, 1.0 / (i + j + 1));
            }
        }
        return hilbert;
    }

    public static Array2DRowRealMatrix generateB(int n) {
        Array2DRowRealMatrix b = new Array2DRowRealMatrix(n, 1);
        for (int i = 0; i < n; i++) {
            b.setEntry(i, 0, Math.pow(0.1, n / 3.0));
        }
        return b;
    }

    public static double computeError(Array2DRowRealMatrix H, Array2DRowRealMatrix x, Array2DRowRealMatrix b) {
        return H.multiply(x).subtract(b).getNorm(); // TODO: norm()
    }

    public static void partA() {
        LU lu = new LU();
        QR qr = new QR();
        Array2DRowRealMatrix hilbert, b;
        for (int i = 2; i <= 20; i++) {
            hilbert = generateHilbert(i);
            b = generateB(i);

            Array2DRowRealMatrix luX = lu.solve_lu(hilbert, b);
            System.out.println("LUX: " + luX.getColumnVector(0).toArray());
            System.out.println("\nLU: n = " + i);
            System.out.println("x = " + luX);
            System.out.println("LU Error: " + lu.computeError(hilbert));
            System.out.println("H Error: " + computeError(hilbert, luX, b));

            Array2DRowRealMatrix qrX = qr.solve_qr_house(hilbert, b);
            System.out.println("\nQR: n = " + i);
            System.out.println("x = " + qrX);
            System.out.println("QR Error: " + qr.computeError(hilbert));
            System.out.println("H Error: " + computeError(hilbert, qrX, b));
        }
    }

    public static void Part1ToCSV() throws FileNotFoundException {
        PrintWriter pw = new PrintWriter(new File("hilbert.csv"));
        StringBuilder sb = new StringBuilder();
        sb.append("n,");
        // all n values
        for (int i = 2; i <= 20; i++) {
            sb.append(i + ",");
        }
        sb.append('\n');

        LU lu = new LU();
        QR qr = new QR();
        Array2DRowRealMatrix hilbert, b;

        sb.append("|LU-H|,");
        for (int i = 2; i <= 20; i++) {
            hilbert = Part1.generateHilbert(i);
            b = Part1.generateB(i);
            lu.solve_lu(hilbert, b);
            sb.append(lu.computeError(hilbert) + ",");
        }
        sb.append('\n');
        sb.append("Householder |QR-H|,");
        // all qr-h values
        for (int i = 2; i <= 20; i++) {
            hilbert = Part1.generateHilbert(i);
            b = Part1.generateB(i);
            qr.solve_qr_house(hilbert, b);
            sb.append(qr.computeError(hilbert) + ",");
        }
        sb.append('\n');
        sb.append("Givens |QR-H|,");
        // all qr-h values
        for (int i = 2; i <= 20; i++) {
            hilbert = Part1.generateHilbert(i);
            b = Part1.generateB(i);
            qr.solve_qr_givens(hilbert, b);
            sb.append(qr.computeError(hilbert) + ",");
        }
        sb.append('\n');
        sb.append("LU: |Hx-b|,");
        // all hx-b;
        for (int i = 2; i <= 20; i++) {
            hilbert = Part1.generateHilbert(i);
            b = Part1.generateB(i);
            Array2DRowRealMatrix x = lu.solve_lu(hilbert, b);
            sb.append(Part1.computeError(hilbert, x, b) + ",");
        }
        sb.append('\n');
        sb.append("QR Householder: |Hx-b|,");
        // all hx-b;
        for (int i = 2; i <= 20; i++) {
            hilbert = Part1.generateHilbert(i);
            b = Part1.generateB(i);
            Array2DRowRealMatrix x = qr.solve_qr_house(hilbert, b);
            sb.append(Part1.computeError(hilbert, x, b) + ",");
        }
        sb.append('\n');
        sb.append("QR Givens: |Hx-b|,");
        // all hx-b;
        for (int i = 2; i <= 20; i++) {
            hilbert = Part1.generateHilbert(i);
            b = Part1.generateB(i);
            Array2DRowRealMatrix x = qr.solve_qr_givens(hilbert, b);
            sb.append(Part1.computeError(hilbert, x, b) + ",");
        }

        pw.write(sb.toString());
        pw.close();
    }

    public static void main(String[] args) {
        partA();
        try {
            Part1ToCSV();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }
}
