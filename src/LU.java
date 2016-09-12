import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.MatrixUtils;

import static java.util.Arrays.stream;

/**
 * Created by Kris on 3/15/2016.
 */
public class LU {
    private Array2DRowRealMatrix L;
    private Array2DRowRealMatrix U;
    private int[] bPivot;

    public void lu_fact(Array2DRowRealMatrix A) {
        int n = A.getColumnDimension();
        L = new Array2DRowRealMatrix(n, n);
        U = new Array2DRowRealMatrix(n, n);
        Array2DRowRealMatrix A2 = pivotise(A).multiply(A);

        for (int j = 0; j < n; j++) {
            L.setEntry(j, j, 1.0);
            for (int i = 0; i < j + 1; i++) {
                double s1 = 0;
                for (int k = 0; k < i; k++)
                    s1 += U.getEntry(k, j) * (L.getEntry(i, k));

                U.setEntry(i, j, A2.getEntry(i, j) - s1);
            }
            for (int i = j; i < n; i++) {
                double s2 = 0;
                for (int k = 0; k < j; k++)
                    s2 += (U.getEntry(k, j)*(L.getEntry(i, k)));

                double setL = (A2.getEntry(i, j) - s2) / (U.getEntry(j, j));
                L.setEntry(i, j, setL);
            }
        }
    }

    public Array2DRowRealMatrix solve_lu(Array2DRowRealMatrix Ab) {
        int rows = Ab.getRowDimension();
        int cols = Ab.getColumnDimension();
        Array2DRowRealMatrix A = (Array2DRowRealMatrix)Ab.getSubMatrix(0, rows - 1, 0, cols - 2);
        Array2DRowRealMatrix b = (Array2DRowRealMatrix)Ab.getSubMatrix(0, rows - 1, cols - 1, cols - 1);

        lu_fact(A);
        int[] bCol = {0};
        b = (Array2DRowRealMatrix) b.getSubMatrix(bPivot, bCol);
        return substitution(b);
    }

    public Array2DRowRealMatrix solve_lu(Array2DRowRealMatrix A, Array2DRowRealMatrix b) {
        lu_fact(A);
        int[] bCol = {0};
        b = (Array2DRowRealMatrix) b.getSubMatrix(bPivot, bCol);
        return substitution(b);
    }

    private Array2DRowRealMatrix substitution(Array2DRowRealMatrix b) {
        Array2DRowRealMatrix y = forwardSub(b);
        return backwardSub(y);
    }

    private Array2DRowRealMatrix forwardSub(Array2DRowRealMatrix b) {
        int n = b.getRowDimension();
        Array2DRowRealMatrix y = new Array2DRowRealMatrix(n, 1);

        // forward substitution for Ly = b
        for (int i = 0; i < n; i++) {
            y.setEntry(i, 0, b.getEntry(i, 0));
            for (int j = 0; j < i; j++) {
                y.setEntry(i, 0, y.getEntry(i, 0) - (L.getEntry(i, j) * y.getEntry(j, 0)));
            }
            y.setEntry(i, 0, y.getEntry(i, 0) / L.getEntry(i , i));
        }
        return y;
    }

    private Array2DRowRealMatrix backwardSub(Array2DRowRealMatrix y) {
        int n = y.getRowDimension();
        Array2DRowRealMatrix x = new Array2DRowRealMatrix(n, 1);

        // backward substitution for Ux = y
        for (int i = n - 1; i >= 0; i--) {
            x.setEntry(i, 0, y.getEntry(i, 0));
            for (int j = i + 1; j < n; j++) {
                x.setEntry(i, 0, x.getEntry(i , 0) - (U.getEntry(i, j) * x.getEntry(j, 0)));
            }
            x.setEntry(i, 0, x.getEntry(i, 0) / U.getEntry(i, i));
        }
        return x;
    }

    private Array2DRowRealMatrix pivotise(Array2DRowRealMatrix A) {
        int n = A.getColumnDimension();
        Array2DRowRealMatrix id = (Array2DRowRealMatrix)MatrixUtils.createRealIdentityMatrix(n);
        bPivot = new int[n];
        for (int i = 0; i < n; i++) {
            bPivot[i] = i;
        }

        for (int i = 0; i < n; i++) {
            double maxRowEntry = A.getEntry(i, i);
            int row = i;
            for (int j = i; j < n; j++) {
                if (A.getEntry(j, i) > maxRowEntry) {
                    maxRowEntry = A.getEntry(j, i);
                    row = j;
                }
            }

            if (i != row) {
                double[] aTmp = id.getRow(i);
                id.setRow(i, id.getRow(row));
                id.setRow(row, aTmp);

                // keep track of pivots so that matrix b can be rearranged appropriately
                int bTmp = bPivot[i];
                bPivot[i] = bPivot[row];
                bPivot[row] = bTmp;
            }
        }
        return id;
    }

    public double computeError(Array2DRowRealMatrix A) { // PA = LUP
        Array2DRowRealMatrix A2 = pivotise(A).multiply(A);
        return L.multiply(U).subtract(A2).getNorm(); // TODO: norm()
    }

    public Array2DRowRealMatrix getL() {
        return L;
    }

    public void setL(Array2DRowRealMatrix l) {
        L = l;
    }

    public Array2DRowRealMatrix getU() {
        return U;
    }

    public void setU(Array2DRowRealMatrix u) {
        U = u;
    }

    //for testing
    public static void main(String[] args) {

        double[][] AbPrep = {{1.0, 3, 5}, {2.0, 4, 7}};
        double[][] AbPrep1 = {{2, 1, 1, 8}, {4, -6, 0, 2}, {2, 7, 2, 3}};
        Array2DRowRealMatrix Ab = new Array2DRowRealMatrix(AbPrep);
        Array2DRowRealMatrix Ab1 = new Array2DRowRealMatrix(AbPrep1);

        // LU SOLUTION TESTS
        LU luAb = new LU();
        Array2DRowRealMatrix solution = luAb.solve_lu(Ab);
        System.out.println("X:  " + solution);

        System.out.println("X1: " + luAb.solve_lu(Ab1));



        // LU TESTS
        double[][] aPrep = {{1.0, 3, 5}, {2.0, 4, 7}, {1.0, 1, 0}};
//        double[][] aPrep = {{1.0, 1, 1}, {2.0, 2, 2}, {1.0, 7, 3}}; // TODO: LU doesn't work with this matrix
//        double[][] aPrep = {{2, 5, 7, 6}, {0, -8, -9, 1}, {7, 5, 4, -3} ,{1, 1, 0, 1}};
        Array2DRowRealMatrix a = new Array2DRowRealMatrix(aPrep);

        double[][] bPrep = {{11.0, 9, 24, 2}, {1.0, 5, 2, 6}, {3.0, 17, 18, 1},
                {2.0, 5, 7, 1}};
//        double[][] bPrep = {{1, 1, -1, 2}, {2, 2, 4, 5}, {1, -1, 1, 7},
//                {2.0, 3, 4, 6}};
        Array2DRowRealMatrix b = new Array2DRowRealMatrix(bPrep);

        LU lu = new LU();
        lu.lu_fact(a);
        System.out.println("L: "  + lu.L.toString());
        System.out.println("U: " + lu.U.toString());
        System.out.println("A Error: " + lu.computeError(a));

        lu.L = new Array2DRowRealMatrix(1, 1);
        lu.U = new Array2DRowRealMatrix(1, 1);
        lu.lu_fact(b);
        System.out.println("L: "  + lu.L.toString());
        System.out.println("U: " + lu.U.toString());
        System.out.println("B Error: " + lu.computeError(b));
    }
}
