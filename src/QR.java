import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;

/**
 * Created by Kris on 4/12/2016.
 */
public class QR {
    private Array2DRowRealMatrix Q;
    private Array2DRowRealMatrix R;

    public void qr_fact_house(Array2DRowRealMatrix A) {
        A = (Array2DRowRealMatrix) A.copy();
        int n = A.getColumnDimension();
        double sum, A0, magnitude;
        byte sign;
        Q = (Array2DRowRealMatrix)MatrixUtils.createRealIdentityMatrix(n);
        R = (Array2DRowRealMatrix)MatrixUtils.createRealIdentityMatrix(n);
        ArrayRealVector v = new ArrayRealVector(n);
        for(int k = 0; k < n; k++) {
            sum = 0;
            sign = 1;
            A0 = A.getEntry(k, k);
            if(A0 < 0){
                sign = -1;
            }
            for(int i = k; i < n; i++) {
                sum += A.getEntry(i, k) * A.getEntry(i, k);
            }
            magnitude = sign * Math.sqrt(sum);
            v.setEntry(k, magnitude + A0);
            magnitude = Math.sqrt(2 * (sum + A0*magnitude));
            v.setEntry(k, v.getEntry(k) / magnitude);
            for(int i = k + 1; i < n; i++) {
                v.setEntry(i, A.getEntry(i, k) / magnitude);
            }
            for(int j = 0; j < n; j++) {
                sum = 0;
                for(int i = k; i < n; i++) {
                    sum += v.getEntry(i) * A.getEntry(i, j);
                }
                for(int i = k; i < n; i++) {
                    A.setEntry(i , j, A.getEntry(i, j) - 2 * v.getEntry(i) * sum);
                }
            }
            for(int j = 0; j < n; j++) {
                sum = 0;
                for(int i = k; i < n; i++) {
                    sum += v.getEntry(i) * Q.getEntry(i, j);
                }
                for(int i = k; i < n; i++) {
                    Q.setEntry(i, j, Q.getEntry(i, j) - 2 * v.getEntry(i) * sum);
                }
            }
        }
        Q = transpose(Q);
        R = A;
    }

    public void qr_fact_givens(Array2DRowRealMatrix A) {
        int n = A.getColumnDimension();
        Array2DRowRealMatrix G = (Array2DRowRealMatrix)MatrixUtils.createRealIdentityMatrix(n);
        Q = (Array2DRowRealMatrix)MatrixUtils.createRealIdentityMatrix(n);
        R = (Array2DRowRealMatrix)MatrixUtils.createRealIdentityMatrix(n);
        double a, b, cos, sin;

        for (int i = 0; i < n; i++) {
            for (int j = (n - 1); j > i; j--) {
                a = A.getEntry(j - 1, i);
                b = A.getEntry(j, i);
                cos = a / (Math.sqrt(a * a + b * b));
                sin = -b / (Math.sqrt(a * a + b * b));

                G.setEntry(j, j, cos);
                G.setEntry(j, j - 1, sin);
                G.setEntry(j - 1, j, -sin);
                G.setEntry(j - 1, j - 1, cos);

                A = G.multiply(A);
                Q = G.multiply(Q);

                G = (Array2DRowRealMatrix)MatrixUtils.createRealIdentityMatrix(n);
            }
        }
        Q = transpose(Q);
        R = A;
    }

    public Array2DRowRealMatrix solve_qr_house(Array2DRowRealMatrix Ab) {
        int rows = Ab.getRowDimension();
        int cols = Ab.getColumnDimension();
        Array2DRowRealMatrix A = (Array2DRowRealMatrix)Ab.getSubMatrix(0, rows - 1, 0, cols - 2);
        Array2DRowRealMatrix b = (Array2DRowRealMatrix)Ab.getSubMatrix(0, rows - 1, cols - 1, cols - 1);

        qr_fact_house(A);
        Array2DRowRealMatrix y = transpose(Q).multiply(b);
        return backwardSub(y);
    }

    public Array2DRowRealMatrix solve_qr_house(Array2DRowRealMatrix A, Array2DRowRealMatrix b) {
        qr_fact_house(A);
        Array2DRowRealMatrix y = transpose(Q).multiply(b);
        return backwardSub(y);
    }

    public Array2DRowRealMatrix solve_qr_givens(Array2DRowRealMatrix A, Array2DRowRealMatrix b) {
        qr_fact_givens(A);
        Array2DRowRealMatrix y = transpose(Q).multiply(b);
        return backwardSub(y);
    }

    private Array2DRowRealMatrix backwardSub(Array2DRowRealMatrix y) {
        int n = y.getRowDimension();
        Array2DRowRealMatrix x = new Array2DRowRealMatrix(n, 1);

        // backward substitution for Ux = y
        for (int i = n - 1; i >= 0; i--) {
            x.setEntry(i, 0, y.getEntry(i, 0));
            for (int j = i + 1; j < n; j++) {
                x.setEntry(i, 0, x.getEntry(i , 0) - (R.getEntry(i, j) * x.getEntry(j, 0)));
            }
            x.setEntry(i, 0, x.getEntry(i, 0) / R.getEntry(i, i));
        }
        return x;
    }

    public static Array2DRowRealMatrix transpose(Array2DRowRealMatrix A) {
        int n = A.getRowDimension();
        Array2DRowRealMatrix result = new Array2DRowRealMatrix(n, n);
        for (int j = 0; j < n; j++)
            for (int i = 0; i < n; i++)
                result.setEntry(i, j, A.getEntry(j, i));
        return result;
    }

    public double computeError(Array2DRowRealMatrix A) {
        return Q.multiply(R).subtract(A).getNorm(); // TODO: norm()
    }

    public Array2DRowRealMatrix getR() {
        return R;
    }

    public void setR(Array2DRowRealMatrix r) {
        R = r;
    }

    public Array2DRowRealMatrix getQ() {
        return Q;
    }

    public void setQ(Array2DRowRealMatrix q) {
        Q = q;
    }

    public static void main(String[] args) {
//        System.out.println(generateHilbert(2));
//        System.out.println(generateB(2));
//        double[][] aPrep = {{-2,1,2}, {0,2,3}, {2,1,-2}};
//        Array2DRowRealMatrix A = new Array2DRowRealMatrix(aPrep);
//        QR qr = new QR();
//        qr.qr_fact_house(A);
//        System.out.println("Q: " + qr.getQ());
//        System.out.println("R: " + qr.getR());
//        System.out.println("HError: " + qr.computeError(A));
//
//        A = new Array2DRowRealMatrix(aPrep);
//        qr. qr_fact_givens(A);
//        System.out.println("QGivens: " + qr.getQ());
//        System.out.println("RGivens: " + qr.getR());
//        System.out.println("GError: " + qr.computeError(A));

//        double[][] aPrep = {{4,1,1}, {3,2,5}, {7,7,1}};
//        double[][] AbPrep = {{1.0, 3, 5}, {2.0, 4, 7}};
//        double[][] AbPrep1 = {{2, 1, 1, 8}, {4, -6, 0, 2}, {2, 7, 2, 3}};
//
//        Array2DRowRealMatrix A = new Array2DRowRealMatrix(aPrep);
//        Array2DRowRealMatrix Ab = new Array2DRowRealMatrix(AbPrep);
//        Array2DRowRealMatrix Ab1 = new Array2DRowRealMatrix(AbPrep1);
//        QR qr = new QR();
//
//        System.out.println("X: " + qr.solve_qr_house(Ab));
//        System.out.println("X: " + qr.solve_qr_house(Ab1));
//
//        qr.qr_fact_house(A);
//        A = new Array2DRowRealMatrix(aPrep); // not sure why this is necessary, but it is
//
//        System.out.println("\nHouse: ");
//        System.out.println("Q After: " + qr.Q.toString());
//        System.out.println("R After: " + qr.R.toString());
//        System.out.println("QR: " + qr.Q.multiply(qr.R));
//        System.out.println("Error: " + qr.computeError(A) + "\n");
//
//        qr.qr_fact_givens(A);
//        System.out.println("Givens: ");
//        System.out.println("Q After: " + qr.Q.toString());
//        System.out.println("R After: " + qr.R.toString());
//        System.out.println("QR: " + qr.Q.multiply(qr.R));
//        System.out.println("Error: " + qr.computeError(A));
    }
}
