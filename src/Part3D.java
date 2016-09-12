import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.linear.MatrixUtils;

/**
 * Created by Kris on 4/19/2016.
 */
public class Part3D {

    public static void main(String[] args) {
        double[][] aPrep = {{-2, 1, 2}, {0, 2, 3}, {2, 1, -2}};
        double[] u0Prep = {1, 0, 0};
        Array2DRowRealMatrix A = new Array2DRowRealMatrix(aPrep);
        Array2DRowRealMatrix I = (Array2DRowRealMatrix) MatrixUtils.createRealIdentityMatrix(3);
        ArrayRealVector u0 = new ArrayRealVector(u0Prep);
        ArrayRealVector w0 = u0.copy();
        double p1 = -(5.0 / 2.0);
        double p2 = (5.0 / 2.0);
        A = (Array2DRowRealMatrix) A.subtract(I.scalarMultiply(p1));
        A = Part3.inverse(A);
        Iterations it = new Iterations();
        System.out.println(it.power_method(A, u0, w0, 0.00005, 100));
    }
}
