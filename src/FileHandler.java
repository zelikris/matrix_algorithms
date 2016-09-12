import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.ArrayRealVector;

import java.io.File;
import java.io.FileNotFoundException;
import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.Scanner;

/**
 * Created by Robert on 4/15/2016.
 */
public class FileHandler
{
    public static Array2DRowRealMatrix readMatrixFromFile(String fileName)
    {
        Array2DRowRealMatrix matrix = new Array2DRowRealMatrix();

        File file = new File(fileName);
        try
        {
            List<Double> matrixFileRead = new LinkedList<>();
            int columnCount = 0;
            int rowCount = 0;
            Scanner fileReader = new Scanner(file);
            while (fileReader.hasNextLine())
            {
                Scanner lineReader = new Scanner(fileReader.nextLine());
                while (lineReader.hasNextDouble())
                {
                    matrixFileRead.add(lineReader.nextDouble());
                }
                if (columnCount == 0)
                {
                    columnCount = matrixFileRead.size();
                }
                rowCount++;
            }
            matrix = (Array2DRowRealMatrix) matrix.createMatrix(rowCount, columnCount); //assume square matrix

            for (int i = 0; i < rowCount; i++)
            {
                for (int j = 0; j < columnCount; j++)
                {
                    //System.out.printf("(%d,%d)=%d\n", i, j, i * columnCount + j);
                    matrix.setEntry(i, j, matrixFileRead.get(i * columnCount + j));
                }
            }
            fileReader.close();
        }
        catch (FileNotFoundException e)
        {
            System.out.println("File Not Found: " + fileName);
            System.out.println(e);
        }

        return matrix;
    }

    public static void main(String[] args)
    {
        if(args.length >=  2)
        {
            System.out.println("Reading matrix from file: \"" + args[0] + "\"");
        }
        else
        {
            System.out.println("Argument error, check readme.txt for proper syntax.");
            System.exit(0);
        }

        Array2DRowRealMatrix matrix = readMatrixFromFile(args[0]);
        System.out.println("Matrix:");
        printMatrix(matrix);

        System.out.println("\nFunction: " + args[1].substring(1));
        if (args[1].equals("-lu_fact"))
        {
            //LU Matrix
            LU decomposedMatrix = new LU();
            decomposedMatrix.lu_fact(matrix);

            System.out.println("\nL:");
            printMatrix(decomposedMatrix.getL());

            System.out.println("\nU:");
            printMatrix(decomposedMatrix.getU());

            System.out.printf("\nError ||LU - A||: %.8f", decomposedMatrix.computeError(matrix));
        }

        if (args[1].equals("-qr_fact_house"))
        {
            //QR Householder matrix
            QR houseMatrix = new QR();
            houseMatrix.qr_fact_givens(matrix);

            System.out.println("\nQ (Householder):");
            printMatrix(houseMatrix.getQ());

            System.out.println("\nR (Householder):");
            printMatrix(houseMatrix.getR());

            System.out.printf("\nError ||QR - A|| (Householder): %.8f", houseMatrix.computeError(matrix));
        }

        if (args[1].equals("-qr_fact_givens"))
        {
            //QR Givens matrix
            QR givensMatrix = new QR();
            givensMatrix.qr_fact_givens(matrix);

            System.out.println("\nQ (Givens):");
            printMatrix(givensMatrix.getQ());

            System.out.println("\nR (Givens):");
            printMatrix(givensMatrix.getR());

            System.out.printf("\nError ||QR - A|| (Givens): %.8f", givensMatrix.computeError(matrix));
        }

        if(args[1].equals("-solve_lu"))
        {
            Array2DRowRealMatrix matrixA = (Array2DRowRealMatrix) matrix.getSubMatrix(0, matrix.getRowDimension()-1, 0, matrix.getColumnDimension()-2);
            Array2DRowRealMatrix matrixB = (Array2DRowRealMatrix) matrix.getSubMatrix(0, matrix.getRowDimension()-1, matrix.getColumnDimension()-1, matrix.getColumnDimension()-1);
            LU decomposedMatrix = new LU();
            decomposedMatrix.lu_fact(matrixA);
            Array2DRowRealMatrix matrixX = decomposedMatrix.solve_lu(matrixA, matrixB);

            System.out.println("Solution x: ");
            printMatrix(matrixX);

            System.out.printf("\nError ||Ax - b|| (LU): %.8f", Part1.computeError(matrixA, matrixX, matrixB));
        }

        if(args[1].equals("-solve_qr_house"))
        {
            Array2DRowRealMatrix matrixA = (Array2DRowRealMatrix) matrix.getSubMatrix(0, matrix.getRowDimension()-1, 0, matrix.getColumnDimension()-2);
            Array2DRowRealMatrix matrixB = (Array2DRowRealMatrix) matrix.getSubMatrix(0, matrix.getRowDimension()-1, matrix.getColumnDimension()-1, matrix.getColumnDimension()-1);
            QR decomposedMatrix = new QR();
            decomposedMatrix.qr_fact_house(matrixA);
            Array2DRowRealMatrix matrixX = decomposedMatrix.solve_qr_house(matrixA, matrixB);

            System.out.println("Solution x: ");
            printMatrix(matrixX);

            System.out.printf("\nError ||Ax - b|| (Householder): %.8f", Part1.computeError(matrixA, matrixX, matrixB));
        }

        if(args[1].equals("-solve_qr_givens"))
        {
            Array2DRowRealMatrix matrixA = (Array2DRowRealMatrix) matrix.getSubMatrix(0, matrix.getRowDimension()-1, 0, matrix.getColumnDimension()-2);
            Array2DRowRealMatrix matrixB = (Array2DRowRealMatrix) matrix.getSubMatrix(0, matrix.getRowDimension()-1, matrix.getColumnDimension()-1, matrix.getColumnDimension()-1);
            QR decomposedMatrix = new QR();
            decomposedMatrix.qr_fact_givens(matrixA);
            Array2DRowRealMatrix matrixX = decomposedMatrix.solve_qr_givens(matrixA, matrixB);

            System.out.println("Solution x: ");
            printMatrix(matrixX);

            System.out.printf("\nError ||Ax - b|| (Givens): %.8f", Part1.computeError(matrixA, matrixX, matrixB));
        }

        if(args[1].equals("-jacobi_iter"))
        {
            Iterations iterativeMatrix = new Iterations();
            Array2DRowRealMatrix startingVector;
            double tolerance = Double.parseDouble(args[3].substring(1));
            int iterationLimit = Integer.parseInt(args[4].substring(1));

            Scanner arrayScanner = new Scanner(args[2].substring(1));
            arrayScanner.useDelimiter(", *");

            List<Double> vectorPoints = new ArrayList<>();
            while (arrayScanner.hasNextDouble())
            {

                vectorPoints.add(arrayScanner.nextDouble());
            }
            double[] startVect = new double[vectorPoints.size()];
            for(int i = 0; i < vectorPoints.size(); i++)
            {
                startVect[i] = vectorPoints.get(i);
            }
            startingVector = new Array2DRowRealMatrix(startVect);

            System.out.printf("Tolerance: %.12f\n", tolerance);
            System.out.println("Max Iterations: " + iterationLimit);
            System.out.println("Starting vector:");
            printMatrix(startingVector);

            Array2DRowRealMatrix result = iterativeMatrix.jacobi_iter(matrix, tolerance, iterationLimit, startingVector);

            System.out.println("\nJacobi iterations used: " + iterativeMatrix.getIterations());

            System.out.println("Output solution x:");
            printMatrix(result);

            System.out.printf("Error ||Ax-b||: %.12f\n", iterativeMatrix.computeError(matrix, result));

            if (!iterativeMatrix.didSucceed())
            {
                System.out.println("***********************************************************");
                System.out.println("FAILURE: Did not converge within specified # of iterations.");
            }
        }

        if(args[1].equals("-gs_iter"))
        {
            Iterations iterativeMatrix = new Iterations();
            Array2DRowRealMatrix startingVector;
            double tolerance = Double.parseDouble(args[3].substring(1));
            int iterationLimit = Integer.parseInt(args[4].substring(1));

            Scanner arrayScanner = new Scanner(args[2].substring(1));
            arrayScanner.useDelimiter(", *");

            List<Double> vectorPoints = new ArrayList<>();
            while (arrayScanner.hasNextDouble())
            {

                vectorPoints.add(arrayScanner.nextDouble());
            }
            double[] startVect = new double[vectorPoints.size()];
            for(int i = 0; i < vectorPoints.size(); i++)
            {
                startVect[i] = vectorPoints.get(i);
            }
            startingVector = new Array2DRowRealMatrix(startVect);

            System.out.printf("Tolerance: %.12f\n", tolerance);
            System.out.println("Max Iterations: " + iterationLimit);
            System.out.println("Starting vector:");
            printMatrix(startingVector);

            Array2DRowRealMatrix result = iterativeMatrix.gs_iter(matrix, tolerance, iterationLimit, startingVector);

            System.out.println("\nGauss-Seidel iterations used: " + iterativeMatrix.getIterations());

            System.out.println("Output solution x:");
            printMatrix(result);

            System.out.printf("Error ||Ax-b||: %.12f\n", iterativeMatrix.computeError(matrix, result));

            if (!iterativeMatrix.didSucceed())
            {
                System.out.println("***********************************************************");
                System.out.println("FAILURE: Did not converge within specified # of iterations.");
            }
        }

        if(args[1].equals("-power_method"))
        {
            Iterations iterativeMatrix = new Iterations();
            Array2DRowRealMatrix startingVector, auxVector;
            double tolerance = Double.parseDouble(args[4].substring(1));
            int iterationLimit = Integer.parseInt(args[5].substring(1));

            Scanner arrayScanner = new Scanner(args[2].substring(1));
            arrayScanner.useDelimiter(", *");

            List<Double> vectorPoints = new ArrayList<>();
            while (arrayScanner.hasNextDouble())
            {

                vectorPoints.add(arrayScanner.nextDouble());
            }
            double[] startVect = new double[vectorPoints.size()];
            for(int i = 0; i < vectorPoints.size(); i++)
            {
                startVect[i] = vectorPoints.get(i);
            }
            startingVector = new Array2DRowRealMatrix(startVect);

            arrayScanner = new Scanner(args[3].substring(1));
            arrayScanner.useDelimiter(", *");
            vectorPoints = new ArrayList<>();
            while (arrayScanner.hasNextDouble())
            {

                vectorPoints.add(arrayScanner.nextDouble());
            }
            double[] auxVect = new double[vectorPoints.size()];
            for(int i = 0; i < vectorPoints.size(); i++)
            {
                auxVect[i] = vectorPoints.get(i);
            }
            auxVector = new Array2DRowRealMatrix(auxVect);

            System.out.printf("Tolerance: %.12f\n", tolerance);
            System.out.println("Max Iterations: " + iterationLimit);
            System.out.println("Starting vector:");
            printMatrix(startingVector);
            System.out.println("\nAuxillary vector:");
            printMatrix(auxVector);

            double result = iterativeMatrix.power_method(matrix,
                                                         (ArrayRealVector) startingVector.getColumnVector(0),
                                                         (ArrayRealVector) auxVector.getColumnVector(0),
                                                         tolerance,
                                                         iterationLimit);
            System.out.println("\nPower Method iterations used: " + iterativeMatrix.getIterations());
            System.out.println("Max Eigenvalue: " + result);
            System.out.println("Corresponding Eigenvector:");
            printMatrix(iterativeMatrix.getEigenvector());
            if (!iterativeMatrix.didSucceed())
            {
                System.out.println("***********************************************************");
                System.out.println("FAILURE: Did not converge within specified # of iterations.");
            }
        }
    }

    public static void printMatrix(Array2DRowRealMatrix matrix)
    {
        double[][] matrixBackingArray = matrix.getData();

        for (int i = 0; i < matrixBackingArray.length; i++) {
            for (int j = 0; j < matrixBackingArray[0].length; j++) {
                System.out.printf("%.8f\t", matrixBackingArray[i][j]);
            }
            System.out.print("\n");
        }
    }
}
