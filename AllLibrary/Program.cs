using System;
using SMAN;
using SMLA;
using Matrixs;

namespace AllLibrary
{
    public class Program
    {
        public static double fun(double x)
        {
            return Math.Exp(x);
        }

        public static double funXY(double x, double y)
        {
            return x + y;
        }

        public static void Main(string[] args)
        {
            string a = Convert.ToString(Integrals.Integrate(0, 1, fun));
            string b = Convert.ToString(Integrals.Integrate(0, 1, 1e-12, fun));
            Console.WriteLine(String.Concat("Integral f(x)=e^x from 0 to 1 = ", a));
            Console.WriteLine(String.Concat("Integral f(x)=e^x from 0 to 1 = ", b));
            double[] c = ODE.Euler(0, 0, 1, 10, funXY);
            double[,] A = new double[3, 3] { { 8, 2, 4 }, { 4, 6, 8 }, { 7, 5, 10 } };
            SquareMatrix MatrixA = new SquareMatrix(A);
            int tst = MatrixA.GetM;
            SquareMatrix Actual = MatrixA * SLEQ.InverseMatrixBiCGStab(MatrixA);
        }
    }
}
