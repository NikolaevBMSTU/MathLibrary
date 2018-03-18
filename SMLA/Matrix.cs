using System;
using Vectors;
using Matrixs;


namespace Vectors
{
    public class Vector
    {
        protected double[] V;

        public double[] GetVector { get { return V; } }
        public int GetLength { get { return V.GetLength(0); } }

        public Vector()
        {
            V = new double[0];
        }

        public Vector(int Length)
        {
            V = new double[Length];
        }

        public Vector(params double[] VectorComponents)
        {
            V = new double[VectorComponents.GetLength(0)];
            for (int i = 0; i < VectorComponents.GetLength(0); i++)
                V[i] = VectorComponents[i];
        }

        public Vector(Vector X)
        {
            V = new double[X.GetVector.GetLength(0)];
            for (int i = 0; i < X.GetVector.GetLength(0); i++)
                V[i] = X.GetVector[i];
        }

        public static Vector operator +(Vector A, Vector B)
        {
            if (A.GetLength != B.GetLength)
            {
                throw new System.ArgumentException("Невозможно сложить вектора различной длины");
            }
            else
            {
                double[] temp = new double[A.GetLength];
                for (int i = 0; i < A.GetLength; i++)
                    temp[i] = A.GetVector[i] + B.GetVector[i];
                return new Vector(temp);
            }
        }

        public static Vector operator -(Vector A, Vector B)
        {
            if (A.GetLength != B.GetLength)
            {
                throw new System.ArgumentException("Невозможно вычитать вектора различной длины");
            }
            else
            {
                double[] temp = new double[A.GetLength];
                for (int i = 0; i < A.GetLength; i++)
                    temp[i] = A.GetVector[i] - B.GetVector[i];
                return new Vector(temp);
            }
        }

        public static double operator *(Vector A, Vector B)
        {
            if (A.GetLength != B.GetLength)
            {
                throw new System.ArgumentException("Невозможно умножить вектора различной длины");
            }
            else
            {
                double temp = 0;
                for (int i = 0; i < A.GetLength; i++)
                    temp += A.GetVector[i] * B.GetVector[i];
                return temp;
            }
        }

        public static Vector operator *(Vector A, double x)
        {
            double[] temp = new double[A.GetLength];
            for (int i = 0; i < A.GetLength; i++)
                temp[i] = A.GetVector[i] * x;
            return new Vector(temp);
        }

        public static Vector operator *(double x, Vector A)
        {
            double[] temp = new double[A.GetLength];
            for (int i = 0; i < A.GetLength; i++)
                temp[i] = A.GetVector[i] * x;
            return new Vector(temp);
        }

        public static double GetNorm(Vector A)
        {
            return Math.Sqrt(A * A);
        }
    }
}

namespace Matrixs
{
    public class Matrix
    {
        protected int n;
        private int m;
        public int GetN { get { return n; } }
        virtual public int GetM { get { return m; } }

        protected double[,] ArrayElements;

        public double[,] GetArray { get { return ArrayElements; } }

        public Matrix()
        {
            n = 0;
            m = 0;
            ArrayElements = new double[n, m];
        }

        public Matrix(int N, int M)
        {
            if (N < 0 || M < 0) { throw new ArgumentException("Размерности матрицы не могут быть отрицательными"); }
            if (N == 0 || M == 0) { throw new ArgumentException("Размерности матрицы не могут быть нулевыми"); }

            n = N;
            m = M;
            ArrayElements = new double[n, m];
        }

        public Matrix(double[,] array)
        {
            n = array.GetLength(0);
            m = array.GetLength(1);
            ArrayElements = new double[n, m];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < m; j++)
                {
                    ArrayElements[i, j] = array[i, j];
                }
        }

        public Matrix(Vector[] Colums)
        {
            bool NotError = true;
            for (int i = 0; i < Colums.GetLength(0); i++)
            {
                NotError = NotError && (Colums[i].GetLength == Colums[0].GetLength);
            }
            if (!NotError) { throw new ArgumentException("Векторы должны быть одной длины"); }

            n = Colums[0].GetLength;
            m = Colums.GetLength(0);

            ArrayElements = new double[n, m];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < m; j++)
                {
                    ArrayElements[i, j] = Colums[j].GetVector[i];
                }
            }
        }

        public double GetCell(int N, int M)
        {
            if (N < 0 || M < 0) { throw new ArgumentException("Номер ячейки матрицы не может быть отрицательным"); }
            if (N >= n || M >= m) { throw new ArgumentException("Номер ячейки превышает размерность матрицы"); }

            return ArrayElements[N, M];
        }

        public Vector GetColumn(int Numer)
        {
            if (Numer >= m)
            {
                throw new System.ArgumentException("Номер столбца превышает размерность матрицы");
            }
            else
            {
                double[] B = new double[n];
                for (int i = 0; i < n; i++)
                {
                    B[i] = ArrayElements[i, Numer];
                }
                return new Vector(B);
            }
        }

        public Vector GetRow(int Numer)
        {
            if (Numer >= n)
            {
                throw new System.ArgumentException("Номер строки превышает размерность матрицы");
            }
            else
            {
                double[] B = new double[m];
                for (int i = 0; i < m; i++)
                {
                    B[i] = ArrayElements[Numer, i];
                }
                return new Vector(B);
            }
        }

        public static Matrix Transpose(Matrix A)
        {
            Vector[] Rows = new Vector[A.GetN];

            for (int i = 0; i < A.GetN; i++)
            {
                Rows[i] = A.GetRow(i);
            }

            Matrix B = new Matrix(Rows);

            return B;
        }

        public static Matrix operator +(Matrix A, Matrix B)
        {
            if (A.GetN != B.GetM || A.GetM != B.GetN)
            {
                throw new System.ArgumentException("Размерности матриц не совпадают");
            }
            else
            {
                int n = A.GetN;
                int m = B.GetM;

                double[,] Y = new double[n, m];
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                    {
                        Y[i, j] = A.GetCell(i, j) + B.GetCell(i, j);
                    }

                return new Matrix(Y);
            }
        }

        public static Matrix operator -(Matrix A, Matrix B)
        {
            if (A.GetN != B.GetM || A.GetM != B.GetN)
            {
                throw new System.ArgumentException("Размерности матриц не совпадают");
            }
            else
            {
                int n = A.GetN;
                int m = B.GetM;

                double[,] Y = new double[n, m];
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                    {
                        Y[i, j] = A.GetCell(i, j) - B.GetCell(i, j);
                    }

                return new Matrix(Y);
            }
        }

        public static Matrix operator *(Matrix A, Matrix B)
        {
            if (A.GetN != B.GetM || A.GetM != B.GetN)
            {
                throw new System.ArgumentException("Размерности матриц не совпадают");
            }
            else
            {
                int n = A.GetN;
                int m = B.GetM;

                double[,] Y = new double[n, m];
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                    {
                        Y[i, j] = A.GetRow(i) * B.GetColumn(j); ;
                    }

                return new Matrix(Y);
            }
        }
    }

    public class SquareMatrix : Matrix
    {
        public new int GetM { get { return n; } }  // Скрытие неиспользуемого наследуемого члена

        public SquareMatrix()
        {
            n = 0;
            ArrayElements = new double[n, n];
        }

        public SquareMatrix(int N)
        {
            if (N < 0) { throw new ArgumentException("Размерность матрицы не может быть отрицательной"); }
            if (N == 0) { throw new ArgumentException("Размерность матрицы не может быть нулевой"); }

            n = N;
            ArrayElements = new double[n, n];
        }

        public SquareMatrix(double[,] a)
        {
            n = a.GetLength(0);
            if (n != a.GetLength(1)) { throw new ArgumentException("Массив должен быть квадратным"); }
            ArrayElements = new double[n, n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                {
                    ArrayElements[i, j] = a[i, j];
                }
        }

        public SquareMatrix(Vector[] Colums)
        {
            bool NotError = true;
            for (int i = 0; i < Colums.GetLength(0); i++)
            {
                NotError = NotError && (Colums[i].GetLength == Colums[0].GetLength);
            }
            if (!NotError) { throw new ArgumentException("Векторы должны быть одной длины"); }

            if (Colums.GetLength(0) != Colums[0].GetLength) { throw new ArgumentException("Количество векторов не соответствует их длине"); }

            n = Colums[0].GetLength;
            ArrayElements = new double[n, n];
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    ArrayElements[i, j] = Colums[j].GetVector[i];
                }
            }
        }

        public new double GetCell(int N, int M)
        {
            if (N < 0 || M < 0) { throw new ArgumentException("Номер ячейки матрицы не может быть отрицательным"); }
            if (N >= n || M >= n) { throw new ArgumentException("Номер ячейки превышает размерность матрицы"); }

            return ArrayElements[N, M];
        }

        public new Vector GetColumn(int Numer)
        {
            if (Numer >= n)
            {
                throw new System.ArgumentException("Номер столбца превышает размерность матрицы");
            }
            else
            {
                double[] B = new double[n];
                for (int i = 0; i < n; i++)
                {
                    B[i] = ArrayElements[i, Numer];
                }
                return new Vector(B);
            }
        }

        public new Vector GetRow(int Numer)
        {
            if (Numer >= n)
            {
                throw new System.ArgumentException("Номер строки превышает размерность матрицы");
            }
            else
            {
                double[] B = new double[n];
                for (int i = 0; i < n; i++)
                {
                    B[i] = ArrayElements[Numer, i];
                }
                return new Vector(B);
            }
        }

        public static SquareMatrix operator *(SquareMatrix A, SquareMatrix B)
        {
            if (A.GetN != B.GetN)
            {
                throw new System.ArgumentException("Невозможно умножить матрицы различной длины");
            }
            else
            {
                int n = A.GetN;
                double[,] C = new double[n, n];
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        C[i, j] = 0;
                        for (int k = 0; k < n; k++)
                        {
                            C[i, j] += A.GetArray[i, k] * B.GetArray[k, j];
                        }
                    }
                }
                return new SquareMatrix(C);
            }
        }

        public static Vector operator *(SquareMatrix A, Vector X)
        {
            if (A.GetN != X.GetLength)
            {
                throw new System.ArgumentException("Невозможно умножить матрицу и вектор различной длины");
            }
            else
            {
                int n = A.GetN;
                double[] Y = new double[n];
                for (int i = 0; i < n; i++)
                {
                    Y[i] = 0;
                    for (int j = 0; j < n; j++)
                    {
                        Y[i] += A.GetArray[i, j] * X.GetVector[j];
                    }
                }
                return new Vector(Y);
            }   
        }

        public static SquareMatrix operator *(double b, SquareMatrix A)
        {
            int n = A.GetN;
            double[,] Y = new double[n,n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                {
                    Y[i, j] = A.GetArray[i, j] * b;
                }
            
            return new SquareMatrix(Y);
        }

        public static SquareMatrix operator *(SquareMatrix A, double b)
        {
            int n = A.GetN;
            double[,] Y = new double[n, n];
            for (int i = 0; i < n; i++)
                for (int j = 0; j < n; j++)
                {
                    Y[i, j] = A.GetArray[i, j] * b;
                }

            return new SquareMatrix(Y);
        }

        public static SquareMatrix operator +(SquareMatrix A, SquareMatrix B)
        {
            if (A.GetN != B.GetN)
            {
                throw new System.ArgumentException("Невозможно сложить матрицы различной размерности");
            }
            else
            {
                int n = A.GetN;
                double[,] Y = new double[n, n];
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                    {
                        Y[i, j] = A.GetArray[i, j] + B.GetArray[i, j];
                    }

                return new SquareMatrix(Y);
            }
        }

        public static SquareMatrix operator -(SquareMatrix A, SquareMatrix B)
        {
            if (A.GetN != B.GetN)
            {
                throw new System.ArgumentException("Невозможно вычесть матрицы различной длины");
            }
            else
            {
                int n = A.GetN;
                double[,] Y = new double[n, n];
                for (int i = 0; i < n; i++)
                    for (int j = 0; j < n; j++)
                    {
                        Y[i, j] = A.GetArray[i, j] - B.GetArray[i, j];
                    }

                return new SquareMatrix(Y);
            }
        }

        public static SquareMatrix GetUnitMatrix(int n)
        {
            double[,] Y = new double[n, n];
            for (int i = 0; i < n; i++)
                Y[i, i] = 1;
            return new SquareMatrix(Y);
        }

        public SquareMatrix GetInverseMatrix()
        {
            double[,] I = GetUnitMatrix(n).GetArray;
            double[,] B = new double[n, n];

            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    B[i, j] = ArrayElements[i, j];
                }
            }

            for (int k = 0; k < n; k++)
            {
                double Bkk = B[k, k];
                for (int j = n - 1; j >= 0; j--)
                {
                    I[k, j] = I[k, j] / Bkk;
                    B[k, j] = B[k, j] / Bkk;
                }
                for (int i = k + 1; i < n; i++)
                {
                    double Bik = B[i, k];
                    for (int j = 0; j < n; j++)
                    {
                        I[i, j] = I[i, j] - I[k, j] * Bik;
                        B[i, j] = B[i, j] - B[k, j] * Bik;
                    }
                }
            }

            for (int k = n-1; k >= 0; k--)
            {
                for (int i = 0; i < k; i++)
                {
                    double Bik = B[i, k];
                    for (int j = 0; j < n; j++)
                    {
                        I[i, j] = I[i, j] - I[k, j] * Bik;
                        B[i, j] = B[i, j] - B[k, j] * Bik;
                    }
                }
            }

            return new SquareMatrix(I);
        }

        public void ILU(out SquareMatrix L, out SquareMatrix U, out SquareMatrix R)
        {
            double[,] l = new double[ArrayElements.GetLength(0), ArrayElements.GetLength(1)];
            double[,] u = new double[ArrayElements.GetLength(0), ArrayElements.GetLength(1)];

            for (int i = 0; i < ArrayElements.GetLength(0); i++)
                for (int j = 0; j < ArrayElements.GetLength(1); j++)
                {
                    l[i, j] = 0/*A[i, j]*/;
                    u[i, j] = 0/*A[i, j]*/;
                }

            //R = new SquareMatrix(A.GetArray);

            for (int k = 0; k < n; k++)
            {
                for (int j = 0; j < k - 1; j++)
                {
                    if (ArrayElements[k, j] == 0) { l[k, j] = 0; }
                    else
                    {
                        double summ = 0;
                        for (int i = 0; i < j - 1; i++)
                            summ += l[k, i] * u[i, j];
                        l[k, j] = (ArrayElements[k, j] - summ) / u[j, j];
                    }
                }
                l[k, k] = 1;
                for (int j = k; j < n; j++)
                {
                    if (ArrayElements[k, j] == 0) { u[k, j] = 0; }
                    else
                    {
                        double summ = 0;
                        for (int i = 0; i < k - 1; i++)
                            summ += l[k, i] * u[i, j];
                        u[k, j] = ArrayElements[k, j] - summ;
                    }
                }
            }

            L = new SquareMatrix(l);
            U = new SquareMatrix(u);
            R = new SquareMatrix(ArrayElements) - L * U;
        }

        public static double GetMaxEigenvalue(SquareMatrix A, double eps)
        {
            double[] x0 = new double[A.GetN];
            for (int i = 0; i < x0.GetLength(0); i++) x0[i] = 1;

            double Lambda = new double();
            double NextLambda = new double();
            Vector X = new Vector(x0);
            Vector NextX = new Vector(x0);

            bool end = false;
            while (!end)
            {
                X = NextX;
                Lambda = NextLambda;
                NextX = A * X;
                NextLambda = NextX.GetVector[0] / X.GetVector[0];

                if (Math.Abs(NextLambda - Lambda) < eps) { end = true; }
            }

            return NextLambda;
        }

    }

}

namespace SMLA
// Some methods of linear algebra
// Некоторые методы линейной алгебры

{
    public class SLEQ
        // System of linear equations
        // Методы решения систем линейных уравнений
    {
        public static Vector MRES(SquareMatrix A, Vector F, Vector X0)
        {
            Vector X = new Vector();
            Vector NextX = new Vector();
            Vector R = new Vector();
            Vector AR = new Vector();
            double tau = 0;

            NextX = new Vector(X0);
            bool end = false;
            bool check = false;
            while (!end)
            {
                X = new Vector(NextX);
                R = A * X - F;
                AR = A * R;
                tau = (AR * R) / Math.Pow(Vector.GetNorm(AR), 2);

                if (!check)
                {
                    SquareMatrix S = SquareMatrix.GetUnitMatrix(A.GetN) - tau * A;
                    double Lambda1 = SquareMatrix.GetMaxEigenvalue(S, 0.01);
                    if (Lambda1 >= 1) throw new System.ArgumentException("Процесс решения не сходится.");
                    check = true;
                }

                NextX = X - (tau * R);
                if (Vector.GetNorm(R) < 0.01) { end = true; }

            }

            return X;
        }

        public static Vector BiCGStab(SquareMatrix A, Vector B, Vector X0, int MaxIter, double eps)
        {
            Vector[] X = new Vector[MaxIter];
            Vector[] V = new Vector[MaxIter];
            Vector[] R = new Vector[MaxIter];
            Vector[] P = new Vector[MaxIter];
            Vector[] S = new Vector[MaxIter];
            Vector[] T = new Vector[MaxIter];

            Vector R0 = B - A * X0;

            double[] alpha = new double[MaxIter];
            double[] beta = new double[MaxIter];
            double[] rho = new double[MaxIter];
            double[] omega = new double[MaxIter];

            double[] NormR = new double[MaxIter];

            R[0] = R0;

            alpha[0] = 1;
            rho[0] = 1;
            omega[0] = 1;

            V[0] = 0 * R0;
            P[0] = 0 * R0;
            X[0] = new Vector(X0);

            for (int k = 1; k < MaxIter; k++)
            {
                rho[k] = R0 * R[k - 1];
                beta[k] = rho[k] / rho[k - 1] * alpha[k - 1] / omega[k - 1];
                P[k] = R[k - 1] + beta[k] * (P[k - 1] - omega[k - 1] * V[k - 1]);
                V[k] = A * P[k];
                alpha[k] = rho[k] / (R0 * V[k]);
                S[k] = R[k - 1] - alpha[k] * V[k];
                T[k] = A * S[k];
                omega[k] = (T[k] * S[k]) / (T[k] * T[k]);
                X[k] = X[k - 1] + omega[k] * S[k] + alpha[k] * P[k];
                R[k] = S[k] - omega[k] * T[k];
                NormR[k] = Vector.GetNorm(R[k]);

                if (Vector.GetNorm(R[k]) < eps) { return X[k]; }
            }

            return X[MaxIter - 1];
        }

        public static Vector BiCGStab(SquareMatrix A, Vector B, Vector X0)
        {
            return BiCGStab(A, B, X0, 100000, 1e-4);
        }

        public static SquareMatrix InverseMatrixBiCGStab(SquareMatrix A)
        {
            SquareMatrix E = SquareMatrix.GetUnitMatrix(A.GetN);
            Vector[] B = new Vector[A.GetN];

            for (int i = 0; i < A.GetN; i++) 
            {
                B[i] = BiCGStab(A, E.GetColumn(i), E.GetColumn(i));
            }

            SquareMatrix IA = new SquareMatrix(B);
            return IA;
        }

        public static SquareMatrix InverseMatrixMRES(SquareMatrix A)
        {
            SquareMatrix E = SquareMatrix.GetUnitMatrix(A.GetN);
            Vector[] B = new Vector[A.GetN];

            for (int i = 0; i < A.GetN; i++)
            {
                B[i] = MRES(A, E.GetColumn(i), E.GetColumn(i));
            }

            SquareMatrix IA = new SquareMatrix(B);
            return IA;
        }
    }
}