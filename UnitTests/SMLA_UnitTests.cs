using System;
using Scalar;
using Vectors;
using Matrixs;
using SMLA;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace UnitTests
{
    [TestClass]
    public class SCAClassTests
    {
        [TestMethod]
        public void SummScalarTest()
        {
            SCA A = new SCA(12, 0.0001);
            SCA B = new SCA(1, 0.1);
            SCA C = A + B;

            double ExpectedValue = 13;
            double ActualValue = C.Value;

            double ExpectedAccuracy = 0.1001;
            double ActualAccuracy = C.Delta;

            Assert.AreEqual(ExpectedValue, ActualValue, 0.00001, "Значение выражения вычисляется неверно");
            Assert.AreEqual(ExpectedAccuracy, ActualAccuracy, 0.00001, "Погрешность выражения вычисляется неверно");
        }
                
        [TestMethod]
        public void ReductionOfTypeScalarTest()
        {
            SCA A = new SCA(12, 0.0001);
            double B = 10;
            double C = A + B;

            double ExpectedValue = 22;
            double ActualValue = C;

            Assert.AreEqual(ExpectedValue, ActualValue, 0.00001, "Значение выражения вычисляется неверно");
        }
    }

    [TestClass]
    public class VectorsClassTests
    {

        [TestMethod]
        public void SummVectorTest()
        {
            Vector A = new Vector(12, 24);
            Vector B = new Vector(28, -4);
            Vector C = A + B;

            double expected = 44.721;
            double actual = Vector.GetNorm(C);
            Assert.AreEqual(expected, actual, 0.001, "Norm not calculated correctly");
        }

        [TestMethod]
        public void SummDifferentVectorsTest()
        {
            Vector A = new Vector(12, 24);
            Vector B = new Vector(28, -4, 0);
            bool actual = false;
            try
            {
                Vector C = A + B;
            }
            catch (System.Exception e)
            {
                if (e.Message == "Невозможно сложить вектора различной длины") actual = true;
            }

            Assert.IsTrue(actual);
        }

        [TestMethod]
        public void SubtractionVectorTest()
        {
            Vector A = new Vector(12, 24, 10);
            Vector B = new Vector(28, -4, 0);
            Vector C = A - B;

            double expected = 33.764;
            double actual = Vector.GetNorm(C);
            Assert.AreEqual(expected, actual, 0.001, "Norm not calculated correctly");
        }

        [TestMethod]
        public void SubtractionDifferentVectorsTest()
        {
            Vector A = new Vector(12, 24, 0, 85);
            Vector B = new Vector(28, -4, 0);
            bool actual = false;
            try
            {
                Vector C = A - B;
            }
            catch (System.Exception e)
            {
                if (e.Message == "Невозможно вычитать вектора различной длины") actual = true;
            }

            Assert.IsTrue(actual);
        }

        [TestMethod]
        public void ScalarMultiplicationVectorsTest()
        {
            Vector A = new Vector(12, 24, 0, 85, 33, 87, 1);
            Vector B = new Vector(28, -4, 0, -23, 8, 17, 65);
            
            double expected = 93;
            double actual = A * B;
            Assert.AreEqual(expected, actual, 0.01, "Scalar Multiplication not calculated correctly");
        }

        [TestMethod]
        public void ScalarMultiplicationDifferentVectorsTest()
        {
            Vector A = new Vector(12, 24, 0, 85, 33, 87, 1);
            Vector B = new Vector(28, -4, 0, -23, 8, 17);

            bool actual = false;
            try
            {
                double ActualValue = A * B;
            }
            catch (System.Exception e)
            {
                if (e.Message == "Невозможно умножить вектора различной длины") actual = true;
            }

            Assert.IsTrue(actual);
        }
    }

    [TestClass]
    public class MatrixClassTests
    {
        [TestMethod]
        public void CreateNegativDimensionMatrixTest()
        {
            bool actual = false;
            try
            {
                Matrix ActualC = new Matrix(-1, -2);
            }
            catch (System.Exception e)
            {
                if (e.Message == "Размерности матрицы не могут быть отрицательными") actual = true;
            }

            Assert.IsTrue(actual, "Создана матрица с отрицательными размерностями");
        }

        [TestMethod]
        public void CreateNegativN_MatrixTest()
        {
            bool actual = false;
            try
            {
                Matrix ActualC = new Matrix(-1, 3);
            }
            catch (System.Exception e)
            {
                if (e.Message == "Размерности матрицы не могут быть отрицательными") actual = true;
            }

            Assert.IsTrue(actual, "Создана матрица с отрицательными n");
        }

        [TestMethod]
        public void CreateNegativM_MatrixTest()
        {
            bool actual = false;
            try
            {
                Matrix ActualC = new Matrix(1, -3);
            }
            catch (System.Exception e)
            {
                if (e.Message == "Размерности матрицы не могут быть отрицательными") actual = true;
            }

            Assert.IsTrue(actual, "Создана матрица с отрицательными m");
        }

        [TestMethod]
        public void MultiplicationMatrix3x3Test()
        {
            double[,] a = new double[3, 3] { { 8, 2, 4 }, { 4, 6, 8 }, { 7, 5, 10 } };
            double[,] b = new double[3, 3] { { 8, 2, 4 }, { 4, 6, 8 }, { 7, 5, 10 } };
            Matrix A = new Matrix(a);
            Matrix B = new Matrix(b);

            double[,] expected = new double[3, 3] { { 100, 48, 88 }, { 112, 84, 144 }, { 146, 94, 168 } };
            Matrix ExpectedC = new Matrix(expected);

            Matrix ActualC = A * B;

            for (int i = 0; i < A.GetN; i++)
                for (int j = 0; j < A.GetN; j++)
                {
                    Assert.AreEqual(ExpectedC.GetArray[i, j], ActualC.GetArray[i, j], 0.1, "Multiplication not calculated correctly");
                }
        }

        [TestMethod]
        public void MultiplicationMatrix2x3and3x2Test()
        {
            double[,] a = new double[2, 3] { { 3, 1, 2 }, { 2, 1, 3 }};
            double[,] b = new double[3, 2] { { 1, 9 }, { 0, 11 }, { 8, 2 } };
            Matrix A = new Matrix(a);
            Matrix B = new Matrix(b);

            double[,] expected = new double[2, 2] { { 19, 42 }, { 26, 35 } };
            Matrix ExpectedC = new Matrix(expected);

            Matrix ActualC = A * B;

            for (int i = 0; i < A.GetN; i++)
                for (int j = 0; j < A.GetN; j++)
                {
                    Assert.AreEqual(ExpectedC.GetArray[i, j], ActualC.GetArray[i, j], 0.1, "Multiplication not calculated correctly");
                }
        }

        [TestMethod]
        public void MultiplicationMatrix3x2and2x3Test()
        {
            double[,] a = new double[2, 3] { { 3, 1, 2 }, { 2, 1, 3 } };
            double[,] b = new double[3, 2] { { 1, 9 }, { 0, 11 }, { 8, 2 } };
            Matrix A = new Matrix(a);
            Matrix B = new Matrix(b);

            double[,] expected = new double[3, 3] { { 21, 10, 29 }, { 22, 11, 33 }, { 28, 10, 22 } };
            Matrix ExpectedC = new Matrix(expected);

            Matrix ActualC = B * A;

            for (int i = 0; i < A.GetN; i++)
                for (int j = 0; j < A.GetN; j++)
                {
                    Assert.AreEqual(ExpectedC.GetArray[i, j], ActualC.GetArray[i, j], 0.1, "Multiplication not calculated correctly");
                }
        }

        [TestMethod]
        public void AdditionMatrixTest()
        {
            double[,] a = new double[3, 3] { { 8, 2, 4 }, { 4, 6, 8 }, { 7, 5, 10 } };
            double[,] b = new double[3, 3] { { 8, 2, 4 }, { 4, 6, 8 }, { 7, 5, 10 } };
            Matrix A = new Matrix(a);
            Matrix B = new Matrix(b);

            double[,] expected = new double[3, 3] { { 16, 4, 8 }, { 8, 12, 16 }, { 14, 10, 20 } };
            Matrix ExpectedC = new Matrix(expected);

            Matrix ActualC = A + B;

            for (int i = 0; i < A.GetN; i++)
                for (int j = 0; j < A.GetN; j++)
                {
                    Assert.AreEqual(ExpectedC.GetArray[i, j], ActualC.GetArray[i, j], 0.1, "Addition not calculated correctly");
                }
        }

        [TestMethod]
        public void AdditionDifferentMatrixTest()
        {
            double[,] a = new double[2, 2] { { 8, 2 }, { 4, 6 } };
            double[,] b = new double[3, 3] { { 8, 2, 4 }, { 4, 6, 8 }, { 7, 5, 10 } };
            Matrix A = new Matrix(a);
            Matrix B = new Matrix(b);
            bool actual = false;

            try
            {
                Matrix ActualC = A + B;
            }
            catch (System.Exception e)
            {
                if (e.Message == "Размерности матриц не совпадают") actual = true;
            }

            Assert.IsTrue(actual);
        }

        [TestMethod]
        public void SubstractMatrixTest()
        {
            double[,] a = new double[3, 3] { { 8, 2, 4 }, { 14, 26, 4 }, { 7, 5, 11 } };
            double[,] b = new double[3, 3] { { 10, 2, 1 }, { 4, 6, 8 }, { 7, 3, 10 } };
            Matrix A = new Matrix(a);
            Matrix B = new Matrix(b);

            double[,] expected = new double[3, 3] { { -2, 0, 3 }, { 10, 20, -4 }, { 0, 2, 1 } };
            Matrix ExpectedC = new Matrix(expected);

            Matrix ActualC = A - B;

            for (int i = 0; i < A.GetN; i++)
                for (int j = 0; j < A.GetN; j++)
                {
                    Assert.AreEqual(ExpectedC.GetArray[i, j], ActualC.GetArray[i, j], 0.1, "Substract not calculated correctly");
                }
        }

        [TestMethod]
        public void SubstractDifferentMatrixTest()
        {
            double[,] a = new double[2, 2] { { 8, 2 }, { 4, 6 } };
            double[,] b = new double[3, 3] { { 8, 2, 4 }, { 4, 6, 8 }, { 7, 5, 10 } };
            Matrix A = new Matrix(a);
            Matrix B = new Matrix(b);
            bool actual = false;

            try
            {
                Matrix ActualC = A - B;
            }
            catch (System.Exception e)
            {
                if (e.Message == "Размерности матриц не совпадают") actual = true;
            }

            Assert.IsTrue(actual);
        }

        [TestMethod]
        public void TransposeMatrixTest()
        {
            double[,] a = new double[3, 3] { { 8, 2, 4 }, { 14, 26, 4 }, { 7, 5, 11 } };
            Matrix A = new Matrix(a);

            double[,] expected = new double[3, 3] { { 8, 14, 7 }, { 2, 26, 5 }, { 4, 4, 11 } };
            Matrix ExpectedC = new Matrix(expected);

            Matrix ActualC = Matrix.Transpose(A);

            for (int i = 0; i < ExpectedC.GetN; i++)
                for (int j = 0; j < ExpectedC.GetM; j++)
                {
                    Assert.AreEqual(ExpectedC.GetCell(i, j), ActualC.GetCell(i, j), 0.1, "Substract not calculated correctly");
                }
        }
    }

    [TestClass]
    public class SquareMatrixClassTests
    {
        [TestMethod]
        public void CreateNegativDimensionMatrix()
        {
            bool actual = false;
            try
            {
                SquareMatrix ActualC = new SquareMatrix(-1);
            }
            catch (System.Exception e)
            {
                if (e.Message == "Размерность матрицы не может быть отрицательной") actual = true;
            }

            Assert.IsTrue(actual, "Создана матрица отрицательной размерности");
        }

        [TestMethod]
        public void CreateZeroDimensionMatrix()
        {
            bool actual = false;
            try
            {
                SquareMatrix ActualC = new SquareMatrix(0);
            }
            catch (System.Exception e)
            {
                if (e.Message == "Размерность матрицы не может быть нулевой") actual = true;
            }

            Assert.IsTrue(actual, "Создана матрица нулевой размерности");
        }

        [TestMethod]
        public void MultiplicationMatrixTest()
        {
            double[,] a = new double[3, 3] { { 8, 2, 4 }, { 4, 6, 8 }, { 7, 5, 10 } };
            double[,] b = new double[3, 3] { { 8, 2, 4 }, { 4, 6, 8 }, { 7, 5, 10 } };
            SquareMatrix A = new SquareMatrix(a);
            SquareMatrix B = new SquareMatrix(b);

            double[,] expected = new double[3, 3] { { 100, 48, 88 }, { 112, 84, 144 }, { 146, 94, 168 } };
            SquareMatrix ExpectedC = new SquareMatrix(expected);

            SquareMatrix ActualC = A * B;

            for (int i = 0; i < A.GetN; i++)
                for (int j = 0; j < A.GetN; j++)
                {
                    Assert.AreEqual(ExpectedC.GetArray[i, j], ActualC.GetArray[i, j], 0.1, "Multiplication not calculated correctly");
                }
        }

        [TestMethod]
        public void AdditionMatrixTest()
        {
            double[,] a = new double[3, 3] { { 8, 2, 4 }, { 4, 6, 8 }, { 7, 5, 10 } };
            double[,] b = new double[3, 3] { { 8, 2, 4 }, { 4, 6, 8 }, { 7, 5, 10 } };
            SquareMatrix A = new SquareMatrix(a);
            SquareMatrix B = new SquareMatrix(b);
            
            double[,] expected = new double[3, 3] { { 16, 4, 8 }, { 8, 12, 16 }, { 14, 10, 20 } };
            SquareMatrix ExpectedC = new SquareMatrix(expected);

            SquareMatrix ActualC = A + B;

            for (int i = 0; i < A.GetN; i++)
                for (int j = 0; j < A.GetN; j++)
                {
                    Assert.AreEqual(ExpectedC.GetArray[i, j], ActualC.GetArray[i, j], 0.1, "Addition not calculated correctly");
                }
        }

        [TestMethod]
        public void AdditionDifferentMatrixTest()
        {
            double[,] a = new double[2, 2] { { 8, 2 }, { 4, 6 } };
            double[,] b = new double[3, 3] { { 8, 2, 4 }, { 4, 6, 8 }, { 7, 5, 10 } };
            SquareMatrix A = new SquareMatrix(a);
            SquareMatrix B = new SquareMatrix(b);
            bool actual = false;

            try
            {
                SquareMatrix ActualC = A + B;
            }
            catch (System.Exception e)
            {
                if (e.Message == "Невозможно сложить матрицы различной размерности") actual = true;
            }

            Assert.IsTrue(actual);
        }

        [TestMethod]
        public void SubstractMatrixTest()
        {
            double[,] a = new double[3, 3] { { 8, 2, 4 }, { 14, 26, 4 }, { 7, 5, 11 } };
            double[,] b = new double[3, 3] { { 10, 2, 1 }, { 4, 6, 8 }, { 7, 3, 10 } };
            SquareMatrix A = new SquareMatrix(a);
            SquareMatrix B = new SquareMatrix(b);

            double[,] expected = new double[3, 3] { { -2, 0, 3 }, { 10, 20, -4 }, { 0, 2, 1 } };
            SquareMatrix ExpectedC = new SquareMatrix(expected);

            SquareMatrix ActualC = A - B;

            for (int i = 0; i < A.GetN; i++)
                for (int j = 0; j < A.GetN; j++)
                {
                    Assert.AreEqual(ExpectedC.GetArray[i, j], ActualC.GetArray[i, j], 0.1, "Substract not calculated correctly");
                }
        }

        [TestMethod]
        public void GetInverseMatrixTest()
        {
            double[,] a = new double[4, 4] { { 8, 2, 4, 3 }, { 14, 26, 4, 12 }, { 7, 5, 11, 7 }, { 8, 2, 8, 9 } };
            SquareMatrix A = new SquareMatrix(a);

            SquareMatrix ExpectedC = new SquareMatrix(SquareMatrix.GetUnitMatrix(A.GetN).GetArray);

            SquareMatrix B = A.GetInverseMatrix();
            SquareMatrix ActualC = A * B;

            for (int i = 0; i < A.GetN; i++)
                for (int j = 0; j < A.GetN; j++)
                {
                    Assert.AreEqual(ExpectedC.GetArray[i, j], ActualC.GetArray[i, j], 0.001, "Вычисление обратной матрицы выполнено неправильно");
                }
        }

        [TestMethod]
        public void SubstractDifferentMatrixTest()
        {
            double[,] a = new double[2, 2] { { 8, 2 }, { 4, 6 } };
            double[,] b = new double[3, 3] { { 8, 2, 4 }, { 4, 6, 8 }, { 7, 5, 10 } };
            SquareMatrix A = new SquareMatrix(a);
            SquareMatrix B = new SquareMatrix(b);

            bool actual = false;
            try
            {
                SquareMatrix ActualC = A - B;
            }
            catch (System.Exception e)
            {
                if (e.Message == "Невозможно вычесть матрицы различной длины") actual = true;
            }

            Assert.IsTrue(actual);
        }

        [TestMethod]
        public void GetMaxEigenvalueTest()
        {
            double[,] a = new double[3, 3] { { 8, 2, 4 }, { 4, 6, 8 }, { 7, 5, 10 } };
            SquareMatrix A = new SquareMatrix(a);

            double expected = 18.143;
            double actual = SquareMatrix.GetMaxEigenvalue(A, 0.01);
            Assert.AreEqual(expected, actual, 0.01, "Max eigenvalue not calculated correctly");
        }

        [TestMethod]
        public void TransposeMatrixTest()
        {
            double[,] a = new double[3, 3] { { 8, 2, 4 }, { 14, 26, 4 }, { 7, 5, 11 } };
            SquareMatrix A = new SquareMatrix(a);

            double[,] expected = new double[3, 3] { { 8, 14, 7 }, { 2, 26, 5 }, { 4, 4, 11 } };
            SquareMatrix ExpectedC = new SquareMatrix(expected);

            SquareMatrix ActualC = SquareMatrix.Transpose(A);

            for (int i = 0; i < ExpectedC.GetN; i++)
                for (int j = 0; j < ExpectedC.GetN; j++)
                {
                    Assert.AreEqual(ExpectedC.GetCell(i, j), ActualC.GetCell(i, j), 0.1, "Substract not calculated correctly");
                }
        }

        [TestMethod]
        public void IsSymmetryMatrixTest()
        {
            double[,] a = new double[3, 3] { { 8, 2, 4 }, { 2, 26, 1 }, { 4, 1, 11 } };
            double[,] b = new double[3, 3] { { 8, 2, 4 }, { 14, 26, 4 }, { 7, 5, 11 } };
            SquareMatrix A = new SquareMatrix(a);
            SquareMatrix B = new SquareMatrix(b);

            bool ActualSymmetryA = A.IsSymmetry();
            Assert.AreEqual(true, ActualSymmetryA, "Симметричная матрица не распознана");

            bool ActualSymmetryB = B.IsSymmetry();
            Assert.AreEqual(false, ActualSymmetryB, "Несимметричная матрица не распознана");
        }

        [TestMethod]
        public void ILU_Test()
        {
            double[,] a = new double[7, 7] { { 9, 0, 0, 3, 1, 0, 1 }, { 0, 11, 2, 1, 0, 0, 2 }, { 0, 1, 10, 2, 0, 0, 0 }, { 2, 1, 2, 9, 1, 0, 0 }, { 1, 0, 0, 1, 12, 0, 1 }, { 0, 0, 0, 0, 0, 8, 0 }, { 2, 2, 0, 0, 3, 0, 8 } };
            SquareMatrix MatrixA = new SquareMatrix(a);
            SquareMatrix ActualL = new SquareMatrix();
            SquareMatrix ActualU = new SquareMatrix();
            SquareMatrix ActualR = new SquareMatrix();

            MatrixA.ILU(out ActualL, out ActualU, out ActualR);

            SquareMatrix LU = ActualL * ActualU + ActualR;
            double expected = 93;
            double actual = 1;
            Assert.AreEqual(expected, actual, 0.01, "Scalar Multiplication not calculated correctly");
            //for (int i = 0; i < ExpectedL.GetN; i++)
            //{
            //    Assert.AreEqual(ExpectedL.GetArray[i,i], ActualL.GetVector[i], 0.01, "Solution not calculated correctly");
            //}
        }
    }

    [TestClass]
    public class SLEQClassTests
    {
        [TestMethod]
        public void MRES_3x3()
        {
            double[,] A = new double[3, 3] { { 8, 2, 4 }, { 4, 6, 8 }, { 7, 5, 10 } };
            SquareMatrix MatrixA = new SquareMatrix(A);
            Vector F = new Vector(2, 1, 3);
            Vector X0 = new Vector(1, 1, 1);
            Vector ExpectedX = new Vector(0.154, -0.577, 0.481);

            Vector ActualX = SLEQ.MRES(MatrixA, F, X0);

            for (int i = 0; i < ExpectedX.GetLength; i++)
            {
                Assert.AreEqual(ExpectedX.GetVector[i], ActualX.GetVector[i], 0.01, "Solution not calculated correctly");
            }
        }

        [TestMethod]
        public void BiCGStab_3x3()
        {
            double[,] A = new double[3, 3] { { 8, 2, 4 }, { 4, 6, 8 }, { 7, 5, 10 } };
            SquareMatrix MatrixA = new SquareMatrix(A);
            Vector F = new Vector(2, 1, 3);
            Vector X0 = new Vector(1, 1, 1);
            Vector ExpectedX = new Vector(0.154, -0.577, 0.481);

            Vector ActualX = SLEQ.BiCGStab(MatrixA, F, X0, 100, 0.1);

            for (int i = 0; i < ExpectedX.GetLength; i++)
            {
                Assert.AreEqual(ExpectedX.GetVector[i], ActualX.GetVector[i], 0.01, "Solution not calculated correctly");
            }
        }

        [TestMethod]
        public void InverseMatrixBiCGStab_3x3()
        {
            double[,] A = new double[3, 3] { { 8, 2, 4 }, { 4, 6, 8 }, { 7, 5, 10 } };
            SquareMatrix MatrixA = new SquareMatrix(A);
            SquareMatrix Expected = SquareMatrix.GetUnitMatrix(MatrixA.GetN);
            SquareMatrix Actual = MatrixA * SLEQ.InverseMatrixBiCGStab(MatrixA);

            for (int i = 0; i < MatrixA.GetN; i++)
            {
                for (int j = 0; j < MatrixA.GetN; j++)
                {
                    Assert.AreEqual(Expected.GetArray[i, j], Actual.GetArray[i, j], 0.01, "Обратная матрица не вычисляется неверно");
                }   
            }
        }

        [TestMethod]
        public void InverseMatrixMRES_3x3()
        {
            double[,] A = new double[3, 3] { { 8, 2, 4 }, { 4, 6, 8 }, { 7, 5, 10 } };
            SquareMatrix MatrixA = new SquareMatrix(A);
            SquareMatrix Expected = SquareMatrix.GetUnitMatrix(MatrixA.GetN);
            SquareMatrix Actual = MatrixA * SLEQ.InverseMatrixMRES(MatrixA);

            for (int i = 0; i < MatrixA.GetN; i++)
            {
                for (int j = 0; j < MatrixA.GetN; j++)
                {
                    Assert.AreEqual(Expected.GetArray[i, j], Actual.GetArray[i, j], 0.01, "Обратная матрица не вычисляется неверно");
                }
            }
        }
    }
}
