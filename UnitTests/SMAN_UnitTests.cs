using System;
using SMAN;
using Microsoft.VisualStudio.TestTools.UnitTesting;

namespace UnitTests
{
    [TestClass]
    public class ODEClassTests
    {
        public double Const(double x, double y)
        {
            return 1;
        }

        public double Line(double x, double y)
        {
            return x;
        }

        public double Exp(double x, double y)
        {
            return Math.Pow(Math.E, x);
        }

        [TestMethod]
        public void Euler_Exp()
        {
            double expected = 1.7183;
            double actual = ODE.Euler(0, 0, 1, 10, Exp)[9];
            Assert.AreEqual(expected, actual, 0.001, "Решение найдено не верно");
        }

        [TestMethod]
        public void Euler_Const()
        {
            double expected = 1;
            double actual = ODE.Euler(0, 0, 1, 10, Const)[9];
            Assert.AreEqual(expected, actual, 0.001, "Решение найдено не верно");
        }

        [TestMethod]
        public void Euler_Line()
        {
            double expected = 0.5;
            double actual = ODE.Euler(0, 0, 1, 10, Line)[9];
            Assert.AreEqual(expected, actual, 0.001, "Решение найдено не верно");
        }

        [TestMethod]
        public void ImplicitEuler_Exp()
        {
            double expected = 1.7183;
            double actual = ODE.ImplicitEuler(0, 0, 1, 10, Exp)[9];
            Assert.AreEqual(expected, actual, 0.001, "Решение найдено не верно");
        }

        [TestMethod]
        public void ImplicitEuler_Const()
        {
            double expected = 1;
            double actual = ODE.ImplicitEuler(0, 0, 1, 10, Const)[9];
            Assert.AreEqual(expected, actual, 0.001, "Решение найдено не верно");
        }

        [TestMethod]
        public void ImplicitEuler_Line()
        {
            double expected = 0.5;
            double actual = ODE.ImplicitEuler(0, 0, 1, 10, Line)[9];
            Assert.AreEqual(expected, actual, 0.001, "Решение найдено не верно");
        }

    }


}
