using System;
using SMAN;

namespace AllLibrary
{
    public class Program
    {
        public static double fun(double x)
        {
            return Math.Exp(x);
        }

        public static void Main(string[] args)
        {
            string a = Convert.ToString(Integrals.Integrate(0, 1, fun));
            string b = Convert.ToString(Integrals.Integrate(0, 1, 1e-12, fun));
            Console.WriteLine(String.Concat("Integral f(x)=e^x from 0 to 1 = ", a));
            Console.WriteLine(String.Concat("Integral f(x)=e^x from 0 to 1 = ", b));
        }
    }
}
