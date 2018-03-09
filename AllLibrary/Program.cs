using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using Dcoder;

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
            string a = Convert.ToString(Mathem.Integrate(0, 1, fun));
            string b = Convert.ToString(Mathem.Integrate(0, 1, 1e-12, fun));
            Console.WriteLine(String.Concat("Integral f(x)=e^x from 0 to 1 = ", a));
            Console.WriteLine(String.Concat("Integral f(x)=e^x from 0 to 1 = ", b));
        }
    }
}
