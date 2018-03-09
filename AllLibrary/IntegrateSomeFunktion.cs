using System;

//Compiler version 4.0, .NET Framework 4.5
//ѕодумаю над снижением версии .NET

namespace SMAN
// Some methods of analysis
// Ќекоторые методы математического анализа
{

    public class Integrals
    {
        public delegate double FunkX(double x);

        public static double Integrate(double a, double b, FunkX y)
        // Method of Gauus-2
        {
            int n = 5;
            int n2 = n * 2;
            double dx = (b - a) / n;
            double ksi = 1 / Math.Sqrt(3);

            double dI = 1;
            double In = 0;
            double I2n = 0;

            while (dI > 1e-4)
            {
                n = n * 2;
                dx = (b - a) / n;
                In = 0;
                for (int i = 0; i <= n - 1; i++)
                {
                    In = In + dx / 2 * (y(a + dx * i + dx / 2 * (-ksi + 1)) + y(a + dx * i + dx / 2 * (ksi + 1)));
                }

                n2 = n * 2;
                dx = (b - a) / n2;
                I2n = 0;
                for (int i = 0; i <= n2 - 1; i++)
                {
                    I2n = I2n + dx / 2 * (y(a + dx * i + dx / 2 * (-ksi + 1)) + y(a + dx * i + dx / 2 * (ksi + 1)));
                }

                //I2n = In;
                dI = Math.Abs(In - I2n) / 15;

            }
            return I2n;
        }

        public static double Integrate(double a, double b, double eps, FunkX y)
        // Method of Gauus-2
        {
            int n = 5;
            int n2 = n * 2;
            double dx = (b - a) / n;
            double ksi = 1 / Math.Sqrt(3);

            double dI = 1;
            double In = 0;
            double I2n = 0;

            while (dI > eps)
            {
                n = n * 2;
                dx = (b - a) / n;
                In = 0;
                for (int i = 0; i <= n - 1; i++)
                {
                    In = In + dx / 2 * (y(a + dx * i + dx / 2 * (-ksi + 1)) + y(a + dx * i + dx / 2 * (ksi + 1)));
                }

                n2 = n * 2;
                dx = (b - a) / n2;
                I2n = 0;
                for (int i = 0; i <= n2 - 1; i++)
                {
                    I2n = I2n + dx / 2 * (y(a + dx * i + dx / 2 * (-ksi + 1)) + y(a + dx * i + dx / 2 * (ksi + 1)));
                }

                dI = Math.Abs(In - I2n) / 15;
            }
            return I2n;
        }
    }

}
    