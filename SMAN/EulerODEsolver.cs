using System;


//Compiler version 4.0, .NET Framework 4.5


namespace SMAN
// Some methods of analysis
// Ќекоторые методы математического анализа
{

    public class ODE
    {
        public delegate double FunkXY(double x, double y);

        public static double[] EulerODEsolve(double y0, double a, double b, FunkXY D)
        // Method of Euler
        {
            int n = 5;
            int n2 = n * 2;

            double[] y = new double[n];
            double[] y2n = new double[n2];
            double dy = 1;
            double dx;

            y[0] = y0;
            //y[0] = y0;

            bool end = false;
            while (!end)
            {
                n = n * 2;
                dx = (b - a) / n;
                y = new double[n];
                for (int i = 1; i <= n; i++)
                {
                    y[i] = y[i - 1] + dx * D(a + dx * i, y[i - 1]);
                }

                dx = (b - a) / n2;
                for (int i = 1; i <= n2; i++)
                {
                    y2n[i] = y2n[i - 1] + dx * D(a + dx * i, y2n[i - 1]);
                }

                end = true;
                for (int i = 0; i <= n; i++)
                {
                    //end = end * (Math.Abs(y[i] - y2n[2 * i]) / 2 < 1e-4);
                }
            }

            return y2n;
        }
    }


}
    
