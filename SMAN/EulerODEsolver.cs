using System;

/*********************************  
* Автор: Николаев Виталий
* 
* Версия: 1.0
*    
**********************************/

//Compiler version 4.0, .NET Framework 4.5
// есть лишнее повторение, когда требуется удвоить шаг

namespace SMAN
{

    public class ODE
    {
        public delegate double FunkXY(double x, double y);

        public static double[] Euler(double y0, double a, double b, int p, FunkXY D)
        // Явный метод Эйлера первого порядка точности
        {
            if (b <= a) { throw new ArgumentException("Отрезок интегрирования должен быть упорядочен"); }
            if (p <= 0) { throw new ArgumentException("Количество точек должно быть более неотрицательным"); }
            if (p < 2) { throw new ArgumentException("Слишком малое количество возвращаемых точек"); }

            int n = p;
            int n2 = 2 * (n - 1) + 1;

            double[] y;
            double[] y2n = new double[n2];
            double dx = 2 * (b - a) / (n - 1);
            int count = -1;

            bool end = false;
            while (!end)
            {
                count += 1;
                dx = dx / 2;
                n = Convert.ToInt32((b - a) / dx) + 1;
                y = new double[n];
                y[0] = y0;
                for (int i = 1; i < n; i++)
                {
                    y[i] = y[i - 1] + dx * D(a + dx * (i - 1), y[i - 1]);
                }

                count += 1;
                dx = dx / 2;
                n2 = Convert.ToInt32((b - a) / dx) + 1;
                y2n = new double[n2];
                y2n[0] = y0;
                for (int i = 1; i < n2; i++)
                {
                    y2n[i] = y2n[i - 1] + dx * D(a + dx * (i - 1), y2n[i - 1]);
                }

                end = true;
                for (int i = 1; i < n; i++)
                {
                    end = end && (Math.Abs(y[i] - y2n[2 * i]) / 2 < 1e-4);
                }
            }

            int k = Convert.ToInt32(Math.Pow(2, count));
            count = 0;
            double[] answer = new double[p];
            for (int i = 1; i < p; i++)
            {
                count = count + k;
                answer[i] = y2n[count];
            }

            return answer;
        }

        public static double[] PCEuler(double y0, double a, double b, int p, FunkXY D)
        // Predictor–corrector method Euler
        // Модифицированный метод Эйлера (предиктор - корректор)
        {
            if (b < a) { throw new ArgumentException("Отрезок интегрирования должен быть упорядочен"); }
            if (b == a) { throw new ArgumentException("Пределы интегрирования не должны совпадать"); }
            if (p <= 0) { throw new ArgumentException("Количество точек должно быть неотрицательным"); }
            if (p < 2) { throw new ArgumentException("Слишком малое количество возвращаемых точек"); }

            int n = p;
            int n2 = 2 * (n - 1) + 1;

            double[] y;
            double[] y2n = new double[n2];
            double dx = 2 * (b - a) / (n - 1);
            int count = -1;

            bool end = false;
            while (!end)
            {
                count += 1;
                dx = dx / 2;
                n = Convert.ToInt32((b - a) / dx) + 1;
                y = new double[n];
                y[0] = y0;
                for (int i = 1; i < n; i++)
                {
                    y[i] = y[i - 1] + dx * D(a + dx * (i - 1), y[i - 1]);
                    y[i] = y[i - 1] + dx * D(a + dx * i, y[i]);
                }

                count += 1;
                dx = dx / 2;
                n2 = Convert.ToInt32((b - a) / dx) + 1;
                y2n = new double[n2];
                y2n[0] = y0;
                for (int i = 1; i < n2; i++)
                {
                    y2n[i] = y2n[i - 1] + dx * D(a + dx * (i - 1), y2n[i - 1]);
                    y2n[i] = y2n[i - 1] + dx * D(a + dx * i, y2n[i]);
                }

                end = true;
                for (int i = 1; i < n; i++)
                {
                    end = end && (Math.Abs(y[i] - y2n[2 * i]) / 2 < 1e-4);
                }
            }

            int k = Convert.ToInt32(Math.Pow(2, count));
            count = 0;
            double[] answer = new double[p];
            for (int i = 1; i < p; i++)
            {
                count = count + k;
                answer[i] = y2n[count];
            }

            return answer;
        }

        public double[] Adams(double y0, double a, double b, int p, FunkXY D)
        {
            double[] с = new double[p];
            return с;
        }

        public double[] BDF(double y0, double a, double b, int p, FunkXY D)
        {
            double[] с = new double[p];
            return с;
        }
    }

}
    
    