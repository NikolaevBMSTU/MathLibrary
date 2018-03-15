using System;
using Matrix

 //Compiler version 4.0, .NET Framework 4.5


 namespace SMLA
 {
 	
 	public class SLEQ
 	{
 		public delegate double FunkXY(double x, double y);
 		
 		public static Vector BiCGStab(Matrix a, Vector b, int p, FunkXY D)
 		// Method of Euler
 		{
 			Vector[] x, r, p, v, t, s;
 			double[] rho, alpha, omega;
 			
 			int end = 0;
 			
 			r[0] = b - A*x[0];
 			r'=r[0];
 			rho[0]=1;
 			alpha[0]=1;
 			omega[0]=1;
 			
 			for (int k; k < 1e4; k++)
 			{
 				
 				beta[k]=(rho[k]/rho[k-1])*(alpha[k-1]/omega[k-1]);
 				
 				
 				if ()
 				{
 					end = k;
 					break;
 				}
 			}
 			
 			return x[end];
 		}
 	}
 	
 	public class Program
 	{
 		public static double fun(double x, double y)
 		{
 			return x;
 		}
 		
 		public static void Main(string[] args)
 		{
 			string a = Convert.ToString(ODESolve.ImplicitEuler(0, 0, 1, 11, fun)[10]);
 			//string b = Convert.ToString(Mathem.Integrate(0, 1, 1e-12, fun));
 			Console.WriteLine(String.Concat("Integral equation f'(x)=x from 0 to 1 = ", a));
 			//Console.WriteLine(String.Concat("Integral f(x)=e^x from 0 to 1 = ", b));
 		}
 	}
 }
    
    