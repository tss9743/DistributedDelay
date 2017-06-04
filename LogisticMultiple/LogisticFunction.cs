using System;
using Accord;
using System.Linq;
using Accord.Math;
using System.Collections.Generic;
using System.Threading.Tasks;

namespace LogisticMultiple
{
    class LogisticFunction
    {
        const double a = 0.75;
        const double K = 17.5;
        const int r = 10;
        const int N = 32;
        const double delt = 0.01;
        
        HelperFunctions helper;
        double[] ti;
        int MCN;
        double[] B0;
        double[,] Q = Matrix.Zeros(N + 1, N + 1);
        double[,] H = Matrix.Zeros(N + 1, N + 1);
        double[] xN0;
        
        public LogisticFunction(double[] _ti, int _MCN)
        {
            // Time axis
            ti = _ti;
            // MonteCarlo Number
            MCN = _MCN;

            // Defining our basis elements
            helper = new HelperFunctions(r, N);
            B0 = helper.B(0);

            // Creating Q Matrix
            for (int i = 0; i < Q.Rows(); i++)
                for (int j = 0; j < Q.Columns(); j++)
                    if (j == i)
                        Q[i,j] = 2.0 / 3.0;
                    else if (j - 1 == i && i != Q.Rows())
                        Q[i,j] = 1.0 / 6.0;
                    else if (j == i - 1 && j != Q.Columns())
                        Q[i,j] = 1.0 / 6.0;

            Q[0,0] = (double)N / (double)r + 1.0 / 3.0;
            Q[N,N] = 1.0 / 3.0;
            Q = Q.Multiply((double)r / (double)N);

            // Solving initial history
            double eta = 0.1;
            double phi = 0.1;

            double[] hInt = Matrix.Zeros(1, N + 1).GetRow(0);
            hInt = hInt.Add((double)r / (double)N);
            hInt[0] = (double)r / (2.0 * (double)N);
            hInt[N] = (double)r / (2.0 * (double)N);
            hInt = hInt.Multiply(phi);
            double[] hN = B0.Multiply(eta).Add(hInt);
            xN0 = Q.Solve(hN);

            // Solving H matrix
            for (int i = 0; i < H.Rows(); i++)
                for (int j = 0; j < H.Columns(); j++)
                    if (j - 1 == i && i != H.Rows())
                        H[i,j] = -1.0 / 2.0;
                    else if (j == i - 1 && j != H.Columns())
                        H[i,j] = 1.0 / 2.0;

            H[0,0] = 1.0 / 2.0;
            H[N,N] = -1.0 / 2.0;
        }

        public Tuple<double[], double[]> RungeKutta(string distributionType, double[] q)
        {
            Console.WriteLine("Started a RK...");
            double[] t = Vector.Interval(0, ti.Last(), delt);
            int nt = t.Length();

            // Check for bad dist params
            //if (q.Any(qe => qe < 0))
            //{
            //    double[] nanArray = Vector.Zeros(ti.Length()).Divide(0.0);
            //    return Tuple.Create(nanArray, nanArray);
            //}

            DistributionSampler dist = new DistributionSampler(distributionType, q);

            List<double[]> x = new List<double[]>();
            Parallel.For(0, MCN, i =>
            {
                double tau = dist.SampleDistribution();
                double[] Bt = helper.B(-tau);
                double[,] xN = Matrix.Zeros(N + 1, nt);
                xN.SetColumn(0, xN0);
                for (int k = 0; k < nt - 1; k++)
                {
                    double[] k1 = delt.Multiply(LogisticRHS(xN.GetColumn(k), Bt));
                    double[] k2 = delt.Multiply(LogisticRHS(xN.GetColumn(k).Add(k1.Multiply(1.0 / 2.0)), Bt));
                    double[] k3 = delt.Multiply(LogisticRHS(xN.GetColumn(k).Add(k2.Multiply(1.0 / 2.0)), Bt));
                    double[] k4 = delt.Multiply(LogisticRHS(xN.GetColumn(k).Add(k3), Bt));
                    k1 = k1.Multiply(1.0 / 6.0);
                    k2 = k2.Multiply(1.0 / 3.0);
                    k3 = k3.Multiply(1.0 / 3.0);
                    k4 = k4.Multiply(1.0 / 6.0);
                    double[] rkTerm = k1.Add(k2).Add(k3).Add(k4);
                    xN.SetColumn(k + 1, xN.GetColumn(k).Add(rkTerm));
                }
                x.Add(xN.GetRow(0));
            });
            Console.WriteLine("Finished a RK...");
            Console.WriteLine("");
            return CalcMeanAndVariance(x.ToArray());
        }

        private double[] LogisticRHS(double[] xN, double[] Bt)
        {
            double x = xN[0];
            double xt = Bt.Dot(xN);

            double dxN = a * x * (1.0 - (xt / K));
            return Q.Solve(B0.Multiply(dxN).Add(H.Dot(xN)));
        }

        private Tuple<double[], double[]> CalcMeanAndVariance(double[][] x)
        {
            double[,] xf = Matrix.Zeros(MCN, ti.Length());
            for (int t = 0; t < ti.Length(); t++)
            {
                double time = ti[t] / delt;
                time = Math.Floor(time);
                for (int i = 0; i < MCN; i++)
                {
                    xf[i, t] = x[i][(int)time];
                }
            }

            double[] meanF = Vector.Zeros(ti.Length());
            double[] varF = Vector.Zeros(ti.Length());
            for (int t = 0; t < ti.Length(); t++)
            {
                double[] values = Vector.Zeros(MCN);
                for (int i = 0; i < MCN; i++)
                {
                    values[i] = xf[i, t];
                }
                meanF[t] = Accord.Statistics.Measures.Mean(values);
                varF[t] = Math.Pow(Accord.Statistics.Measures.StandardDeviation(values), 2.0);
            }
            return Tuple.Create(meanF, varF);
        }
    }
}
