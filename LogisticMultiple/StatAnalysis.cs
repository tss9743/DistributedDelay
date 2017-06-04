using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using Accord.Math;
using Accord.IO;
using Accord.Statistics;
using Microsoft.SolverFoundation.Common;
using Microsoft.SolverFoundation.Solvers;

namespace LogisticMultiple
{
    class StatAnalysis
    {
        LogisticFunction logisticFunction;

        public Tuple<double[], double> RunStatAnalysis(double[] ti, string distributionType, double[] pInit)
        {
            // Load test data
            //string path1 = @"C:\Users\Chase\Google Drive\Temi Research\Logistic Multiple\meanF.csv";
            //string path2 = @"C:\Users\Chase\Google Drive\Temi Research\Logistic Multiple\varF.csv";
            string path1 = @"C:\Users\Administrator\Documents\meanF.csv";
            string path2 = @"C:\Users\Administrator\Documents\varF.csv";

            CsvReader r1 = new CsvReader(path1, false);
            CsvReader r2 = new CsvReader(path2, false);

            double[] meanF = r1.ToJagged().GetRow(0);
            double[] varF = r2.ToJagged().GetRow(0);

            // Setup logistic
            int MCN = 200;
            logisticFunction = new LogisticFunction(ti, MCN);

            // Minimize f with NelderMead
            Func<Double[], Double> f = q => MinFunction(distributionType, q, meanF, varF);
            var solution = NelderMeadSolver.Solve(f, pInit, new[] { 0.0, 0.0 }, new[] { 50.0, 50.0 });

            return Tuple.Create(new[] { solution.GetValue(1), solution.GetValue(2) }, solution.GetSolutionValue(0));
        }

        double MinFunction(string distributionType, double[] q, double[] meanF, double[] varF)
        {
            Tuple<double[], double[]> meanAndVar = logisticFunction.RungeKutta(distributionType, q);
            double[] xMean = meanAndVar.Item1;
            double[] xVar = meanAndVar.Item2;

            double[] meanTerms = Elementwise.Pow(meanF.Subtract(xMean), 2.0);
            double[] varTerms = Elementwise.Pow(Elementwise.Sqrt(varF).Subtract(Elementwise.Sqrt(xVar.Add(1.0))), 2.0);

            return meanTerms.Add(varTerms).Sum();
        }
    }
}
