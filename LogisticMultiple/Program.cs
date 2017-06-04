using System;
using System.IO;
using System.Collections.Generic;
using System.Linq;
using Accord.Math;
using Accord.IO;
using Accord.Statistics;

namespace LogisticMultiple
{
    class Program
    {
        static void Main(string[] args)
        {
            double[] ti = Vector.Interval(0.0, 80.0, 2.0);
            string distributionType = "gamma";
            double[] pInit = new[] { 15.0, 30.0 };

            StatAnalysis statAnalysis = new StatAnalysis();
            Tuple<double[], double> qAndJopt = statAnalysis.RunStatAnalysis(ti, distributionType, pInit);

            Console.WriteLine(String.Format("q = [ {0}, {1} ]", qAndJopt.Item1[0], qAndJopt.Item1[1]));
            Console.WriteLine(String.Format("Jopt = {0}", qAndJopt.Item2));
            Console.ReadLine();
        }
    }
}
