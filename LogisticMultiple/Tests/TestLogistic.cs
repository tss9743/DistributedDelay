using System;
using System.IO;
using Accord;
using NUnit.Framework;
using System.Linq;
using Accord.Math;
using Accord.IO;
using System.Collections.Generic;
using System.Threading.Tasks;

namespace LogisticMultiple.Tests
{
    [TestFixture]
    public class LogisticTests
    {
        [Test]
        public void TestOptimalMcn()
        {
            double[] ti = Vector.Interval(0.0, 80.0, 2.0);
            int[] MCNs = new[] { 100, 200, 300, 400, 500 };
            string distributionType = "gamma";
            double[] q = new[] { 10.0, 20.0 };

            int fileNum = 0;
            foreach (var MCN in MCNs)
            {
                LogisticFunction logistic = new LogisticFunction(ti, MCN);
                Tuple<double[], double[]> meanAndVar = logistic.RungeKutta(distributionType, q);

                string path1 = String.Format(@"C:\Users\Chase\Desktop\mean{0}.txt", fileNum);
                TextWriter tw1 = File.CreateText(path1);
                CsvWriter writer1 = new CsvWriter(tw1, ' ');
                writer1.Write(meanAndVar.Item1.ToMatrix());
                tw1.Close();

                fileNum++;
            }
        }
    }
}