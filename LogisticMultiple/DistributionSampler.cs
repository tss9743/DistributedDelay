using System;
using Accord;
using System.Linq;
using Accord.Math;
using Accord.Statistics;

namespace LogisticMultiple
{
    class DistributionSampler
    {
        Accord.Statistics.Distributions.Univariate.UnivariateContinuousDistribution distribution;

        public DistributionSampler(string distributionType, double[] q)
        {
            if (distributionType == "gamma")
            {
                distribution = new Accord.Statistics.Distributions.Univariate.GammaDistribution(1/q[0], q[1]);
            }
        }

        public double SampleDistribution()
        {
            return distribution.Generate(1)[0];
        }
    }
}
