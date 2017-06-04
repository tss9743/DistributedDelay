using System;
using Accord;
using Accord.Math;

namespace LogisticMultiple
{
    class HelperFunctions
    {
        private double r;
        private double N;

        public HelperFunctions(int _r, int _N)
        {
            r = _r;
            N = _N;
        }

        private double TrianglePulse(double a, double b, double c, double x)
        {
            double tp;
            if (a <= x && x <= b)
            {
                tp = (x - a) / (b - a);
            }
            else if (b <= x && x <= c)
            {
                tp = (c - x) / (c - b);
            }
            else
            {
                tp = 0.0;
            }
            return tp;
        }

        public double HatFunction(int i, double t)
        {
            return TrianglePulse(i * (-r / N) + (-r / N), i * (-r / N), i * (-r / N) - (-r / N), t);
        }

        public double[] B(double t)
        {
            double[] b = new double[(int)N + 1];
            for (int i = 0; i < (int)N + 1; i++)
            {
                b[i] = HatFunction(i, t);
            }
            return b;
        }
    }
}
