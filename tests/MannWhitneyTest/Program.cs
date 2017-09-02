using System;
using System.Collections.Generic;
using System.Linq;
using ASELib;
using MathNet.Numerics;

namespace MannWhitneyTest
{
	class ComparableElement: IComparer<ComparableElement>
	{
		public double value;
		public bool group;

		public ComparableElement(double value_, bool group_)
		{
			value = value_;
			group = group_;
		}

		public int Compare(ComparableElement a, ComparableElement b)
		{
			if (a.value > b.value) return 1;
			if (a.value < b.value) return -1;
			return 0;
		}
	}

	class Program
	{


		static void Main(string[] args)
		{
			double[] values =
			{   21,22.8,10.4,15.2,19.2,30.4,21.4
			};

			int[] groups =
			{ 1,0,0,0,0,1,1 };

			List<ComparableElement> test1 = new List<ComparableElement>();

			for (var i = 0; i < values.Count(); i++)
			{
				test1.Add(new ComparableElement(values[i], groups[i] == 1));
			}

			ASETools.MannWhitney<ComparableElement>.GetValue getValue = new ASETools.MannWhitney<ComparableElement>.GetValue(m => m.value);
			ASETools.MannWhitney<ComparableElement>.WhichGroup whichGroup = new ASETools.MannWhitney<ComparableElement>.WhichGroup(m => m.group);

			bool enoughData;
			bool reversed;
			double nFirstGroup;
			double nSecondGroup;
			double U;
			double z;

			var result = ASETools.MannWhitney<ComparableElement>.ComputeMannWhitney(test1, test1[0], whichGroup, getValue,
				out enoughData, out reversed, out nFirstGroup, out nSecondGroup, out U, out z);

			if (U != 2)
			{
				throw new ArgumentException("Mann Whitney statistic differs from R package");
			}

			List<ComparableElement> v1 = new List<ComparableElement>();
			List<ComparableElement> v2 = new List<ComparableElement>();

			for (int i = 0; i < 10000; i++)
			{
				v1.Add(new ComparableElement(MathNet.Numerics.Distributions.Normal.Sample(2, 1.0), true));
				v2.Add(new ComparableElement(MathNet.Numerics.Distributions.Normal.Sample(2, 4.0), false));
			}

			List<ComparableElement> test2 = new List<ComparableElement>();
			test2.AddRange(v1);
			test2.AddRange(v2);

			result = ASETools.MannWhitney<ComparableElement>.ComputeMannWhitney(test2, test2[0], whichGroup, getValue,
				out enoughData, out reversed, out nFirstGroup, out nSecondGroup, out U, out z);

			// This result should not be significant
			if (result < 0.1)
			{
				throw new ArgumentException("distributions should not differ. P-value too low");
			}

            var random = new Random();
            var pValueHistogram = new ASETools.Histogram("p values for uniform distribution");
            for (int i = 0; i < 100000; i++)
            {
                int group0Size = random.Next(50000) + 10;
                int group1Size = random.Next(5000) + 10;

                var randomValues = new List<ComparableElement>(); 
                for (int j = 0; j < group0Size + group1Size; j++)
                {
                    randomValues.Add(new ComparableElement(random.NextDouble(), j < group0Size));
                }

                var p = ASETools.MannWhitney<ComparableElement>.ComputeMannWhitney(randomValues, randomValues[0], x => x.group, x => x.value, out enoughData, out reversed, out nFirstGroup, out nSecondGroup, out U, out z);

                if (enoughData)
                {
                    pValueHistogram.addValue(p);
                }

                if (i % 1000 == 999) Console.Write(".");
            }

            Console.WriteLine();

            var outputFile = ASETools.CreateStreamWriterWithRetry("MannWhitneyTestPValueHistograms.txt");

            outputFile.WriteLine("p value histogram (should be flat line because the input's a uniform distribution)");
            outputFile.WriteLine(ASETools.HistogramResultLine.Header());
            pValueHistogram.ComputeHistogram(0, 1, .01).ToList().ForEach(x => outputFile.WriteLine(x));

            //
            // Now generate a p-value histogram for a distribution more like what we're measuring with ASE.  The values are means of randomly sized uniform distributions.
            //
            pValueHistogram = new ASETools.Histogram("");
            for (int iteration = 0; iteration < 100000; iteration++)
            {

            }

            outputFile.Close();
		}

	}
}
