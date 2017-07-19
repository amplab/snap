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

			List<ComparableElement> ComparableElements = new List<ComparableElement>();

			for (var i = 0; i < values.Count(); i++)
			{
				ComparableElements.Add(new ComparableElement(values[i], groups[i] == 1));
			}

			ASETools.MannWhitney<ComparableElement>.GetValue getValue = new ASETools.MannWhitney<ComparableElement>.GetValue(m => m.value);
			ASETools.MannWhitney<ComparableElement>.WhichGroup whichGroup = new ASETools.MannWhitney<ComparableElement>.WhichGroup(m => m.group);

			bool enoughData;
			bool reversed;
			double nFirstGroup;
			double nSecondGroup;
			double U;
			double z;

			var result = ASETools.MannWhitney<ComparableElement>.ComputeMannWhitney(ComparableElements, ComparableElements[0], whichGroup, getValue,
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

			List<ComparableElement> cars_est = new List<ComparableElement>();
			cars_est.AddRange(v1);
			cars_est.AddRange(v2);

			result = ASETools.MannWhitney<ComparableElement>.ComputeMannWhitney(cars_est, cars_est[0], whichGroup, getValue,
				out enoughData, out reversed, out nFirstGroup, out nSecondGroup, out U, out z);

			// This result should not be significant
			if (result < 0.1)
			{
				throw new ArgumentException("distributions should not differ. P-value too low");
			}
		}

	}
}
