using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;

namespace QNDMannWhitney
{
    class Program
    {
        class MannWhitneyItem : IComparer<MannWhitneyItem>
        {
            public double value;
            public bool isMyeloid;

            public MannWhitneyItem(double value_, bool isMyeloid_)
            {
                value = value_;
                isMyeloid = isMyeloid_;
            }

            public int Compare(MannWhitneyItem one, MannWhitneyItem two)
            {
                return one.value.CompareTo(two.value);
            }

        }
        static void Main(string[] args)
        {
            double[] neutrophils = { .67, 1, .78, .54, .77, .33, .75, .58, .28, .83 };
            double[] t_cells = { 0, 0, 0, .12, .21, .1, .03, .03, .5, .025, 0, 0 };
            double[] myeloid_blasts = { .93, .83, .64, .73, .75, .82, .6, .67, .87, .7, .88, .98, .55, .15, .89, .95 };
            double[] monocytes = { .9, .84, .61, .73, .95, .64, .79 };


            List<MannWhitneyItem> items = new List<MannWhitneyItem>();
            neutrophils.ToList().ForEach(x => items.Add(new MannWhitneyItem(x, true)));
            t_cells.ToList().ForEach(x => items.Add(new MannWhitneyItem(x, false)));

            bool enoughData;
            bool reversed;
            double nFirstGroup, nSecondGroup, U, z; 
            var p = ASETools.MannWhitney<MannWhitneyItem>.ComputeMannWhitney(items, items[0], x => x.isMyeloid, x => x.value, out enoughData, out reversed, out nFirstGroup, out nSecondGroup, out U, out z);

            Console.WriteLine("P neutrophils vs. T " + p);

            items = new List<MannWhitneyItem>();
            neutrophils.ToList().ForEach(x => items.Add(new MannWhitneyItem(x, true)));
            myeloid_blasts.ToList().ForEach(x => items.Add(new MannWhitneyItem(x, false)));

            p = ASETools.MannWhitney<MannWhitneyItem>.ComputeMannWhitney(items, items[0], x => x.isMyeloid, x => x.value, out enoughData, out reversed, out nFirstGroup, out nSecondGroup, out U, out z);

            Console.WriteLine("P neutrophils vs. myeloid blasts " + p);

            items = new List<MannWhitneyItem>();
            monocytes.ToList().ForEach(x => items.Add(new MannWhitneyItem(x, true)));
            myeloid_blasts.ToList().ForEach(x => items.Add(new MannWhitneyItem(x, false)));

            p = ASETools.MannWhitney<MannWhitneyItem>.ComputeMannWhitney(items, items[0], x => x.isMyeloid, x => x.value, out enoughData, out reversed, out nFirstGroup, out nSecondGroup, out U, out z);

            Console.WriteLine("P monocytes vs. myeloid blasts " + p);

            items = new List<MannWhitneyItem>();
            monocytes.ToList().ForEach(x => items.Add(new MannWhitneyItem(x, true)));
            t_cells.ToList().ForEach(x => items.Add(new MannWhitneyItem(x, false)));

            p = ASETools.MannWhitney<MannWhitneyItem>.ComputeMannWhitney(items, items[0], x => x.isMyeloid, x => x.value, out enoughData, out reversed, out nFirstGroup, out nSecondGroup, out U, out z);

            Console.WriteLine("P monocytes vs. t cells " + p);

            items = new List<MannWhitneyItem>();
            myeloid_blasts.ToList().ForEach(x => items.Add(new MannWhitneyItem(x, true)));
            t_cells.ToList().ForEach(x => items.Add(new MannWhitneyItem(x, false)));

            p = ASETools.MannWhitney<MannWhitneyItem>.ComputeMannWhitney(items, items[0], x => x.isMyeloid, x => x.value, out enoughData, out reversed, out nFirstGroup, out nSecondGroup, out U, out z);

            Console.WriteLine("P myeloid blasts vs. t cells " + p);

            items = new List<MannWhitneyItem>();
            neutrophils.ToList().ForEach(x => items.Add(new MannWhitneyItem(x, true)));
            monocytes.ToList().ForEach(x => items.Add(new MannWhitneyItem(x, false)));

            p = ASETools.MannWhitney<MannWhitneyItem>.ComputeMannWhitney(items, items[0], x => x.isMyeloid, x => x.value, out enoughData, out reversed, out nFirstGroup, out nSecondGroup, out U, out z);

            Console.WriteLine("P neutrophils vs. t monocytes " + p);
        }
    }
}
