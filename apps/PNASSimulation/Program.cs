using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ASELib;

namespace PNASSimulation
{
    class Program
    {
        delegate int GrowthCount();
        delegate List<int> DifferentiationCount(); // List ordered the same as the allCompartments list.

        class Compartment
        {
            public Compartment(string name_, int initialN)
            {
                name = name_;
                N = initialN;
            }

            public void setFunctions(GrowthCount growthCount_, DifferentiationCount differentiationCount_)
            {
                growthCount = growthCount_;
                differentiationCount = differentiationCount_;

            }

            public int getN()
            {
                return N;
            }

            public readonly string name;
            GrowthCount growthCount;
            DifferentiationCount differentiationCount;
            int N;
        }

        static List<int> FixedPercentDifferentiation(Compartment thisCompartment, int nCompartments, int outputCompartment, double percentage)
        {
            var retVal = new List<int>();
            for (int i = 0; i < nCompartments; i++)
            {
                if (i == outputCompartment)
                {
                    retVal.Add((int)(thisCompartment.getN() * percentage)); // Maybe should round better
                } else
                {
                    retVal.Add(0);
                }
            }

            return retVal;
        }

        static void Main(string[] args)
        {
            var listOfCompartments = new List<Compartment>();

            listOfCompartments.Add(new Compartment("A", 1000));
            listOfCompartments.Add(new Compartment("B", 2000));
            listOfCompartments.Add(new Compartment("C", 4000));
            listOfCompartments.Add(new Compartment("D", 8000));

            int nCompartments = listOfCompartments.Count();

            var compartments = listOfCompartments.GroupByToDictUnique(_ => _.name);

            compartments["A"].setFunctions(() => compartments["A"].getN(), () => FixedPercentDifferentiation(compartments["A"], nCompartments, 1 /*B*/, 1.0));    // 100% Differentiation because it's based on the old count, which doubles, so 100% of the old ones



        }
    }
}
