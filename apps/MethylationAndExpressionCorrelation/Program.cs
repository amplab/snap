using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.Diagnostics;
using System.IO;

namespace MethylationAndExpressionCorrelation
{
    class Program
    {
        class MethylationMap
        {
            Dictionary<string, Dictionary<int, double>> methylationBetaMap = new Dictionary<string, Dictionary<int, double>>(); // contig -> locus -> beta
        }
        static void Main(string[] args)
        {
        }
    }
}
