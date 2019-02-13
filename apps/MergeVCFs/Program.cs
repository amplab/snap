using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ASELib;
using System.Diagnostics;

namespace MergeVCFs
{
    class Program
    {
        class VCFSamplePair
        {
            StreamReader 
        }
        static void Main(string[] args)
        {
            if (args.Count() < 4 || args.Count() % 1 == 1)
            {
                Console.WriteLine("usage: MergeVCFs [vcfFilename SampleToUse]");
                Console.WriteLine("you must have at least two input VCF/sample pairs.");
                return;
            }


        }
    }
}
