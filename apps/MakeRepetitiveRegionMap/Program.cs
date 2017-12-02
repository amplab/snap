using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Diagnostics;
using ASELib;

namespace MakeRepetitiveRegionMap
{
    class Program
    {
        static ASETools.Genome genome;
        static 

        static void Main(string[] args)
        {
            if (args.Count() != 3)
            {
                Console.WriteLine("usage: GenomeIndex MakeRepetitiveRegionMap inputDirectory outputFilename");
            }
            genome = new ASETools.Genome();

            var timer = new Stopwatch();
            timer.Start();
            Console.Write("Loading genome...");
            if (!genome.load(args[0]))
            {
                Console.WriteLine();
                Console.WriteLine("Unable to load genome from directory " + args[0]);
                return;
            }
            Console.WriteLine(ASETools.ElapsedTimeInSeconds(timer));

            var queue = new List<int>();
            for (int i = 1; i <= ASETools.nHumanAutosomes; i++)
            {
                queue.Add(i);
            }

            var threading = new ASETools.WorkerThreadHelper<int, int>(queue, HandleOneChromosome, null, null, 1);

            threading.run();
        } // Main

        static 
    }
}
