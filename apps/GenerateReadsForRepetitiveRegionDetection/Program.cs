using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;
using System.Diagnostics;

namespace GenerateReadsForRepetitiveRegionDetection
{
    class Program
    {
        static void Main(string[] args)
        {

            if (args.Count() != 4)
            {
                Console.WriteLine("usage: GenerateReadsForRepetitiveRegionDetection snapIndex contig outputFile readLength");
                Console.WriteLine("Generates a set of reads (in SAM format) that correspond to the starting location of each place in the contig");
                Console.WriteLine("This is intended to let SNAP align them to see if there are good secondary matches that would indicate");
                Console.WriteLine("repetitive regions.");
                return;
            }

            int readLength = Convert.ToInt32(args[3]);
            string contigName = args[1];

            var genome = new ASETools.Genome();

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

            if (!genome.isAContig(contigName))
            {
                Console.WriteLine(contigName + " does not appear to be a contig.");
                return;
            }

            var outputFile = ASETools.CreateStreamWriterWithRetry(args[2]);
            if (null == outputFile)
            {
                Console.WriteLine("Unable to open output file " + args[2]);
                return;
            }

            timer.Reset();
            timer.Start();
            
            Console.Write("Generating output (1 dot/10MB): ");
            outputFile.WriteLine("@HD\tVN:1.4");
            outputFile.Write("@PG\tID:GRRD.1\tPG:GenerateReadsForRepetitiveRegionDetection\tCL:");
            for (int i = 0; i < args.Count(); i++)
            {
                outputFile.Write(args[i]);
                if (i != args.Count() - 1)
                {
                    outputFile.Write(" ");
                }
            }
            outputFile.WriteLine();

            string qual = "";
            for (int i = 0; i < readLength; i++)
            {
                qual += "2";
            }

            
            for (int locus = 1; locus <= genome.getContigLength(contigName) - readLength; locus++ )
            {
                if ((locus + 1) % 10000000 == 0)
                {
                    Console.Write(".");
                }
                string seq = "";
                for (int i = 0; i < readLength; i++)
                {
                    seq += char.ToUpper(genome.getBase(contigName, locus + i));
                }
                outputFile.WriteLine(contigName + "." + locus + "\t2\t" + contigName + "\t" + locus + "\t70\t" + readLength + "M\t*\t*\t" + readLength + "\t" + seq + "\t" + qual);
            }
            Console.WriteLine(ASETools.ElapsedTimeInSeconds(timer));

            outputFile.Close();

        }
    }
}
