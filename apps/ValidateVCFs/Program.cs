using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Diagnostics;
using System.Threading;
using ASELib;

namespace ValidateVCFs
{
    class Program
    {

        static Stopwatch timer = new Stopwatch();
        static List<KeyValuePair<string, ASETools.Case>> casesToVerify;
        static int nFailed = 0;

        static void usage()
        {
            Console.WriteLine("usage: ValidateVCFs -a | <caseIDs>");
            Console.WriteLine("-a means to validate them all.");
        }

        static void ValidateVCFs()
        {
            while (true)
            {
                ASETools.Case case_;

                lock (casesToVerify)
                {
                    if (casesToVerify.Count() == 0)
                    {
                        return;
                    }

                    case_ = casesToVerify[0].Value;
                    casesToVerify.RemoveAt(0);
                } // lock 

                var vcfFile = ASETools.CreateStreamReaderWithRetry(case_.vcf_filename);

                string line;

                //
                // Skip the header.
                //
                while (null != (line = vcfFile.ReadLine()) && line.Count() != 0 && line[0] == '#')
                {
                    // Skip the header lines
                }

                if (null == line || line.Count() == 0)
                {
                    Console.WriteLine("Corrupt vcf: missing body: " + case_.vcf_filename);
                    continue;
                }

                var chromosomesLeftToSee = new List<string>();
                for (int i = 1; i <= 22; i++)
                {
                    chromosomesLeftToSee.Add("chr" + i);
                }

                chromosomesLeftToSee.Add("chrX");
                //
                // Don't do Y, because of women.  :-)
                //

                for (; line != null; line = vcfFile.ReadLine())
                {
                    var fields = line.Split('\t');

                    if (fields.Count() != 10)
                    {
                        Console.WriteLine(case_.vcf_filename + " has the wrong number of fields (" + fields.Count() + " != 10 on line: '" + line + "'.");
                        chromosomesLeftToSee = new List<string>(); // Just to suppress the error
                        break;
                    }

                    if (chromosomesLeftToSee.Contains(fields[0]))
                    {
                        chromosomesLeftToSee.Remove(fields[0]);
                    }
                }

                if (chromosomesLeftToSee.Count() != 0)
                {
                    lock (casesToVerify)
                    {
                        nFailed++;

                        Console.Write(case_.vcf_filename + " has no variants called for chromosome(s):");
                        foreach (var chromosome in chromosomesLeftToSee)
                        {
                            Console.Write(" " + chromosome);
                        }
                        Console.WriteLine();
                    }
                }

            } // while (true)
        } // ValidateVCFs

        static void ProgressThread()
        {
            long lastPrintAt = 0;

            while (true)
            {
                var now = timer.ElapsedMilliseconds;
                if (now < lastPrintAt + 60000)
                {
                    Thread.Sleep((int)(lastPrintAt + 60000 - now));
                }

                lastPrintAt = timer.ElapsedMilliseconds;

                lock (casesToVerify)
                {
                    if (casesToVerify.Count() == 0)
                    {
                        return;
                    }
                    Console.WriteLine(ASETools.ElapsedTimeInSeconds(timer) + ": " + casesToVerify.Count() + " remain queued.");
                }
            }
        }

        static void Main(string[] args)
        {
            timer.Start();

            var configuration = ASETools.ASEConfirguation.loadFromFile(args);

            if (configuration.commandLineArgs.Count() == 0)
            {
                usage();
                return;
            }

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

            if (configuration.commandLineArgs.Count() == 1 && configuration.commandLineArgs[0] == "-a")
            {
                casesToVerify = cases.Where(x => x.Value.vcf_filename != null && x.Value.vcf_filename != "").ToList();
            }
            else
            {
                casesToVerify = cases.Where(x => x.Value.vcf_filename != null && x.Value.vcf_filename != "" && configuration.commandLineArgs.Contains(x.Value.case_id)).ToList();
                if (casesToVerify.Count() != configuration.commandLineArgs.Count()) {
                    Console.WriteLine("Some args don't correspond to case IDs that have VCFs.");
                    usage();
                    return;
                }
            }

            var nToVerify = casesToVerify.Count();
            Console.WriteLine("Validating " + nToVerify + " VCFs.");

            var threads = new List<Thread>();
            for (int i = 0; i < Environment.ProcessorCount; i++)
            {
                threads.Add(new Thread(() => ValidateVCFs()));
            }

            threads.ForEach(t => t.Start());

            new Thread(() => ProgressThread()).Start();

            threads.ForEach(t => t.Join());

            Console.WriteLine("Validated " + nToVerify + " VCFs (" + nFailed + " failed) in " + ASETools.ElapsedTimeInSeconds(timer));
        } // Main
    }
}
