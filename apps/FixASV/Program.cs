using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;
using System.Threading;
using System.Diagnostics;

namespace FixASV
{
    class Program
    {
        //
        // One version of AnnotateVariants didn't write the trailing tabs in lines for which there was no normal RNA samples.  This causes HeaderizedFile no end of greif, so rather than rerunning
        // AnnotateVariants (it's already fixed), this program looks for short lines and just adds in the trailing tabs.
        //
        static List<ASETools.Case> casesToProcess = new List<ASETools.Case>();
        static int TotalProcessed = 0;


        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            var configuration = ASETools.ASEConfirguation.loadFromFile(args);

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

            foreach (var caseEntry in cases)
            {
                var case_ = caseEntry.Value;

                if (case_.annotated_selected_variants_filename != "")
                {
                    casesToProcess.Add(case_);
                }
            }

            Console.Write("Processing " + casesToProcess.Count() + " cases (1 dot/100 cases): ");

            var threads = new List<Thread>();
            for (int i = 0; i < Environment.ProcessorCount; i++)
            {
                threads.Add(new Thread(() => WorkerThread()));
            }

            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());

            Console.WriteLine();
            Console.Write("Total runtime " + ASETools.ElapsedTimeInSeconds(timer));

        } // Main

        static void WorkerThread()
        {
            ASETools.Case case_;

            lock (casesToProcess)
            {
                if (casesToProcess.Count() == 0)
                {
                    return;
                }

                case_ = casesToProcess[0];
                casesToProcess.RemoveAt(0);
            } // lock

            var inputLines = ASETools.ReadAllLinesWithRetry(case_.annotated_selected_variants_filename);

            var outputFile = ASETools.CreateStreamWriterWithRetry(case_.annotated_selected_variants_filename);

            for (int i = 0; i < inputLines.Count(); i++)
            {
                if (inputLines[i] == "**done**")
                {
                    outputFile.WriteLine("**done**");
                    continue;
                }

                var fields = inputLines[i].Split('\t');
                if (fields.Count() == 24)
                {
                    outputFile.WriteLine(inputLines[i]);
                } else if (fields.Count() == 20)
                {
                    outputFile.WriteLine(inputLines[i] + "\t\t\t\t");
                } else
                {
                    Console.WriteLine("Strange input line in file " + case_.annotated_selected_variants_filename + ": " + inputLines[i] + ".  Field count " + fields.Count());
                    outputFile.WriteLine(inputLines[i]);
                }
            } // for all input lines
            outputFile.Close();

            if (Interlocked.Increment(ref TotalProcessed) % 100 == 0)
            {
                Console.Write(".");
            }
        }
    } // Program
} // FixASV
