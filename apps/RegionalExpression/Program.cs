using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.IO.Compression;
using System.Diagnostics;
using System.Threading;
using ASELib;


namespace RegionalExpression
{
	class Program
	{

	

        static int regionSize;

        class OneRun
        {
            public string allcountFilename;
            public string tumor_rna_file_id;
            public string case_id;
        }

        static void ProcessRuns(List<OneRun> runs, ASETools.ExpressionFile expression)
        {
            while (true) {

                string allcountFilename;
                string tumor_rna_file_id;
                string case_id;
                lock (runs)
                {
                    if (runs.Count() == 0)
                    {
                        //
                        // No more work, we're done.
                        //
                        return;
                    }

                    allcountFilename = runs[0].allcountFilename;
                    tumor_rna_file_id = runs[0].tumor_rna_file_id;
                    case_id = runs[0].case_id;

                    runs.RemoveAt(0);
                }

                //
                // Run through the allcount file
                //
                var allcountTimer = new Stopwatch();
                allcountTimer.Start();

                var allcountReader = new ASETools.AllcountReader(allcountFilename);
                long mappedHQNuclearReads;
                int numContigs;

                bool fileCorrupt = !allcountReader.openFile(out mappedHQNuclearReads, out numContigs);

                if (fileCorrupt)
                {
                    Console.WriteLine("Ignoring corrupt allcount file " + allcountFilename);
                    continue;
                }

                //
                // Now process the body of the file.
                //

                int indexOfLastSlash = allcountFilename.LastIndexOf('\\');
                if (-1 == indexOfLastSlash) {
                    Console.WriteLine("Couldn't find a backslash in allcount pathname, which is supposed to be absolute: " + allcountFilename);
                    continue;
                }

                string directory = ASETools.GetDirectoryFromPathname(allcountFilename);
                var outputFilename = directory + @"\" + tumor_rna_file_id + ASETools.regionalExpressionExtension;
                var outputFile = new StreamWriter(outputFilename);

                outputFile.WriteLine("RegionalExpression v3.0\t" + tumor_rna_file_id + "\t" + allcountFilename + "\t" + regionSize);
                outputFile.WriteLine("NumContigs: " + numContigs);
                ASETools.Region.printHeader(outputFile);

				ASETools.Region region = new ASETools.Region(regionSize, expression, expression.higestOffsetForEachContig, mappedHQNuclearReads, outputFile);

                ASETools.AllcountReader.ProcessBase processBase = (x, y, z) => region.processBase(x, y, z);

                fileCorrupt = !allcountReader.ReadAllcountFile(processBase);

                region.closeRegion();

                if (fileCorrupt)
                {
                    outputFile.Close();
                    File.Delete(outputFilename);
                }
                else
                {
                    outputFile.WriteLine("**done**");
                    outputFile.Close();

                    allcountTimer.Stop();
                    lock (runs)
                    {
                        Console.WriteLine("Processed " + case_id + " in " + (allcountTimer.ElapsedMilliseconds + 500) / 1000 + "s, " + runs.Count() + " remain" +  ((runs.Count() == 1) ? "s" : "") + " in the queue.");
                    }
                }
            }
        }

        static void Main(string[] args)
        {
            var configuration = ASETools.Configuration.loadFromFile(args);

            if (configuration.commandLineArgs.Count() < 3 || configuration.commandLineArgs[2] == "-f" && configuration.commandLineArgs.Count() != 4)
            {
                Console.WriteLine("usage: RegionalExpression expression_disease_file regionSize <one or more caseIDs with the same disease | -f inputFilename> ");
                return;
            }

            try
            {
                regionSize = Convert.ToInt32(configuration.commandLineArgs[1]);
            }
            catch (FormatException)
            {
                Console.WriteLine("Couldn't parse region size from command line");
                return;
            }

            if (regionSize < 1)
            {
                Console.WriteLine("Bogus regionSize from command line");
                return;
            }

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

            ASETools.ExpressionFile expression;
            

            Stopwatch timer;
            try
            {
                Console.Write("Loading expression file " + configuration.commandLineArgs[0] + "...");
                timer = new Stopwatch();
                timer.Start();

                expression = new ASETools.ExpressionFile();

                expression.LoadFromFile(configuration.commandLineArgs[0]);   // This can take forever!

                timer.Stop();
                Console.WriteLine((timer.ElapsedMilliseconds + 500) / 1000 + "s");
            }
            catch (OutOfMemoryException)
            {
                Console.WriteLine("Out of memory exception loading the expression file (I'm really not sure why this happens when there's plenty of memory).  Running in the visual studio debugger seems to help, though.");
                return;
            }

            var runs = new List<OneRun>();

            List<string> case_ids;

            if (configuration.commandLineArgs[2] == "-f")
            {
                case_ids = File.ReadAllLines(configuration.commandLineArgs[3]).ToList();
            }
            else
            {
                case_ids = new List<string>();
                for (int i = 2; i < configuration.commandLineArgs.Count(); i++)
                {
                    case_ids.Add(configuration.commandLineArgs[i]);
                }
            }
           
            foreach (var case_id in case_ids)  // for each person we're processing
            {
                if (!cases.ContainsKey(case_id))
                {
                    Console.WriteLine("Couldn't find case_id " + case_id + ", are you sure it's a correct case ID? Ignoring");
                    continue;
                }

                var case_ = cases[case_id];

                string allcountFilename = case_.tumor_rna_allcount_filename;
                if (allcountFilename == "")
                {
                    Console.WriteLine("Case " + case_id + " doesn't appear to have a tumor RNA allcount file");
                    continue;
                }

                var run = new OneRun();
                run.allcountFilename = allcountFilename;
                run.tumor_rna_file_id = case_.tumor_rna_file_id;
                run.case_id = case_id;

                runs.Add(run);
            }

            //
            // Process the runs in parallel
            //
            timer.Reset();
            timer.Start();

            var threads = new List<Thread>();
            for (int i = 0; i < Environment.ProcessorCount; i++) {
                threads.Add(new Thread(() => ProcessRuns(runs, expression)));
            }

            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());

            timer.Stop();
            Console.WriteLine("Processed " + (configuration.commandLineArgs.Count() - 2) + " cases in " + ASETools.ElapsedTimeInSeconds(timer));
         }
    }
}
