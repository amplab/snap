using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using ASELib;
using System.Diagnostics;

namespace MeasureRNADepthDistributionAtMutantSites
{
    class Program
    {
        static ASETools.Configuration configuration;
        static List<ASETools.Case> cases;

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            configuration = ASETools.Configuration.loadFromFile(args);

            if (null == configuration)
            {
                Console.WriteLine("Giving up because we were unable to load configuration.");
                return;
            }

            if (configuration.commandLineArgs.Count() == 0)
            {
                Console.WriteLine("usage: MeasureRNAReadDepthDistributonAtMutantSites {OneOrMoreDiseaseIdentifiers}");
                return;
            }

            cases = ASETools.Case.LoadCases(configuration.casesFilePathname).Select(x => x.Value).ToList();
            if (null == cases)
            {
                Console.WriteLine("Unable to load cases from file " + configuration.casesFilePathname);
                return;
            }

            var nDiseases = args.Count();

            var threading = new ASETools.WorkerThreadHelper<string, int>(args.ToList(), HandleOneDisease, null, null, 0);

            threading.run();

            Console.WriteLine("Processed " + nDiseases + " in " + ASETools.ElapsedTimeInSeconds(timer));
        }

        static void HandleOneDisease(string disease, int unused)
        {
            var ReadDepths = new Dictionary<string, Dictionary<int, List<double>>>();

            var interestingCases = cases.Where(x => x.disease() == disease).ToList();

            var nCases = interestingCases.Count();

            if (nCases == 0)
            {
                Console.WriteLine("Found no cases for disease " + disease + ".  Is it really a disease identifier?");
                return;
            }

            if (interestingCases[0].maf_filename == "")
            {
                Console.WriteLine("Must download the MAF file before you can run this.");
                return;
            }

            if (interestingCases.Any(x => x.maf_filename != interestingCases[0].maf_filename))
            {
                Console.WriteLine("Huh?  Not all cases for disase " + disease + " have the same MAF file name.  Skipping disease " + disease + ".");
                return;
            }

            if (interestingCases.Any(x => x.tumor_rna_allcount_filename == ""))
            {
                Console.WriteLine("At least one case doesn't have a tumor RNA allcount file.  Skipping disease " + disease + ".");
                return;
            }

            var mafLines = ASETools.MAFLine.ReadFile(interestingCases[0].maf_filename, interestingCases[0].maf_file_id, true, configuration.isBeatAML);

            int nLoci = 0;

            foreach (var mafLine in mafLines)
            {
                var chr = mafLine.Chromosome.ToLower();
                if (!ReadDepths.ContainsKey(chr))
                {
                    ReadDepths.Add(chr, new Dictionary<int, List<double>>());
                }

                if (!ReadDepths[chr].ContainsKey(mafLine.Start_Position)) {
                    ReadDepths[chr].Add(mafLine.Start_Position, new List<double>());
                    nLoci++;
                }
            }

            Console.Write("Disease " + disease + " has " + nLoci + " mutant loci.  Processing " + nCases + " cases, 1 dot/10 cases: ");

            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, int>(interestingCases, (x, y) => ProcessOneCase(x, ReadDepths), null, null, 10); // Yes, n^2 threads.  They'll all be blocked on IO anyway.  Plus, we only run one disease/machine at a time anyway
            threading.run();
            Console.WriteLine();

            var outputFile = ASETools.CreateStreamWriterWithRetry(ASETools.Configuration.expression_distribution_directory + ASETools.Expression_distribution_filename_base + disease);
            outputFile.Write("Chromosome\tLocus\tmin");
            for (int i = 10; i < 100; i += 10)
            {
                outputFile.Write("\t" + i + "th %ile");
            }
            outputFile.WriteLine("\tmax");

            foreach (var chromEntry in ReadDepths)
            {
                var chromosome = chromEntry.Key;

                foreach (var locusEntry in chromEntry.Value)
                {
                    var locus = locusEntry.Key;

                    var values = locusEntry.Value;

                    int nValues = values.Count();

                    for (int i = nValues; i < nCases; i++)
                    {
                        values.Add(0);  // No reads at all in the allcount file don't get reported.  Add them here.
                    }
                    values.Sort();

                    outputFile.Write(chromosome + "\t" + locus);

                    outputFile.Write("\t" + values[0]);
                    for (int i = 10; i < 100; i+=10)
                    {
                        outputFile.Write("\t" + values[(int)Math.Round(((double)i / 100) * (nCases - 1), 1)]);
                    }
                    outputFile.WriteLine("\t" + values[nCases - 1]);

                } // For each mutant locus 
            } // for each contig

            outputFile.WriteLine("**done**");
            outputFile.Close();
        } // HandleOneDisease

        static void ProcessOneCase(ASETools.Case case_, Dictionary<string, Dictionary<int, List<double>>> ReadDepths)
        {
            var allcountReader = new ASETools.AllcountReader(case_.tumor_rna_allcount_filename);

            long numMappedHQNuclearReads;
            int numContigs;

            var numMappedHQNuclearBases = ASETools.MappedBaseCount.readFromFile(case_.tumor_rna_mapped_base_count_filename).mappedBaseCount;

            if (!allcountReader.openFile(out numMappedHQNuclearReads, out numContigs))
            {
                Console.WriteLine("Failed to open allcount file " + case_.tumor_rna_allcount_filename);
                return;
            }

            allcountReader.ReadAllcountFile((x, y, z) => ProcessOneBase(ReadDepths, numMappedHQNuclearBases, x, y, z));

            allcountReader.Close();
        }

        static void ProcessOneBase(Dictionary<string, Dictionary<int, List<double>>> ReadDepths, long numMappedHQNuclearBases, string contig, int locus, int mappedReadCount)
        {
            contig = contig.ToLower();
            if (!ReadDepths.ContainsKey(contig) || !ReadDepths[contig].ContainsKey(locus))
            {
                return;
            }

            lock (ReadDepths[contig][locus])
            {
                ReadDepths[contig][locus].Add((double)mappedReadCount / numMappedHQNuclearBases);
            }
        }


    }
}
