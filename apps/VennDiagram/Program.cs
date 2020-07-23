using System;
using ASELib;
using System.IO;
using System.Diagnostics;
using System.Threading;
using System.Collections.Generic;
using System.Linq;
using System.Diagnostics.Contracts;

namespace VennDiagram
{
    class Program
    {
        static ASETools.Configuration configuration;
        static Dictionary<string, ASETools.Case> cases;
        static List<ASETools.Case> listOfCases;

        class ConcordancePair
        {
            public int both = 0;
            public int neither = 0;
            public int AButNotB = 0;
            public int BButNotA = 0;

            public void merge(ConcordancePair peer)
            {
                both += peer.both;
                neither += peer.neither;
                AButNotB += peer.AButNotB;
                BButNotA += peer.BButNotA;
            }
        } // ConcordancePair

        class ConcordanceResults
        {
            public ConcordanceResults()
            {
                foreach (var alignerPair in ASETools.alignerPairs)
                {
                    pairs.Add(alignerPair, new ConcordancePair());
                }
            }

            public readonly Dictionary<ASETools.AlignerPair, ConcordancePair> pairs = new Dictionary<ASETools.AlignerPair, ConcordancePair>();

            public void merge(ConcordanceResults peer)
            {
                foreach (var alignerPair in ASETools.alignerPairs)
                {
                    pairs[alignerPair].merge(peer.pairs[alignerPair]);
                }
            }
        } // ConcordanceResults

        static ConcordanceResults globalResults = new ConcordanceResults();

        static void Main(string[] args)
        {
            Stopwatch timer = new Stopwatch();
            timer.Start();

            configuration = ASETools.Configuration.loadFromFile(args);
            if (null == configuration)
            {
                return;
            }

            cases = ASETools.Case.LoadCases(configuration.casesFilePathname);
            if (null == cases)
            {
                Console.WriteLine("Unable to load cases");
                return;
            }

            listOfCases = cases.Select(_ => _.Value).ToList();

            foreach (var variantCaller in ASETools.EnumUtil.GetValues<ASETools.VariantCaller>())
            {
                foreach (var tumor in ASETools.BothBools)
                {
                    if (tumor) continue;  // For now, just ignore the tumor case.
                    var casesToProcess = listOfCases.Where(c => c.concordance.All(_ => _.Value[variantCaller][tumor].concordance_tarball_filename != "")).ToList();

                    Stopwatch oneRunTimer = new Stopwatch();
                    oneRunTimer.Start();

                    int nPerDot;
                    ASETools.PrintMessageAndNumberBar("Processing", "cases (" + ASETools.variantCallerName[variantCaller] + " " + ASETools.tumorToString[tumor] + ")", casesToProcess.Count(), out nPerDot);
                    var threading = new ASETools.WorkerThreadHelper<ASETools.Case, ConcordanceResults>(listOfCases, (c,v) => HandleOneCase(variantCaller, tumor, c, v), FinishUp, null, nPerDot);
                    Console.WriteLine(ASETools.ElapsedTimeInSeconds(oneRunTimer));
                } // tumor/normal
            } // variant caller
            
        } // Main

        static void FinishUp(ConcordanceResults state)
        {
            lock (globalResults)
            {
                globalResults.merge(state);
            } // lock (globalResults)
        } // FinishUp


        struct Locus
        {
            public readonly string contig;
            public readonly int pos;

            public Locus(string contig_, int pos_)
            {
                contig = contig_;
                pos = pos_;
            }
        }
        class Variant
        {
            public readonly Locus locus;
            public readonly string refAllele;
            public readonly string altAllele;
            public readonly bool fp, fn;

            Variant(Locus locus_, string refAllele_, string altAllele_, bool fp_, bool fn_)
            {
                locus = locus_;
                refAllele = refAllele_;
                altAllele = altAllele_;
                fp = fp_;
                fn = fn_;

                if (fp && fn)
                {
                    throw new Exception("Variant is both false positive and false negative");
                }
            }

            public static Variant FromVCFLine(string vcfLine)
            {
                var fields = vcfLine.Split('\t');
                if (fields.Count() < 11)
                {
                    throw new Exception("VCF line with too few fields: " + vcfLine);
                }

                if (fields[9].ToLower().Contains("fp") || fields[10].ToLower().Contains("fn"))
                {
                    throw new Exception("Variant with FP or FN on wrong allele: " + vcfLine);
                }

                return new Variant(new Locus(fields[0], Convert.ToInt32(fields[1])), fields[3], fields[4], fields[9].ToLower().Contains("fn"), fields[10].ToLower().Contains("fp"));
            } // FromVCFLine
        }

        static Dictionary<Locus, List<Variant>> LoadVariantsForOneVCF(ASETools.VariantCaller variantCaller, bool tumor, ASETools.Case case_, ASETools.AlignerPair alignerPair)
        {
            var retVal = new Dictionary<Locus, List<Variant>>();

            if (case_.concordance[alignerPair][variantCaller][tumor].concordance_tarball_filename == "")
            {
                return retVal;  // This is just here for debugging when we run without all the data
            }
            string tempDir = @"d:\temp\" + case_.case_id + @"_VennDiagram\";
            Directory.CreateDirectory(tempDir);

            string VCFName = case_.getDNAFileIdByTumor(tumor) + "." + alignerPair + ".concordance.vcf";
            string gzippedVCFName = VCFName + ".gz";

            ASETools.ExtractFileFromTarball(case_.concordance[alignerPair][variantCaller][tumor].concordance_tarball_filename, "./" + gzippedVCFName, tempDir);
            ASETools.RunAndWaitForProcess("gzip.exe", "-d " + gzippedVCFName);


            string inputLine;
            var inputFile = ASETools.CreateStreamReaderWithRetry(tempDir + VCFName);
            if (inputFile == null)
            {
                throw new Exception("Unable to open VCF file " + tempDir + VCFName);
            }
            while (null != (inputLine = inputFile.ReadLine()))
            {
                if (inputLine.StartsWith("#") || inputLine == "")
                {
                    continue;
                }

                var variant = Variant.FromVCFLine(inputLine);
                if (!retVal.ContainsKey(variant.locus))
                {
                    retVal.Add(variant.locus, new List<Variant>());
                }
                retVal[variant.locus].Add(variant);

            }

            inputFile.Close();
            Directory.Delete(tempDir, true);    // True recursively removes the directory and any contents

            return retVal;
        }

        static void HandleOneCase(ASETools.VariantCaller variantCaller, bool tumor, ASETools.Case case_, ConcordanceResults state)
        {
            var variants = new Dictionary<ASETools.AlignerPair, Dictionary<Locus, List<Variant>>>();
            var allLociWithVariants = new HashSet<Locus>();

            foreach (var alignerPair in ASETools.alignerPairs)
            {
                variants.Add(alignerPair, LoadVariantsForOneVCF(variantCaller, tumor, case_, alignerPair));

                foreach (var locus in variants[alignerPair].Select(_ => _.Key).ToList())
                {
                    allLociWithVariants.Add(locus);
                }
            }

            foreach (var alignerPair in ASETools.alignerPairs)
            {
                foreach (var locus in allLociWithVariants)
                {
                    if (!variants[alignerPair].ContainsKey(locus))
                    {
                        //
                        // Neither one has it.
                        //
                        state.pairs[alignerPair].neither++;
                    } 
                    else if (variants[alignerPair][locus].All(_ => _.fp))
                    {
                        state.pairs[alignerPair].BButNotA++;
                    }
                    else if (variants[alignerPair][locus].All(_ => _.fn))
                    {
                        state.pairs[alignerPair].AButNotB++;
                    } else
                    {
                        state.pairs[alignerPair].both++;
                    }
                } // foreach locus
            } // foreach aligner pair
        }
    } // Program
} // namespace VennDiagram
