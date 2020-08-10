using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.Diagnostics;
using System.IO;
using System.Threading;
using System.Runtime.CompilerServices;
using System.CodeDom;

namespace Venn
{
    class Program
    {
        static ASETools.Configuration configuration;
        static Dictionary<string, ASETools.Case> cases;
        static List<ASETools.Case> listOfCases;

        enum VariantType { SNV, Indel };

        class ConcordanceResults
        {
            public ConcordanceResults()
            {
                foreach (var alignerSet in ASETools.allAlignerSets)
                {
                    count.Add(alignerSet, new Dictionary<VariantType, ASETools.RunningMeanAndStdDev>());
                    foreach (var variantType in ASETools.EnumUtil.GetValues<VariantType>())
                    {
                        count[alignerSet].Add(variantType, new ASETools.RunningMeanAndStdDev());
                    }
                }
            } // ctor

            public readonly Dictionary<ASETools.AlignerSet, Dictionary<VariantType, ASETools.RunningMeanAndStdDev>> count = 
                new Dictionary<ASETools.AlignerSet, Dictionary<VariantType, ASETools.RunningMeanAndStdDev>>();
            public void merge(ConcordanceResults peer)
            {
                foreach (var alignerSet in ASETools.allAlignerSets)
                {
                    foreach (var variantType in ASETools.EnumUtil.GetValues<VariantType>())
                    {
                        count[alignerSet][variantType].merge(peer.count[alignerSet][variantType]);
                    }
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
                    int nCases = casesToProcess.Count();

                    Stopwatch oneRunTimer = new Stopwatch();
                    oneRunTimer.Start();

                    int nPerDot;
                    ASETools.PrintMessageAndNumberBar("Processing", "cases (" + ASETools.variantCallerName[variantCaller] + " " + ASETools.tumorToString[tumor] + ")", casesToProcess.Count(), out nPerDot);
                    var threading = new ASETools.WorkerThreadHelper<ASETools.Case, ConcordanceResults>(casesToProcess, (c, v) => HandleOneCase(variantCaller, tumor, c, v), FinishUp, null, nPerDot);
                    threading.run(1 /*BJB*/);
                    Console.WriteLine(ASETools.ElapsedTimeInSeconds(oneRunTimer));


                    foreach (var variantType in ASETools.EnumUtil.GetValues<VariantType>()) {
                        Console.WriteLine(ASETools.tumorToString[tumor] + " " + ASETools.variantCallerName[variantCaller] + " " +(variantType == VariantType.SNV ? "SNVs" : "Indels"));
                        foreach (var alignerSet in ASETools.allAlignerSets)
                        {
                            var meanAndStdDev = globalResults.count[alignerSet][variantType].getMeanAndStdDev();
                            Console.WriteLine(alignerSet + "\t" + meanAndStdDev.mean  + "\t" + meanAndStdDev.stddev);
                        }

                        Console.WriteLine();
                        Console.WriteLine("Total area");
                        foreach (var alignerSet in ASETools.allAlignerSets)
                        {
                            var n = globalResults.count.Where(_ => alignerSet.isSubsetOf(_.Key)).Sum(_ => _.Value[variantType].getMeanAndStdDev().mean);
                            Console.WriteLine(alignerSet + " = " + n);
                        }

                        Console.WriteLine();
                    }
                    Console.WriteLine();
                } // tumor/normal
            } // variant caller

        } // Main


        struct Locus
        {
            public readonly string contig;
            public readonly int pos;

            public Locus(string contig_, int pos_)
            {
                contig = contig_;
                pos = pos_;
            }

            public override int GetHashCode()
            {
                return contig.GetHashCode() ^ pos;
            }

            static public bool operator==(Locus us, Locus them)
            {
                return us.pos == them.pos && us.contig == them.contig;
            }

            static public bool operator!=(Locus us, Locus them)
            {
                return !(us == them);
            }

            public override bool Equals(object o)
            {
                return this == (Locus)o;
            }
        }
        class Variant
        {
            public readonly Locus locus;
            public readonly string refAllele;
            public readonly string altAllele;
            public readonly bool fp, fn;
            public readonly ASETools.AlignerPair alignerPair;

            Variant(Locus locus_, string refAllele_, string altAllele_, bool fp_, bool fn_, ASETools.AlignerPair alignerPair_)
            {
                locus = locus_;
                refAllele = refAllele_;
                altAllele = altAllele_;
                fp = fp_;
                fn = fn_;
                alignerPair = alignerPair_;
            }

            public bool isIndel()
            {
                return refAllele.Length != 1 || altAllele.Length != 1;
            }

            public VariantType getVariantType()
            {
                return isIndel() ? VariantType.Indel : VariantType.SNV;
            }

            public static Variant FromVCFLine(string vcfLine, ASETools.AlignerPair alignerPair)
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

                return new Variant(new Locus(fields[0], Convert.ToInt32(fields[1])), fields[3], fields[4], fields[9].ToLower().Contains("fn"), fields[10].ToLower().Contains("fp"), alignerPair);
            } // FromVCFLine

            public bool presentInAligner(ASETools.Aligner aligner)
            {
                if (alignerPair.firstAligner == aligner)
                {
                    return !fn;
                } else if (alignerPair.secondAligner == aligner)
                {
                    return !fp;
                } else
                {
                    return false;
                }
            } // presentInAligner

            public bool equivalentTo(Variant peer)
            {
                return locus == peer.locus && refAllele == peer.refAllele && altAllele == peer.altAllele;
            }
        } // Variant

        static Dictionary<Locus, List<Variant>> LoadVariantsForOneVCF(ASETools.VariantCaller variantCaller, bool tumor, ASETools.Case case_, ASETools.AlignerPair alignerPair)
        {
            var retVal = new Dictionary<Locus, List<Variant>>();

            if (case_.concordance[alignerPair][variantCaller][tumor].concordance_tarball_filename == "")
            {
                return retVal;  // This is just here for debugging when we run without all the data
            }
            string tempDir = @"d:\temp\" + case_.case_id + @"_VennDiagram\";
            ASETools.TryToRecursivelyDeleteDirectory(tempDir);

            Directory.CreateDirectory(tempDir);

            string VCFName = case_.getDNAFileIdByTumor(tumor) + "." + alignerPair + ".concordance.vcf";
            string gzippedVCFName = VCFName + ".gz";

            ASETools.ExtractFileFromTarball(case_.concordance[alignerPair][variantCaller][tumor].concordance_tarball_filename, "./" + gzippedVCFName, tempDir);

            ASETools.RunAndWaitForProcess(@"c:\bolosky\bin\gzip.exe", "-d " + tempDir + gzippedVCFName);

            var inputFile = ASETools.CreateStreamReaderWithRetry(tempDir + VCFName);
            if (inputFile == null)
            {
                throw new Exception("Couldn't open unzipped input file " + tempDir + VCFName);
            }
            string inputLine;
            while (null != (inputLine = inputFile.ReadLine()))
            {
                if (inputLine.StartsWith("#") || inputLine == "")
                {
                    continue;
                }

                var variant = Variant.FromVCFLine(inputLine, alignerPair);
                if (!retVal.ContainsKey(variant.locus))
                {
                    retVal.Add(variant.locus, new List<Variant>());
                }

                retVal[variant.locus].Add(variant);

            }

            inputFile.Close();

            ASETools.TryToRecursivelyDeleteDirectory(tempDir);

            return retVal;
        }

        static void HandleOneCase(ASETools.VariantCaller variantCaller, bool tumor, ASETools.Case case_, ConcordanceResults state)
        {
            var variantsByLocus = new Dictionary<Locus, List<Variant>>();

            foreach (var alignerPair in ASETools.allAlignerPairs)
            {
                var variantsForThisPair = LoadVariantsForOneVCF(variantCaller, tumor, case_, alignerPair);

                foreach (var locus in variantsForThisPair.Select(_ => _.Key))
                {
                    if (!variantsByLocus.ContainsKey(locus))
                    {
                        variantsByLocus.Add(locus, new List<Variant>());
                    } else { }

                    variantsByLocus[locus].AddRange(variantsForThisPair[locus]);
                } // locus
            } // aligner pair

            var variantCountByAlignerSetAndType = new Dictionary<ASETools.AlignerSet, Dictionary<VariantType, int>>();
            foreach (var alignerSet in ASETools.allAlignerSets)
            {
                variantCountByAlignerSetAndType.Add(alignerSet, new Dictionary<VariantType, int>());
                foreach (var variantType in ASETools.EnumUtil.GetValues<VariantType>())
                {
                    variantCountByAlignerSetAndType[alignerSet].Add(variantType, 0);
                }
            }

            foreach (var locus in variantsByLocus.Select(_ => _.Key))
            {
                var variantsToProcess = variantsByLocus[locus].Where(_ => true).ToList();    // just a funky to copy the list

                while (variantsToProcess.Count() > 0)
                {
                    var equivalentVariants = variantsToProcess.Where(_ => _.equivalentTo(variantsToProcess[0])).ToList();

                    var presentInAligners = new ASETools.AlignerSet(ASETools.EnumUtil.GetValues<ASETools.Aligner>().Where(aligner => equivalentVariants.Any(variant => variant.presentInAligner(aligner))).ToList());

                    var variantType = equivalentVariants[0].getVariantType();

                    variantCountByAlignerSetAndType[presentInAligners][variantType]++;

                    if (equivalentVariants.Count() >= variantsToProcess.Count())
                    {
                        break;
                    }
                    equivalentVariants.ForEach(_ => variantsToProcess.Remove(_));
                }
            } // for each locus with variants

            foreach (var variantType in ASETools.EnumUtil.GetValues<VariantType>())
            {
                int nVariants = variantCountByAlignerSetAndType.Sum(_ => _.Value[variantType]);

                foreach (var alignerSet in ASETools.allAlignerSets)
                {
                    state.count[alignerSet][variantType] += ((double)variantCountByAlignerSetAndType[alignerSet][variantType]) / nVariants;
                } // alignerSet
            } // variantType
        } // HandleOneCase
        static void FinishUp(ConcordanceResults state)
        {
            lock (globalResults)
            {
                globalResults.merge(state);
            } // lock (globalResults)
        } // FinishUp

    } // Program

}
