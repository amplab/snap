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

        static ASETools.VariantCaller variantCaller;
        static bool tumor;

        static ASETools.BedFile lowCompexityRegions;


         static void Main(string[] args)
        {
            Stopwatch timer = new Stopwatch();
            timer.Start();

            configuration = ASETools.Configuration.loadFromFile(args);
            if (null == configuration)
            {
                return;
            }

            if (configuration.commandLineArgs.Count() < 4)
            {
                Console.WriteLine("usage: Venn variantCaller tumorOrNormal lowComplexityBEDFileName {caseID}");
                return;
            }

            if (!ASETools.tumorToString.ContainsValue(configuration.commandLineArgs[1]))
            {
                Console.WriteLine("Second parameter must be either " + ASETools.tumorToString[true] + " or " + ASETools.tumorToString[false]);
                return;
            }

            if (!ASETools.variantCallerName.ContainsValue(configuration.commandLineArgs[0]))
            {
                Console.Write("First parameter must be a variant caller name, one of:");
                ASETools.alignerName.Select(_ => _.Value).ToList().ForEach(_ => Console.Write(" " + _));
                Console.WriteLine();
                return;
            }

            tumor = configuration.commandLineArgs[1] == ASETools.tumorToString[true];
            variantCaller = ASETools.variantCallerName.Where(_ => _.Value == configuration.commandLineArgs[0]).ToList()[0].Key;

            lowCompexityRegions = new ASETools.BedFile(args[2]);

            cases = ASETools.Case.LoadCases(configuration.casesFilePathname);
            if (null == cases)
            {
                Console.WriteLine("Unable to load cases");
                return;
            }

            for (int i = 3; i < configuration.commandLineArgs.Count(); i++)
            {
                if (!cases.ContainsKey(configuration.commandLineArgs[i]))
                {
                    Console.WriteLine(configuration.commandLineArgs[i] + " isn't a case ID");
                    return;
                }
            }

            listOfCases = cases.Select(_ => _.Value).ToList();

            var casesToRun = listOfCases.Where(_ => configuration.commandLineArgs.ToList().Contains(_.case_id)).ToList();

            if (casesToRun.Any(case_ => case_.concordance.Any(alignerPair => case_.concordance[alignerPair.Key][variantCaller][tumor].concordance_tarball_filename == ""))) 
            {
                Console.WriteLine("One or more cases is missing a concordance tarball for this variant caller and tumor/normal pair.  Skipping them.");
                casesToRun = casesToRun.Where(case_ => case_.concordance.All(alignerPair => case_.concordance[alignerPair.Key][variantCaller][tumor].concordance_tarball_filename != "")).ToList();
            }

            if (casesToRun.Any(case_ => case_.concordance.Any(alignerPair => case_.concordance[alignerPair.Key][variantCaller][tumor].concordance_tarball_size < 512 * 1024)))
            {
                Console.Write("The following concordance tarballs are suspiciously small and will not be processed:");
                foreach (var caseWithSmallTarball in casesToRun.Where(case_ => case_.concordance.Any(alignerPair => case_.concordance[alignerPair.Key][variantCaller][tumor].concordance_tarball_size < 512 * 1024)))
                {
                    foreach (var alignerPair in ASETools.allAlignerPairs)
                    {
                        if (caseWithSmallTarball.concordance[alignerPair][variantCaller][tumor].concordance_tarball_size < 512 * 1024) 
                        {
                            Console.Write(" " + caseWithSmallTarball.concordance[alignerPair][variantCaller][tumor].concordance_tarball_filename + " (" + caseWithSmallTarball.concordance[alignerPair][variantCaller][tumor].concordance_tarball_size + ")");
                        } // bad tarball
                    } // aligner pair
                } // bad case

                Console.WriteLine();
                casesToRun = casesToRun.Where(case_ => case_.concordance.All(alignerPair => case_.concordance[alignerPair.Key][variantCaller][tumor].concordance_tarball_size >= 512 * 1024)).ToList();
            } // if any tarballs are too small

            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "cases", casesToRun.Count(), out nPerDot);
            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, int>(casesToRun, HandleOneCase, null, null, nPerDot);
            threading.run();
            Console.WriteLine(ASETools.ElapsedTimeInSeconds(timer));
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

            public ASETools.VariantType getVariantType()
            {
                return isIndel() ? ASETools.VariantType.Indel : ASETools.VariantType.SNV;
            }

            public static Variant FromVCFLine(string vcfLine, ASETools.AlignerPair alignerPair)
            {
                var fields = vcfLine.Split('\t');
                if (fields.Count() < 11)
                {
                    throw new Exception("VCF line with too few fields: " + vcfLine);
                }

#if false // too slow
                if (fields[9].Contains("FP") || fields[10].Contains("FN"))
                {
                    throw new Exception("Variant with FP or FN on wrong allele: " + vcfLine);
                }
#endif


                return new Variant(new Locus(fields[0], Convert.ToInt32(fields[1])), fields[3], fields[4], fields[9].Contains("FN"), fields[10].Contains("FP"), alignerPair);
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

                if (lowCompexityRegions.isRangeIncluded(variant.locus.contig, variant.locus.pos))
                {
                    // Skip the ones in low-complexity regions.
Console.WriteLine("Skipping variant in LCR " + variant.locus.contig + ":" + variant.locus.pos);
                    continue;
                }
Console.WriteLine("Keeping variant " + variant.locus.contig + ":" + variant.locus.pos);

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

        static void HandleOneCase(ASETools.Case case_, int state)
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
                    }

                    variantsByLocus[locus].AddRange(variantsForThisPair[locus]);
                } // locus
            } // aligner pair

            var variantCountByAlignerSetAndType = new Dictionary<ASETools.AlignerSet, Dictionary<ASETools.VariantType, int>>();
            foreach (var alignerSet in ASETools.allAlignerSets)
            {
                variantCountByAlignerSetAndType.Add(alignerSet, new Dictionary<ASETools.VariantType, int>());
                foreach (var variantType in ASETools.EnumUtil.GetValues<ASETools.VariantType>())
                {
                    variantCountByAlignerSetAndType[alignerSet].Add(variantType, 0);
                }
            }

            foreach (var locus in variantsByLocus.Select(_ => _.Key))
            {
                var variantsToProcess = variantsByLocus[locus].Where(_ => true).ToList();    // just a funky way to copy the list

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

            var outputFilename = ASETools.GetDirectoryFromPathname(case_.concordance[ASETools.allAlignerPairs[0]][variantCaller][tumor].concordance_tarball_filename) + @"\" +
                case_.getDNAFileIdByTumor(tumor) + "." + ASETools.tumorToString[tumor] + "-" + ASETools.variantCallerName[variantCaller] + "-Venn.txt";
            var outputFile = ASETools.CreateStreamWriterWithRetry(outputFilename);

            if (outputFile == null)
            {
                Console.WriteLine("Unable to open output file " + outputFilename);
                return;
            }

            bool anyFieldsWritten = false;
            foreach (var variantType in ASETools.EnumUtil.GetValues<ASETools.VariantType>())
            {
                if (anyFieldsWritten)
                {
                    outputFile.Write("\t");
                } else
                {
                    anyFieldsWritten = true;
                }

                outputFile.Write(" nVariants" + " " + ASETools.variantTypeToName[variantType]);

                foreach (var alignerSet in ASETools.allAlignerSets)
                {
                    outputFile.Write("\t" + alignerSet.ToString() + " " + ASETools.variantTypeToName[variantType]);
                }
            }

            outputFile.WriteLine();

            anyFieldsWritten = false;
            foreach (var variantType in ASETools.EnumUtil.GetValues<ASETools.VariantType>())
            {
                if (anyFieldsWritten)
                {
                    outputFile.Write("\t");
                }
                else
                {
                    anyFieldsWritten = true;
                }

                outputFile.Write(variantCountByAlignerSetAndType.Sum(_ => _.Value[variantType]));

                foreach (var alignerSet in ASETools.allAlignerSets)
                {
                    outputFile.Write("\t" + variantCountByAlignerSetAndType[alignerSet][variantType]);
                }
            }

            outputFile.WriteLine();
            outputFile.WriteLine("**done**");
            outputFile.Close();
        } // HandleOneCase


    } // Program

}
