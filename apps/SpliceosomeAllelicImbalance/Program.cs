using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Diagnostics;
using ASELib;

namespace SpliceosomeAllelicImbalance
{
    class Program
    {
        static ASETools.CommonData commonData;
        static ASETools.IsoformMap isoformMap;
        // hugo_symbol->isoform->spliceosomemutant->frac of all reads mapped to this gene that are in this isoform
        static Dictionary<string, Dictionary<string, Dictionary<bool, List<double>>>> isoformBalanceBySpliceosome = new Dictionary<string, Dictionary<string, Dictionary<bool, List<double>>>>();   


        class DoubleAndBool : IComparer<DoubleAndBool>
        {
            public bool whichClass;
            public double value;

            public DoubleAndBool(bool whichClass_, double value_)
            {
                whichClass = whichClass_;
                value = value_;
            }

            public int Compare(DoubleAndBool x, DoubleAndBool y)
            {
                return x.value.CompareTo(y.value);
            }
        }

        class NameAndPValue
        {
            public string name;
            public double pValue;
        }

        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);

            if (null == commonData)
            {
                return;
            }

            if (commonData.cases.Select(_ => _.Value).Any(_ => _.tumor_rna_allcount_filename == ""))
            {
                Console.WriteLine("Not all cases have tumor RNA allcount files.");
                return;
            }

            isoformMap = new ASETools.IsoformMap(commonData.geneLocationInformation.genesByName);

            // initialize the data structure
            foreach (var geneInfo in commonData.geneLocationInformation.genesByName.Select(x => x.Value).ToList())
            {
                isoformBalanceBySpliceosome.Add(geneInfo.hugoSymbol, new Dictionary<string, Dictionary<bool, List<double>>>());
                foreach (var isoformName in geneInfo.isoforms.Select(x => x.ucscId))
                {
                    isoformBalanceBySpliceosome[geneInfo.hugoSymbol].Add(isoformName, new Dictionary<bool, List<double>>());
                    foreach (var mutantSpliceosome in ASETools.BothBools)
                    {
                        isoformBalanceBySpliceosome[geneInfo.hugoSymbol][isoformName].Add(mutantSpliceosome, new List<double>());
                    } // mutant
                } // isoforms
            } // genes

            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, int>(commonData.cases.Select(_ => _.Value).ToList(), HandleOneCase, null, null, 100);
            Console.WriteLine("Processing allcount files, 1 dot/100 cases:");
            ASETools.PrintNumberBar(commonData.cases.Count() / 100);
            threading.run();

            int minSamplesPerClass = 10;

            int nTooFewSamples = 0;
            int nConsidered = 0;
            var results = new List<NameAndPValue>();

            foreach (var geneInfo in commonData.geneLocationInformation.genesByName.Select(x => x.Value).ToList())
            {
                foreach (var isoformName in geneInfo.isoforms.Select(x => x.ucscId))
                {
                    if (isoformBalanceBySpliceosome[geneInfo.hugoSymbol][isoformName].Select(_ => _.Value.Count()).Min() < minSamplesPerClass)
                    {
                        nTooFewSamples++;
                        continue;
                    }

                    var items = new List<DoubleAndBool>();
                    foreach (var mutant in ASETools.BothBools)
                    {
                        isoformBalanceBySpliceosome[geneInfo.hugoSymbol][isoformName][mutant].ForEach(_ => items.Add(new DoubleAndBool(mutant, _)));
                    }

                    nConsidered++;

                    bool enoughData, reversed;
                    double nFirstGroup, nSecondGroup, U, z;

                    ASETools.MannWhitney<DoubleAndBool>.ComputeMannWhitney(items, items[0], _ => _.whichClass, _ => _.value, out enoughData, out reversed, out nFirstGroup, out nSecondGroup, out U, out z);
                    


                }
            }
        } // Main

        static void HandleOneCase(ASETools.Case case_, int state)
        {
            var allcountReader = new ASETools.AllcountReader(case_.tumor_rna_allcount_filename);

            if (!allcountReader.openFile())
            {
                Console.WriteLine("Unable to open allcount file " + case_.tumor_rna_allcount_filename);
                return;
            }

            var valuesForThisCase = new Dictionary<string, Dictionary<string, int>>();  // gene->isoform->totalReadCount
            foreach (var geneInfo in commonData.geneLocationInformation.genesByName.Select(x => x.Value).ToList())
            {
                valuesForThisCase.Add(geneInfo.hugoSymbol, new Dictionary<string, int>());
                foreach (var isoformName in geneInfo.isoforms.Select(x => x.ucscId))
                {
                    valuesForThisCase[geneInfo.hugoSymbol].Add(isoformName, 0);
                } // isoform
            } // gene

            allcountReader.ReadAllcountFile((contigName, location, currentMappedReadCount) => processBase(valuesForThisCase, contigName, location, currentMappedReadCount));

            var mafLines = ASETools.MAFLine.ReadFile(case_.all_maf_lines_filename, "", false);
            bool anySpliceosomeMutations = mafLines.Any(_ => !_.IsSilent() && ASETools.spliceosome_genes.Contains(_.Hugo_Symbol));

            lock (isoformBalanceBySpliceosome)
            {
                foreach (var hugo_symbol in valuesForThisCase.Select(_ => _.Key).ToList())
                {
                    int totalIsoformReads = valuesForThisCase[hugo_symbol].Select(_ => _.Value).Sum();  // This is more than the actual mapped reads, since some exons are in more than one isoform.  Hence "isoformReads"
                    foreach (var isoformName in valuesForThisCase[hugo_symbol].Select(_ => _.Key))
                    {
                        isoformBalanceBySpliceosome[hugo_symbol][isoformName][anySpliceosomeMutations].Add((double)valuesForThisCase[hugo_symbol][isoformName] / totalIsoformReads);
                    }
                }
            }

        }   // HandleOneCase

        static void processBase(Dictionary<string, Dictionary<string, int>> valuesForThisCase, string contigName, int location, int currentMappedReadCount)
        {
            foreach (var isoform in isoformMap.getIsoformsMappedTo(contigName, location))
            {
                valuesForThisCase[isoform.hugo_symbol][isoform.ucscId] += currentMappedReadCount;
            } // foreach isoform
        } // processBase
    }
}
