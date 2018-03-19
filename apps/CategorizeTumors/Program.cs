using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Diagnostics;
using ASELib;

namespace CategorizeTumors
{
    class Program
    {

        class Result : IComparable<Result>
        {
            public Result(string hugo_symbol_)
            {
                hugo_symbol = hugo_symbol_;
            }
            public readonly string hugo_symbol;
            public int nZero = 0;
            public int nMultiple = 0;
            public int nLossOfHeterozygosity = 0;
            public int nOnePlusASE = 0;
            public int nOnePlusReverseASE = 0;
            public int nMinorSubclone = 0;
            public int nNonsenseMediatedDecay = 0;
            public int nSingle = 0;
            public int nTooFewReads = 0;

            public int totalEvaluated()
            {
                return nZero + nMultiple + nLossOfHeterozygosity + nOnePlusASE + nOnePlusReverseASE + nMinorSubclone + nNonsenseMediatedDecay + nSingle;
            }

            public int totalMutant()    // Excluding minor subclones
            {
                return nMultiple + nLossOfHeterozygosity + nOnePlusASE + nOnePlusReverseASE + nSingle + nNonsenseMediatedDecay;
            }

            string frac(int numerator)
            {
                if (totalEvaluated() == 0)
                {
                    return "*";
                }

                return Convert.ToString((double)numerator / totalEvaluated());
            }

            string fracOfAllMutant(int numerator)
            {
                if (totalMutant() == 0)
                {
                    return "*";
                }

                return Convert.ToString((double)numerator / totalMutant());
            }
            public string fracZero() { return frac(nZero); }
            public string fracMultiple() { return frac(nMultiple); }
            public string fracLossOfHeterozygosity() { return frac(nLossOfHeterozygosity); }
            public string fracOnePlusASE() { return frac(nOnePlusASE); }
            public string fracOnePlusReverseASE() { return frac(nOnePlusReverseASE); }
            public string fracMinorSubclone() { return frac(nMinorSubclone); }
            public string fracNonsenseMediatedDecay() { return frac(nNonsenseMediatedDecay); }
            public string fracSingle() { return frac(nSingle); }

            public string fracOfAllMutantMultiple() { return fracOfAllMutant(nMultiple); }
            public string fracOfAllMutantLossOfHeterozygosity() { return fracOfAllMutant(nLossOfHeterozygosity); }
            public string fracOfAllMutantOnePlusASE() { return fracOfAllMutant(nOnePlusASE); }
            public string fracOfAllMutantOnePlusReverseASE() { return fracOfAllMutant(nOnePlusReverseASE); }
            public string fracOfAllMutantSingle() { return fracOfAllMutant(nSingle); }
            public string fracOfAllMutantNonsenseMediatedDecay() { return fracOfAllMutant(nNonsenseMediatedDecay); }

            public int CompareTo(Result peer)
            {
                return hugo_symbol.CompareTo(peer.hugo_symbol);
            }
        }

        static void Main(string[] args)
        {
            var stopwatch = new Stopwatch();
            stopwatch.Start();

            var configuration = ASETools.Configuration.loadFromFile(args);
            if (null == configuration)
            {
                Console.WriteLine("Unable to load configuration.");
                return;
            }

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname).Select(x => x.Value).ToList();

            if (null == cases)
            {
                Console.WriteLine("Unable to load cases.  You must generate it before running this tool.");
                return;
            }

            var results = new List<Result>();

            var scatterGraphLinesByGene = ASETools.GeneScatterGraphLine.LoadAllGeneScatterGraphLines(configuration.geneScatterGraphsDirectory, true, "*").
                Where(x => x.Variant_Classification != "Silent" && x.Chromosome != "chrX" && x.Chromosome != "chrY" && x.Chromosome != "chrM" && x.Chromosome != "chrMT").
                GroupByToDict(x => x.Hugo_Symbol);    // true->from unfiltered

            Console.WriteLine("Loaded " + scatterGraphLinesByGene.Count() + " scatter graph genes with " + scatterGraphLinesByGene.Select(x => x.Value.Count()).Sum() + " lines in " + ASETools.ElapsedTimeInSeconds(stopwatch));


            Console.Write("Processing " + scatterGraphLinesByGene.Count() + " genes, 1 dot/100 genes: ");
            int nGenesProcessed = 0;

            foreach (var geneEntry in scatterGraphLinesByGene)
            {
                var hugo_symbol = geneEntry.Key;
                var scatterGraphLines = geneEntry.Value;

                var result = new Result(hugo_symbol);

                foreach (var case_ in cases)
                {
                    var linesForThisCase = scatterGraphLines.Where(x => x.case_id == case_.case_id).ToList();

                    if (linesForThisCase.Count() == 0)
                    {
                        result.nZero++;
                        continue;
                    }

                    if (linesForThisCase.Count() > 1)
                    {
                        result.nMultiple++;
                        continue;
                    }

                    var line = linesForThisCase[0];

                    if (line.tumorDNAReadCounts.usefulReads() < configuration.minDNAReadCoverage)
                    {
                        result.nTooFewReads++;
                        continue;
                    }

                    if (line.tumorDNAReadCounts.AltFraction() < 0.4)
                    {
                        result.nMinorSubclone++;
                        continue;
                    }

                    if (line.tumorDNAReadCounts.AltFraction() > 0.6)
                    {
                        result.nLossOfHeterozygosity++;
                        continue;
                    }

                    if (line.tumorRNAReadCounts.usefulReads() < configuration.minRNAReadCoverage)
                    {
                        result.nTooFewReads++;
                        continue;
                    }

                    if (ASETools.NonsenseMediatedDecayCausingVariantClassifications.Contains(line.Variant_Classification))
                    {
                        result.nNonsenseMediatedDecay++;
                        continue;
                    }

                    if (line.tumorRNAReadCounts.AltFraction() > 0.6)
                    {
                        result.nOnePlusASE++;
                        continue;
                    }

                    if (line.tumorRNAReadCounts.AltFraction() < 0.4)
                    {
                        result.nOnePlusReverseASE++;
                        continue;
                    }

                    result.nSingle++;                    
                } // foreach case

                results.Add(result);

                nGenesProcessed++;
                if (nGenesProcessed % 100 == 0)
                {
                    Console.Write(".");
                }
            } // foreach gene

            Console.WriteLine(ASETools.ElapsedTimeInSeconds(stopwatch));

            results.Sort();

            var outputFile = ASETools.CreateStreamWriterWithRetry(configuration.finalResultsDirectory + ASETools.Configuration.geneCategorizationFilename);

            outputFile.WriteLine("Hugo Symbol\tn Too Few Reads\tn Evaluated\tn Mutant (ignoring minor subclones)\tn No Mutations\tfrac No Mutations\tn Minor Subclone\tfrac Minor Subclone\tn Multiple Mutations\tfrac Multiple Mutations\tfrac of all mutant multiple mutations"
                + "\tn Loss of Heterozygosity\tfrac Loss of Heterozygosity\tfrac of all mutant loss of heterozygosity"
                + "\tn nonsense mediated decay\tfrac nonsense mediated decay\tfrac of all mutant nonsense medidated decay"
                + "\tn Reverse ASE\tfrac Reverse ASE\tfrac of all mutant reverse ASE\tn ASE\tfrac ASE\tfrac of all mutant ASE\tn Single\tfrac Single\tfrac of all mutant single\tgraph title\theading1\theading2\theading3\theading4\theading5\theading6\theading7\theading8");

            foreach (var result in results)
            {
                outputFile.WriteLine(ASETools.ConvertToExcelString(result.hugo_symbol) + "\t" + result.nTooFewReads + "\t" + result.totalEvaluated() + "\t" + result.totalMutant()
                    + "\t" + result.nZero + "\t" + result.fracZero()
                    + "\t" + result.nMinorSubclone + "\t" + result.fracMinorSubclone()
                    + "\t" + result.nMultiple + "\t" + result.fracMultiple() + "\t" + result.fracOfAllMutantMultiple()
                    + "\t" + result.nLossOfHeterozygosity + "\t" + result.fracLossOfHeterozygosity() +"\t" + result.fracOfAllMutantLossOfHeterozygosity()
                    + "\t" + result.nNonsenseMediatedDecay + "\t" + result.fracNonsenseMediatedDecay() + "\t" + result.fracOfAllMutantNonsenseMediatedDecay()
                    + "\t" + result.nOnePlusReverseASE + "\t" + result.fracOnePlusReverseASE() + "\t" + result.fracOfAllMutantOnePlusReverseASE()
                    + "\t" + result.nOnePlusASE + "\t" + result.fracOnePlusASE() + "\t" + result.fracOfAllMutantOnePlusASE()
                    + "\t" + result.nSingle + "\t" + result.fracSingle() + "\t" + result.fracOfAllMutantSingle()
                    + "\tBreakdown of " + result.hugo_symbol + " mutant tumors excluding minor subclones\tNo mutations\tMinor subclone\tMultiple mutations\tLoss of heterozygosity\tNonsense mediated decay\t> 60% Wild type\t> 60% Mutant\tEven expression"
                    );
            }

            outputFile.Close();


        } // main
    }
}
