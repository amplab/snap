using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;

namespace PhylogenyFromASE
{
    class Program
    {

        enum ExpressionState { Unknown, Reference, Alt, Both };

        static Dictionary<string, Dictionary<string, ExpressionState>> expressionStateByCell = new Dictionary<string, Dictionary<string, ExpressionState>>();   // Barcode->locus->ExpressionState
        static Dictionary<string, Dictionary<string, ExpressionState>> expressionStateByLocus = new Dictionary<string, Dictionary<string, ExpressionState>>();  // Locus->Barcode->ExpressionState

        static int ChromosomeFromLocus(string locus)
        {
            return Convert.ToInt32(locus.Split(':')[0]);
        }

        struct BarcodeAndChromosome
        {
            public BarcodeAndChromosome(string barcode_, string chromosome_)
            {
                barcode = barcode_;
                chromosome = chromosome_;
            }

            public readonly string barcode;
            public readonly string chromosome;
        }

        static void Main(string[] args)
        {
            var configuration = ASETools.Configuration.loadFromFile(args);

            var inputFilename = configuration.baseDirectory + "DJA20_WESConfSNV_PassHardFilt_Tile.txt";
            var inputFile = ASETools.CreateStreamReaderWithRetry(inputFilename);

            if (null == inputFile)
            {
                Console.WriteLine("Unable to open input file: " + inputFilename);
                return;
            }

            //
            // Format is, with all columns tab separated:
            // Header line: blank column, list of cell barcodes
            // Other lines: genomic locus in the form n:l (where n is the chromosome number and l is the locus within the chromosome in hg19), then one column per cell with 1 for no data, 3 ref only, 4 both, 5 alt only
            //

            var header = inputFile.ReadLine();
            var headerColumns = header.Split('\t');

            if (headerColumns.Count() < 2 || headerColumns[0] != "")
            {
                Console.WriteLine("Header didn't look right: " + header);
                return;
            }

            var barcodes = headerColumns.ToList();
            barcodes.RemoveAt(0);

            if (barcodes.Any(_ => _ == ""))
            {
                Console.WriteLine("Found null barcode in header: " + header);
                return;
            }

            foreach (var barcode in barcodes)
            {
                if (expressionStateByCell.ContainsKey(barcode))
                {
                    Console.WriteLine("Duplicate barcode: " + barcode);
                    return;
                }

                expressionStateByCell.Add(barcode, new Dictionary<string, ExpressionState>());
            }

            //
            // Now read the rest of the lines.
            //
            string line;
            var loci = new List<string>();

            while (null != (line = inputFile.ReadLine()))
            {
                if (line == "")
                {
                    continue;   // There's a blank line at the end, skip it.
                }

                var fields = line.Split('\t');
                if (fields.Count() != headerColumns.Count())
                {
                    Console.WriteLine("Found line with wrong number of columns: " + line);
                    return;
                }

                var locus = fields[0];
                if (expressionStateByLocus.ContainsKey(locus))
                {
                    Console.WriteLine("Duplicate locus " + locus);
                    return;
                }

                loci.Add(locus);

                expressionStateByLocus.Add(locus, new Dictionary<string, ExpressionState>());

                for (int i = 0; i < barcodes.Count(); i++)
                {
                    ExpressionState expressionState;
                    var stateInText = fields[i + 1];
                    var barcode = barcodes[i];
                    if (stateInText == "1")
                    {
                        expressionState = ExpressionState.Unknown;
                    } else if (stateInText == "3")
                    {
                        expressionState = ExpressionState.Reference;
                    } else if (stateInText == "4")
                    {
                        expressionState = ExpressionState.Both;
                    } else if (stateInText == "5")
                    {
                        expressionState = ExpressionState.Alt;
                    } else
                    {
                        Console.WriteLine("Invalid expression state cell value for locus " + locus + ", barcode " + barcode + ": " + stateInText);
                        return;
                    }

                    expressionStateByCell[barcode].Add(locus, expressionState);
                    expressionStateByLocus[locus].Add(barcode, expressionState);
                }
            } // For every non-header line 

            int nLoci = loci.Count();

            inputFile.Close();

            int nDifferentiating = 0;
            int nNoData = 0;
            int totalLociBoth = 0;
            int totalLociSingle = 0;

            int totalMeasurementsRef = 0;
            int totalMeasurementsAlt = 0;
            int totalMeasurementsBoth = 0;
            int totalMeasurementsUnknown = 0;
 

            var cellCountHistogramForDifferentiatingLoci = new ASETools.PreBucketedHistogram(0, 100, 1);
            var differentiatingPowerByLocusHistogram = new ASETools.PreBucketedHistogram(0, 100, 1);

            var outputFile = ASETools.CreateStreamWriterWithRetry(@"f:\temp\phylogeny_from_ase_output.txt");

            foreach (var locus in loci)
            {
                int nRef = expressionStateByLocus[locus].Where(_ => _.Value == ExpressionState.Reference).Count();
                int nAlt = expressionStateByLocus[locus].Where(_ => _.Value == ExpressionState.Alt).Count();
                int nBoth = expressionStateByLocus[locus].Where(_ => _.Value == ExpressionState.Both).Count();

                totalMeasurementsAlt += nAlt;
                totalMeasurementsRef += nRef;
                totalMeasurementsBoth += nBoth;
                totalMeasurementsUnknown += expressionStateByLocus[locus].Where(_ => _.Value == ExpressionState.Unknown).Count();

                int nStates = 0;
                if (nRef > 0) nStates++;
                if (nAlt > 0) nStates++;
                if (nBoth > 0) nStates++;


                if (nStates > 1)
                {
                    cellCountHistogramForDifferentiatingLoci.addValue(nRef + nAlt + nBoth);

                    int[] values = { nRef, nAlt, nBoth };
                    Array.Sort(values);

                    var differentiatingPower = Math.Min(values[1] + values[0], values[2]); // The smaller of the most common one and not the most common one
                    differentiatingPowerByLocusHistogram.addValue(differentiatingPower);

                    if (differentiatingPower >= 50)
                    {
                        outputFile.WriteLine(locus + " has differentiating power " + differentiatingPower);
                    }

                    nDifferentiating++;
                }
                else if (nStates == 0)
                {
                    nNoData++;
                }
                else if (nBoth > 0)
                {
                    totalLociBoth++;
                }
                else
                {
                    totalLociSingle++;
                }
            } // locus

            outputFile.WriteLine(totalMeasurementsAlt + " alt measurements, " + totalMeasurementsRef + " ref measurements, " + totalMeasurementsBoth + " both measurements, and " + totalMeasurementsUnknown + " with no data");

            outputFile.WriteLine("Distribution of differentiating loci by cells with values for that locus");
            cellCountHistogramForDifferentiatingLoci.WriteHistogram(outputFile);

            outputFile.WriteLine();
            outputFile.WriteLine("Distribution of differentiating loci by cells excluded by it");
            differentiatingPowerByLocusHistogram.WriteHistogram(outputFile);

            //
            // Look for cnLOH-type events, which here we characterize by a run of all single allele from the beginning or end of the chromosome of at least 5 loci.
            //

            outputFile.WriteLine();
            outputFile.WriteLine("Barcode\tchromosome\tBeginning/End/Whole\tCount");

            int nLociToDetermineCNLOH = 5;

            var wholeChromsomeLOH = new List<BarcodeAndChromosome>();

            int nLOH = 0;
            foreach (var barcode in barcodes)
            {
                var currentChromosome = "";
                int nFromBeginning = 0;
                int nToEnd = 0;
                bool seenAnyNonLOH = true;
                for (int i = 0; i < nLoci; i++)
                {
                    var locus = loci[i];
                    var chromosome = locus.Split(':')[0];

                    if (currentChromosome != chromosome)
                    {
                        if (!seenAnyNonLOH)
                        {
                            if (nFromBeginning >= nLociToDetermineCNLOH)
                            {
                                //
                                // Whole chromosome LOH
                                //
                                nLOH++;

                                outputFile.WriteLine(barcode + "\t" + currentChromosome + "\tWhole\t" + nFromBeginning);

                                wholeChromsomeLOH.Add(new BarcodeAndChromosome(barcode, currentChromosome));
                            }
                        } else
                        {
                            if (nFromBeginning >= nLociToDetermineCNLOH)
                            {
                                nLOH++;
                                outputFile.WriteLine(barcode + "\t" + currentChromosome + "\tBeginning\t" + nFromBeginning);
                            }
                            if (nToEnd >= nLociToDetermineCNLOH)
                            {
                                nLOH++;
                                outputFile.WriteLine(barcode + "\t" + currentChromosome + "\tEnd\t" + nToEnd);
                            }
                        }

                        nFromBeginning = 0;
                        seenAnyNonLOH = false;
                        nToEnd = 0;

                        currentChromosome = chromosome;
                    }

                    var expressionState = expressionStateByCell[barcode][locus];
                    if (expressionState == ExpressionState.Unknown)
                    {
                        continue;
                    }

                    if (expressionState == ExpressionState.Both)
                    {
                        seenAnyNonLOH = true;
                        nToEnd = 0;
                    } else
                    {
                        if (!seenAnyNonLOH)
                        {
                            nFromBeginning++;
                        }
                        nToEnd++;
                    }
                } // loci
            } // barcode

            for (int chromosomeNumber = 1; chromosomeNumber <= ASETools.nHumanAutosomes; chromosomeNumber++)
            {
                var candidateBarcodes = wholeChromsomeLOH.Where(_ => _.chromosome == chromosomeNumber.ToString()).ToList();
                if (candidateBarcodes.Count() > 1)
                {
                    var candidateLoci = loci.Where(locus => ChromosomeFromLocus(locus) == chromosomeNumber && candidateBarcodes.Select(x => x.barcode).Where(barcode => expressionStateByCell[barcode][locus] != ExpressionState.Unknown).Count() > 1).ToList();

                    if (candidateLoci.Count() > 0)
                    {
                        outputFile.WriteLine();
                        outputFile.WriteLine("Cells with whole chromosome LOH for chromosome " + chromosomeNumber + " and loci in that chromosome with at least two of them being determined");
                        foreach (var barcode in candidateBarcodes.Select(_ => _.barcode))
                        {
                            outputFile.Write("\t" + barcode);
                        }
                        outputFile.WriteLine();
                        foreach (var locus in candidateLoci)
                        {
                            outputFile.Write(locus);
                            foreach (var barcode in candidateBarcodes.Select(_ => _.barcode))
                            {
                                if (expressionStateByCell[barcode][locus] == ExpressionState.Unknown)
                                {
                                    outputFile.Write("\t");
                                } else if (expressionStateByCell[barcode][locus] == ExpressionState.Alt)
                                {
                                    outputFile.Write("\tA");
                                } else if (expressionStateByCell[barcode][locus] == ExpressionState.Reference)
                                {
                                    outputFile.Write("\tR");
                                } else
                                {
                                    outputFile.Write("\tB");    // This shouldn't happen, since these are all whole chromosome LOH cases.
                                }
                            } // barcode
                            outputFile.WriteLine();
                        } // locus
                    } // if we have enough loci


                } // if we have enough cells
            } // chromosome

            var perChromosomeMesauredSitesHistogram = new ASETools.PreBucketedHistogram(0, 100, 1);
            foreach (var barcode in barcodes)
            {
                for (int chromosomeNumber = 1; chromosomeNumber <= ASETools.nHumanAutosomes; chromosomeNumber++)
                {
                    perChromosomeMesauredSitesHistogram.addValue(expressionStateByCell[barcode].Where(_ => chromosomeNumber == ChromosomeFromLocus(_.Key) && _.Value != ExpressionState.Unknown).Count());
                } // chromosome
            } // cell

            outputFile.WriteLine("Distribution of chromosomes by count of measured sites with known values.");
            perChromosomeMesauredSitesHistogram.WriteHistogram(outputFile);


            outputFile.WriteLine("**done**");
            outputFile.Close();

        } // Main
    }
}
