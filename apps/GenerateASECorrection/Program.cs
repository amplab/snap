using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;
using System.Diagnostics;
using MathNet.Numerics;

namespace GenerateASECorrection
{
    class Program
    {
        static double generateASE(Random random, double[] cdf)
        {
            var randomValue = random.NextDouble();

            for (int i = 0; i < cdf.Count(); i++)
            {
                if (cdf[i] > randomValue)
                {
                    return (double)i / (cdf.Count() - 1);
                }
            }
            return 1;
        } // generateASE

        static double generateMeasuredASE(Random random, double trueASE, int readCount)
        {
            int a = 0, b = 0;

            for (int i = 0; i < readCount; i++)
            {
                if (random.NextDouble() <= trueASE / 2 + .5)
                {
                    a++;
                }
                else
                {
                    b++;
                }
            }

            return (double)Math.Abs(a - b) / (a + b);
        }

        static int generateReadDepth(Random random, int minReadDepth, double[] cdf)
        {
            int result;

            do
            {
                var value = random.NextDouble();

                for (result = 0; result < cdf.Count() - 1; result++)
                {
                    if (cdf[result] > value)
                    {
                        break;
                    }
                }
            } while (result < minReadDepth);

            return result;
        }

        //
        // Bucketizes a ASE p into one of g buckets.
        //
        static int Dindex(double p, int g)
        {
            return (int)Math.Round(p * (g - 1));
        }

        static double D(double p, int g)
        {
            return (double)Dindex(p, g) / (g - 1);
        }


        static double[,] CorrectionArray; // Indexed by read depth then ASE index
        static int nDiscreteASEValues;
        static int maxReadDepth;
        static int nIterations = 1000000;

        static void ComputeOneRow(int readDepth, int unused)
        {
            var random = new Random();
            var values = new ASETools.RunningMeanAndStdDev[nDiscreteASEValues];
            for (int i = 0; i < nDiscreteASEValues; i++)
            {
                values[i] = new ASETools.RunningMeanAndStdDev();
            }

            for (int iteration = 0; iteration < nIterations; iteration++)
            {
                var realASE = generateASE(random, assumedRealASEDistribution);

                var measuredASE = generateMeasuredASE(random, realASE, readDepth);

                values[Dindex(measuredASE, nDiscreteASEValues)].addValue(realASE);
            }

            lock (CorrectionArray)
            {
                for (var measuredASEIndex = 0; measuredASEIndex < nDiscreteASEValues; measuredASEIndex++)
                {
                    if (values[measuredASEIndex].getCount() == 0)
                    {
                        CorrectionArray[readDepth, measuredASEIndex] = -1;
                    }
                    else
                    {
                        CorrectionArray[readDepth, measuredASEIndex] = (double)measuredASEIndex / (nDiscreteASEValues - 1) - values[measuredASEIndex].getMeanAndStdDev().mean;
                    }
                }
            }
        }
        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            var configuration = ASETools.Configuration.loadFromFile(args);
            if (null == configuration)
            {
                Console.WriteLine("Unable to load configuration");
                return;
            }

            assumedRealASEDistribution = measuredASEDistribution = ASETools.Histogram.ReadFileToCDFValues(configuration.finalResultsDirectory + ASETools.TumorGermlineASEDistributionFilename);
            empericalReadDepthPDF = ASETools.Histogram.ReadFileToPDFValues(configuration.finalResultsDirectory + ASETools.TumorRNAReadDepthDistributionFilename);


            var pdfMatrix = ASETools.MeasuredASEMatrix.readFromFile(configuration.finalResultsDirectory + ASETools.UncorrectedOverallASEFilename, true, false); // Gets tumor germline RNA read depths

            if (null == pdfMatrix || null == assumedRealASEDistribution || null == empericalReadDepthPDF)
            {
                Console.WriteLine("Exiting because failed to read input data.");
                return;
            }

            maxReadDepth = empericalReadDepthPDF.Count() - 1;
            nDiscreteASEValues = assumedRealASEDistribution.Count();    // This is probably 101 (0->1 by 0.01)

            //
            // The read depth PDF includes weights for depths with too few reads.  Set them to zero and normalize the rest of the PDF.
            //
            double totalWeightForTooFewReads = 0;
            for (int i = 0; i < configuration.getMinReadCoverage(false); i++)
            {
                totalWeightForTooFewReads += empericalReadDepthPDF[i];
            }

            for (int i = 0; i < configuration.getMinReadCoverage(false); i++)
            {
                empericalReadDepthPDF[i] = 0;
            }

            for (int i = configuration.getMinReadCoverage(false); i < maxReadDepth;  i++)
            {
                empericalReadDepthPDF[i] /= (1 - totalWeightForTooFewReads);
            }


            for (int i = 0; i < configuration.commandLineArgs.Count(); i++)
            {
                if (configuration.commandLineArgs[i] == "-maxReadDepth" && i < configuration.commandLineArgs.Count() - 1) {
                    maxReadDepth = Convert.ToInt32(configuration.commandLineArgs[i + 1]);
                    i++;
                }
                else
                {
                    Console.WriteLine("Usage: GenerateASECorrection {-iterationsPerReadDepth iterationsPerReadDepth} {-totalSteps totalSteps} {-maxReadDepth maxReadDepth}");
                    return;
                }
            } // for each arg

            Console.WriteLine("" + nDiscreteASEValues + " discrete ASE values, " + maxReadDepth + " max read depth, " + empericalReadDepthPDF.Sum() + " total of the emperical read depth pdf.");

            for (int iteration = 0; iteration < /*10*/1; iteration++)
            {

                CorrectionArray = new double[maxReadDepth + 1, nDiscreteASEValues];

                var readDepths = new List<int>();

                for (int readDepth = configuration.getMinReadCoverage(false); readDepth <= maxReadDepth; readDepth++)
                {
                    readDepths.Add(readDepth);
                }

                var threading = new ASETools.WorkerThreadHelper<int, int>(readDepths, ComputeOneRow, null, null, 0);
                threading.run();



                if (false)
                {
                    //
                    // Now apply the correction to the assumed real ASE distribution to get the new assumed real distribution.
                    //
                    var newAssumedRealASEPDF = new double[nDiscreteASEValues];

                    var totalARASE = assumedRealASEDistribution.Sum();
                    var totalReadPDF = empericalReadDepthPDF.Sum();

                    for (int row = 0; row < nDiscreteASEValues; row++)
                    {
                        for (int readDepth = configuration.getMinReadCoverage(false); readDepth <= maxReadDepth; readDepth++)
                        {
                            var measuredASE = (double)row / (nDiscreteASEValues - 1);
                            var realASE = measuredASE - CorrectionArray[readDepth, row];
                            var realASEIndex = Dindex(realASE, nDiscreteASEValues);
                            var probabilityDensity = (row == 0) ? assumedRealASEDistribution[row] : (assumedRealASEDistribution[row] - assumedRealASEDistribution[row - 1]);    // AssumedRealASEDistribution is a cdf, but we want the pdf, so this converts.

                            newAssumedRealASEPDF[realASEIndex] += probabilityDensity * empericalReadDepthPDF[readDepth];
                        }
                    }
                    var totalNew = newAssumedRealASEPDF.Sum();

                    double difference = 0;
                    double total = 0;
                    for (int i = 0; i < nDiscreteASEValues; i++)
                    {
                        total += newAssumedRealASEPDF[i];
                        difference += Math.Abs(total - assumedRealASEDistribution[i]);
                        assumedRealASEDistribution[i] = total;
                    }

                    Console.WriteLine("" + ASETools.ElapsedTimeInSeconds(timer) + ": iteration " + iteration + " had difference of " + difference + " and new total of " + total);
                    assumedRealASEDistribution = newAssumedRealASEPDF;
                }



            } // iterations

            var outputFile = ASETools.CreateStreamWriterWithRetry(configuration.finalResultsDirectory + ASETools.ASECorrectionFilename);

            for (int readDepth = configuration.getMinReadCoverage(false); readDepth <= maxReadDepth; readDepth++)
            {
                outputFile.Write(readDepth);

                for (int aseIndex = 0; aseIndex < nDiscreteASEValues; aseIndex++)
                {
                    if (CorrectionArray[readDepth, aseIndex] == -1)
                    {
                        outputFile.Write("\t*");
                    }
                    else
                    {
                        outputFile.Write("\t" + (-1 * CorrectionArray[readDepth, aseIndex]));   // Reverse the sign because that's what the other code expects.
                    }
                }

                outputFile.WriteLine();
            }

            outputFile.Close();

            Console.WriteLine("Took " + ASETools.ElapsedTimeInSeconds(timer));

        } // Main

        static double[] measuredASEDistribution;

        static double[] assumedRealASEDistribution;

        static double[] empericalReadDepthPDF;

    } // Program
} // GenerateASECorrection
