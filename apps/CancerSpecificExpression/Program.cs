using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace CancerSpecificExpression
{
    class Program
    {
        class Region
        {
            public double meanDividedByStdDev = 0;
            public double meanMinusStdDev = 0;
            public double meanMean = 0;
            public double meanStdDev = 0;
            public string chromosome = "";
            public int locusGroup = 0;
        }


        static void Main(string[] args)
        {
            const int locusGroupSize = 10000;    // # of bases to group together

            // diseaseName -> (chromosome -> (locus rounded down -> Region))
            var allDiseases = new Dictionary<string, Dictionary<string, Dictionary<int, Region>>>();

            foreach (var disease in args)
            {
                var reader = new StreamReader(@"\sequence\reads\tcga\expression\expression_" + disease);

                allDiseases.Add(disease, new Dictionary<string, Dictionary<int, Region>>());

                string line;
                int currentLocusGroup = -1;
                Dictionary<string, Dictionary<int, Region>> diseaseRecord = allDiseases[disease];
                Dictionary<int, Region> chromosomeRecord = null;
                string chromosomeName = "";
                Region region = null;


                while (null != (line = reader.ReadLine()))
                {
                    var fields = line.Split('\t');
                    if (fields.Count() == 1)
                    {
                        //
                        // Changing chromosome
                        //
                        chromosomeName = fields[0].ToLower();
                        diseaseRecord.Add(chromosomeName, new Dictionary<int,Region>());
                        chromosomeRecord = diseaseRecord[chromosomeName];
                        currentLocusGroup = -1;

                        continue;
                    }

                    int locus = Convert.ToInt32(fields[0]);
                    int n = Convert.ToInt32(fields[1]);
                    double mean = Convert.ToDouble(fields[2]);
                    double stdDev = Convert.ToDouble(fields[3]);

                    int locusGroup = (locus / locusGroupSize) * locusGroupSize;

                    if (locusGroup != currentLocusGroup)
                    {
                        region = new Region();
                        chromosomeRecord.Add(locusGroup, region);
                        region.chromosome = chromosomeName;
                        region.locusGroup = locusGroup;
                        currentLocusGroup = locusGroup;
                    }

                    region.meanMinusStdDev += mean - stdDev;
                    region.meanMean += mean;
                    region.meanStdDev += stdDev;
                    if (stdDev <= 0)
                    {
                        Console.WriteLine("Zero (or negative!) std dev " + disease + " " + chromosomeName + " " + locus);
                    }
                    else
                    {
                        region.meanDividedByStdDev += mean / stdDev;
                    }
                }

            } // for each disease listed on the command line

            //
            // Build the aggregate
            //
            var aggregate =  new Dictionary<string, Dictionary<int, Region>>();

            int nNonMainDiseases = args.Count() - 1;
            var mainDisease = allDiseases[args[0]];

            foreach (var diseaseEntry in allDiseases)
            {
                if (diseaseEntry.Key == args[0])
                {
                    //
                    // Don't add in the main disease
                    //
                    continue;
                }
                foreach (var chromosomeEntry in diseaseEntry.Value)
                {
                    string chromosomeName = chromosomeEntry.Key;
                    if (!aggregate.ContainsKey(chromosomeName))
                    {
                        aggregate.Add(chromosomeName, new Dictionary<int, Region>());
                    }

                    foreach (var regionEntry in chromosomeEntry.Value) {
                        Region region;
                        if (!aggregate[chromosomeName].ContainsKey(regionEntry.Key))
                        {
                            region = new Region();
                            aggregate[chromosomeName].Add(regionEntry.Key, region);
                            region.locusGroup = regionEntry.Value.locusGroup;
                            region.chromosome = regionEntry.Value.chromosome;

                        }
                        else
                        {
                            region = aggregate[chromosomeName][regionEntry.Key];
                        }

                        region.meanDividedByStdDev += regionEntry.Value.meanDividedByStdDev / nNonMainDiseases;
                        region.meanMean += regionEntry.Value.meanMean / nNonMainDiseases;
                        region.meanStdDev += regionEntry.Value.meanStdDev / nNonMainDiseases;
                        region.meanMinusStdDev += regionEntry.Value.meanMinusStdDev / nNonMainDiseases;
                    }
                }
            }

            var excessRegions = new List<Region>();
            foreach (var chromosomeEntry in mainDisease)
            {
                if (aggregate.ContainsKey(chromosomeEntry.Key))
                {
                    foreach (var regionEntry in chromosomeEntry.Value)
                    {
                        if (aggregate[chromosomeEntry.Key].ContainsKey(regionEntry.Key))
                        {
                            var excessRegion = new Region();
                            excessRegion.chromosome = chromosomeEntry.Key;
                            excessRegion.locusGroup = regionEntry.Key;
                            excessRegion.meanDividedByStdDev = regionEntry.Value.meanDividedByStdDev / aggregate[chromosomeEntry.Key][regionEntry.Key].meanDividedByStdDev;
                            excessRegion.meanMinusStdDev = regionEntry.Value.meanMinusStdDev / aggregate[chromosomeEntry.Key][regionEntry.Key].meanMinusStdDev;
                            excessRegions.Add(excessRegion);
                        }
                    }
                }
            }


            excessRegions.Sort(delegate(Region x, Region y)
            {
                if (x.meanDividedByStdDev > y.meanDividedByStdDev) return -1;
                if (x.meanDividedByStdDev < y.meanDividedByStdDev) return 1;
                return 0;
            });

            Console.WriteLine("Results for mean / std dev");
            for (int i = 0; i < 200; i++)
            {
                Console.WriteLine(excessRegions[i].chromosome + ":" + excessRegions[i].locusGroup);
            }


            excessRegions.Sort(delegate(Region x, Region y)
            {
                if (x.meanMinusStdDev > y.meanMinusStdDev) return -1;
                if (x.meanMinusStdDev < y.meanMinusStdDev) return 1;
                return 0;
            });

            Console.WriteLine("Results for mean - std dev");
            for (int i = 0; i < 200; i++)
            {
                Console.WriteLine(excessRegions[i].chromosome + ":" + excessRegions[i].locusGroup);
            }


        }
    }
}
