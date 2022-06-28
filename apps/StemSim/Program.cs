using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Management.Instrumentation;
using System.Text;
using System.Threading.Tasks;
using ASELib;

namespace StemSim
{
    internal class Program
    {
        class Cell
        {
            public Cell(int cellNumberOfOrigin_, bool leukemic_)
            {
                cellNumberOfOrigin = cellNumberOfOrigin_;
                leukemic = leukemic_;
            } // ctor

            public Cell(Cell peer)
            {
                cellNumberOfOrigin = peer.cellNumberOfOrigin;
                leukemic = peer.leukemic;
            } // ctor

            public int cellNumberOfOrigin;
            public bool leukemic;
        } // Cell

        static int[] leukemicPercentagesToTrack = { 5, 10, 25, 50, 75, 100 };
        static int[] cellPercentilesToTrack = { 5, 10, 25, 50, 75, 90, 95 };

        class Status
        {
            public long totalCells = 0;
            public long totalSurvivingLineages = 0;
            public int nReports = 0;
            public long totalLeukemic = 0;  // toal leukemic stem cells
            public int nLeukemic = 0; // n runs with any leukemic stem cells
            public Dictionary<int, int> nPercentLeukemic = new Dictionary<int, int>();
            public List<double> allCellCounts = new List<double>(); // All of the data, used for computing percentiles

            public Status()
            {
                foreach (var percentLeukemic in leukemicPercentagesToTrack)
                {
                    nPercentLeukemic.Add(percentLeukemic, 0);
                }
            } // ctor

            public int minCells = 1000000;
            public int maxCells = 0;

            public double percentileCellCount(int percentile)
            {
                if (allCellCounts.Count() == 0)
                {
                    return 0;
                }

                allCellCounts.Sort(); // Doing this every time is redundant, but it's not close to performance critical.

                return allCellCounts[allCellCounts.Count() * percentile / 100];
            }
        } // Status


        static void Main(string[] args)
        {
            var random = new Random();

            int nInitialCells = 11000;
            int nInitialLeukemic = nInitialCells * 5 / 100;
            int cellDivideTimeInHours = 30 * 24;    // One division/cell/month
            double differentiationRate = 0.5;

            //
            // Bounds for applying regulation.
            //
            // We use a simple model: keep the cell count between min and max by
            // biasing the differentiation decision by linearly reducing the
            // chance of in/decereasing the count and applying that probability evenly
            // to the other two possible outcomes.  Between lower and upperRegulationStart
            // it's unbiased.
            //
            int minCellTarget = nInitialCells * 95 / 100;
            int maxCellTarget = nInitialCells * 105 / 100;
            int lowerRegulationStart = nInitialCells * 975 / 1000;  
            int upperRegulationStart = nInitialCells * 1025 / 1000;

            var timer = new Stopwatch();
            timer.Start();

            double finishTime = 80 * 365 * 24;  // 80 years in hours
            double statusPeriod = 30 * 24;

            var outputStreams = new List<StreamWriter>();
            outputStreams.Add(new StreamWriter(Console.OpenStandardOutput()));
            outputStreams.Add(ASETools.CreateStreamWriterWithRetry(@"c:\temp\StemSimOut3.txt"));


            var statuses = new List<Status>();
            for (int i = 0; i < finishTime / statusPeriod; i++)
            {
                statuses.Add(new Status());
            }

            int nIterations = 500;

            for (int iteration = 0; iteration < nIterations; iteration++)
            {
                double timeInHours = 0;
                int whichStatus = 0;
                int nextCellToDivide = 0;

                var cells = new List<Cell>();

                for (int i = 0; i < nInitialCells; i++)
                {
                    cells.Add(new Cell(i, i < nInitialLeukemic));
                }

                while (timeInHours < finishTime)
                {
                    var cellCount = cells.Count();

                    if (timeInHours >= whichStatus * statusPeriod)
                    {
                        var survivingCells = new HashSet<int>();
                        foreach (var cell in cells)
                        {
                            survivingCells.Add(cell.cellNumberOfOrigin);
                        }

                        int nSurvivingOriginalCells = survivingCells.Count();
                        int nLeukemicCells = cells.Where(_ => _.leukemic).Count();

                        statuses[whichStatus].nReports++;
                        statuses[whichStatus].totalCells += cellCount;
                        statuses[whichStatus].totalSurvivingLineages += nSurvivingOriginalCells;
                        statuses[whichStatus].totalLeukemic += nLeukemicCells;
                        if (nLeukemicCells != 0)
                        {
                            statuses[whichStatus].nLeukemic++;
                        }

                        foreach (var percentLeukemic in leukemicPercentagesToTrack)
                        {
                            if (nLeukemicCells >= (double)cellCount * percentLeukemic / 100.0)
                            {
                                statuses[whichStatus].nPercentLeukemic[percentLeukemic]++;
                            }
                        }

                        statuses[whichStatus].minCells = Math.Min(statuses[whichStatus].minCells, cellCount);
                        statuses[whichStatus].maxCells = Math.Max(statuses[whichStatus].maxCells, cellCount);
                        statuses[whichStatus].allCellCounts.Add(cellCount);

                        whichStatus++;
                    } // if we're reporting

                    timeInHours += (double)cellDivideTimeInHours / cellCount;
                    int cellToDivide = nextCellToDivide;// random.Next(cells.Count());

                    double differentiateOut = 0.25; // If the random number is less than this, then both daughters differentiate out
                    double bothStay = 0.75;         // If it's greater than this then neither daughter differentiates.

                    if (cellCount < lowerRegulationStart)
                    {
                        double newDifferentiateOut = ((double)cellCount - minCellTarget) / (lowerRegulationStart - minCellTarget) * differentiateOut;
                        bothStay -= (differentiateOut - newDifferentiateOut) / 3.0; // Allocate 1/3 of the probability given up by reduction to incerease, leaving 2/3 of it to staying the same (keeping the 2-1 ratio in the unregulated setting)
                        differentiateOut = newDifferentiateOut;
                    } else if (cellCount > upperRegulationStart)
                    {
                        double bothStayProbability = 1 - bothStay;  // Since both staying happens if the random number is >= bothStay
                        double newBothStayProbability = ((double)maxCellTarget - cellCount) / (maxCellTarget - upperRegulationStart) * bothStayProbability;

                        differentiateOut -= (bothStayProbability - newBothStayProbability) / 3.0; // See above for comment
                        bothStay = 1.0 - newBothStayProbability;
                    }


                    double randomDouble = random.NextDouble();
                    if (randomDouble < differentiateOut)
                    {
                        cells.RemoveRange(cellToDivide, 1);
                    } 
                    else if (randomDouble >= bothStay)
                    {
                        cells.Add(new Cell(cells[cellToDivide]));
                    }

                    if (cells.Count() == 0)
                    {
                        //
                        // No blood.  Dead.
                        //
                        break;
                    }

                    nextCellToDivide = (nextCellToDivide + 1) % cells.Count();

                } // while we have time left to simulate
            } // iteration

            Console.WriteLine("Run took " + ASETools.ElapsedTime(timer));

            ASETools.WriteToMultipleStreams(outputStreams, "Time (years)\tn stem cells\tmin %stem cells\tmax %stem cells\tn original cells with descendents\tfraction original with descendents\tfraction dead\tfraction leukemic cells\tfraction with any leukemia");
            foreach (var percentLeukemic in leukemicPercentagesToTrack)
            {
                ASETools.WriteToMultipleStreams(outputStreams, "\tfraction " + percentLeukemic + "% leukemic");
            }

            foreach (var percentile in cellPercentilesToTrack)
            {
                ASETools.WriteToMultipleStreams(outputStreams, "\t" + percentile + "th percentile overall cell count");
            }


            ASETools.WriteLineToMultipleStreams(outputStreams);

            for (int i = 0; i < statuses.Count(); i++)
            {
                ASETools.WriteToMultipleStreams(outputStreams, (double)i / 12 + "\t" + statuses[i].totalCells / statuses[i].nReports + "\t" + (double)statuses[i].minCells / nInitialCells + "\t" +
                                                    (double)statuses[i].maxCells / nInitialCells + "\t" + statuses[i].totalSurvivingLineages / statuses[i].nReports + "\t" +
                                                    (double)statuses[i].totalSurvivingLineages / statuses[i].nReports / nInitialCells + "\t" + ((double)nIterations - statuses[i].nReports) / nIterations + "\t" +
                                                    (double)statuses[i].totalLeukemic / statuses[i].totalCells + "\t" + (double)statuses[i].nLeukemic / statuses[i].nReports);

                foreach (var percentLeukemic in leukemicPercentagesToTrack)
                {
                    ASETools.WriteToMultipleStreams(outputStreams, "\t" + (double)statuses[i].nPercentLeukemic[percentLeukemic] / statuses[i].nReports);
                }

                foreach (var percentile in cellPercentilesToTrack)
                {
                    ASETools.WriteToMultipleStreams(outputStreams, "\t" + statuses[i].percentileCellCount(percentile) / nInitialCells);
                }

                ASETools.WriteLineToMultipleStreams(outputStreams);
            } // for each status

            outputStreams.ForEach(_ => _.Close());

        } // Main
    } // Program
} // namespace StemSim
