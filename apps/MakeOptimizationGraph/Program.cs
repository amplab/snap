using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Reflection.Emit;
using System.Text;
using System.Threading.Tasks;
using ASELib;

namespace MakeOptimizationGraph
{
    internal class Program
    {
        class Optimization
        {
            public Optimization(string shortName_, string longName_, string runtimeFlag_)
            {
                shortName = shortName_;
                longName = longName_;
                runtimeFlag = runtimeFlag_;
            }

            public readonly string shortName;
            public readonly string longName;
            public readonly string runtimeFlag;
        } // Optimization

        class OptimizationSet
        {
            public OptimizationSet(List<Optimization> optimizations_, List<bool> optimizationsDisabled_)
            {
                optimizations = optimizations_;
                optimizationsDisabled = optimizationsDisabled_;
                hasG = false;
            } // ctor


            public OptimizationSet(List<Optimization> optimizations_) // ctor for -G-
            {
                optimizations = optimizations_;
                optimizationsDisabled = new List<bool>();
                foreach (var optimization in optimizations)
                {
                    optimizationsDisabled.Add(false);
                }

                hasG = true;
            }

            public string getShortName()
            {
                if (hasG)
                {
                    return "NoAffine";
                }

                if (optimizationsDisabled.All(_ => !_))
                {
                    return "None";
                }

                var retVal = "";
                for (int i = 0; i < optimizationsDisabled.Count; i++)
                {
                    if (optimizationsDisabled[i])
                    {
                        retVal += optimizations[i].shortName;
                    }
                }

                return retVal;
            } // getShortName

            public List<string> getFlags()
            {
                var retVal = new List<string>();

                if (hasG)
                {
                    retVal.Add("-G-");
                    return retVal;
                }

                for (int i = 0; i < optimizations.Count; i++)
                {
                    if (optimizationsDisabled[i])
                    {
                        retVal.Add(optimizations[i].runtimeFlag);
                    }
                }

                return retVal;
            } // getFlags

            public string getFlagsAsString()    // This is a printable sting suitable for graphs
            {
                if (hasG)
                {
                    return "NoAffine";
                }

                if (!optimizationsDisabled.Any(_ => _))
                {
                    return "None";
                }

                var flags = getFlags();
                string retVal = "";
                for (int i = 0; i < flags.Count(); i++)
                {
                    retVal += flags[i].Replace("-n","").ToUpper();
                }

                return retVal;

            } // getFlagsAsString

            public long getSortKey()
            {
                if (hasG)
                {
                    return -1;
                }

                long key = 0;

                for (int i = 0; i < optimizations.Count; i++)
                {
                    if (optimizationsDisabled[i]) key |= ((long)1 << i);
                }

                return key;
            } // getSortKey

            public readonly bool hasG;  // -G-, in which case optimizationsDisabled must be empty
            public readonly List<Optimization> optimizations;
            public readonly List<bool> optimizationsDisabled;
        } // OptimizationSet

        class Sample
        {
            public Sample(string name_, string shortName_, List<Optimization> optimizations_, List<OptimizationSet> allOptimizationSets)
            {
                name = name_;
                shortName = shortName_;
                optimizations = optimizations_;
                foreach (var optimizationsDisabled in allOptimizationSets.Select(_ => _.getSortKey()))
                {
                    runs.Add(optimizationsDisabled, new List<int>());
                }
            } // ctor

            public bool AddRun(string commandLine, string stats)
            {
                if (!commandLine.Contains(name))
                {
                    return false;
                }

                char[] justASpace = { ' ' };
                var fields = stats.Split(justASpace, StringSplitOptions.RemoveEmptyEntries);

                if (fields.Count() < 12)
                {
                    throw new FormatException("Sample.AddRun: not enough fields in stats line " + stats);
                }

                var runTime = Convert.ToInt32(fields[11].Replace(",", "")); // The replace() is because SNAP prints its times with commas, i.e., 11,000

                long optimizationsDisabled = 0;
                bool containsG = commandLine.Contains(" -G-");
                for (int i = 0; i < optimizations.Count(); i++)
                {
                    if (commandLine.Contains(optimizations[i].runtimeFlag)) {
                        optimizationsDisabled |= ((long)1 << i);
                    }
                }

                if (containsG) {
                    if (optimizationsDisabled != 0)
                    {
                        Console.WriteLine("Found run with -G- and some optimizations disabled: " + commandLine);
                        return false;
                    }
                    optimizationsDisabled = -1;
                } // containsG

                runs[optimizationsDisabled].Add(runTime);

                return true;

            } // AddRun

            public string getNormalized(long sortKey)   // returns a string so we can return the empty string if there's missing data
            {
                if (runs[sortKey].Count == 0 || runs[0].Count == 0 || runs[0].Average() == 0)
                {
                    return "";
                }

                return Convert.ToString((double)runs[sortKey].Average() / runs[0].Average());
            } // getNormalized

            public string getRunTime(long sortKey) // returns a string so we can return the empty string if there's missing data
            {
                if (runs[sortKey].Count == 0)
                {
                    return "";
                }

                return Convert.ToString(runs[sortKey].Average());
            } // getRunTime

            public int getNRuns(long sortKey)
            {
                return runs[sortKey].Count();
            } // getNRuns

            public string getNormalizedErrorBar(long sortKey) // returns a string so we can return the empty string if there's missing data
            {
                if (runs[sortKey].Count() < 1 || runs[0].Count() < 1)  // Change the constant to have a minimum number to compute standard error
                {
                    return "";  // fill this in
                }

                var meanAndStdDev = new ASETools.RunningMeanAndStdDev();
                runs[sortKey].ForEach(_ => meanAndStdDev.addValue(_ / runs[sortKey].Average()));

                return Convert.ToString(meanAndStdDev.getStandardError());

            } // getNormalizedErrorBar

            readonly List<Optimization> optimizations;
            public readonly string name;
            public readonly string shortName;
            Dictionary<long, List<int>> runs = new Dictionary<long, List<int>>();   // Maps optimizations disabled->set of runs.  Key -1 is for -G-



        } // Sample
        static void Main(string[] args)
        {
            var optimizations = new List<Optimization>();

            optimizations.Add(new Optimization("U", "Ukkonen", "-nu"));
            optimizations.Add(new Optimization("O", "Ordering", "-no"));
            optimizations.Add(new Optimization("T", "Truncation", "-nt"));
            optimizations.Add(new Optimization("B", "Banded Affine", "-nb"));
            optimizations.Add(new Optimization("E", "Edit distance", "-ne"));
            optimizations.Add(new Optimization("I", "Indel MaxK", "-ni"));

            var allOptimizationSets = new List<OptimizationSet>();

            //
            // Build the set of all possible combinations of optimizations.  There are 2^n+1 of these, and
            // we get them all by counting from 0..2^n - 1 and looking at the bits in the resulting numbers,
            // as well as adding in the special -1 set, which is -G- (no affine gap).
            //
            allOptimizationSets.Add(new OptimizationSet(optimizations));

            for (int i = 0; i < Math.Pow(2, optimizations.Count); i++)
            {
                var optimizationsDisabled = new List<bool>();
                for (int j = 0; j < optimizations.Count; j++)
                {
                    optimizationsDisabled.Add((i & (1 << j)) != 0);
                } // all optimizations

                allOptimizationSets.Add(new OptimizationSet(optimizations, optimizationsDisabled));
            } // all sets


            var inputFilename = @"f:\temp\optimizations_raw.txt";
            var inputFile = ASETools.CreateStreamReaderWithRetry(inputFilename);

            if (inputFile == null)
            {
                Console.WriteLine("Unable to open " + inputFilename);
                return;
            }

            string lastCommandLine = null;
            bool justSawStatusHeader = false;

            var allSamples = new List<Sample>();
            //
            // This is in descending order of read length.
            //
            allSamples.Add(new Sample("hg002", "hg2", optimizations, allOptimizationSets));
            allSamples.Add(new Sample("hg003", "hg3", optimizations, allOptimizationSets));
            allSamples.Add(new Sample("hg004", "hg4", optimizations, allOptimizationSets));
            allSamples.Add(new Sample("hg005", "hg5", optimizations, allOptimizationSets));
            allSamples.Add(new Sample("hg001", "hg1", optimizations, allOptimizationSets));
            allSamples.Add(new Sample("hg006", "hg6", optimizations, allOptimizationSets));
            allSamples.Add(new Sample("hg007", "hg7", optimizations, allOptimizationSets));
            allSamples.Add(new Sample("ERR194146", "146", optimizations, allOptimizationSets));
            allSamples.Add(new Sample("ERR194147", "147", optimizations, allOptimizationSets));

            string line;
            while (null != (line = inputFile.ReadLine()))
            {
                if (line.Contains("snap "))
                {
                    if (justSawStatusHeader)
                    {
                        Console.WriteLine("Saw command line right after status header.");
                        //return;
                    }

                    if (lastCommandLine != null)
                    {
                        Console.WriteLine("Saw two command lines in a row.");
                        //return;
                    }

                    lastCommandLine = line;
                    justSawStatusHeader = false;
                } else if (line.StartsWith("Total Reads"))
                {
                    if (justSawStatusHeader)
                    {
                        Console.WriteLine("Saw two status headers in a row");
                        //return;
                    }

                    justSawStatusHeader = true;
                } else if (justSawStatusHeader)
                {
                    bool anyWorked = false;
                    foreach (var sample in allSamples)
                    {
                        bool worked = sample.AddRun(lastCommandLine, line);
                        if (worked && anyWorked)
                        {
                            Console.WriteLine("Two samples claimed the same run");
                            return;
                        }
                        anyWorked |= worked;
                    }

                    if (!anyWorked)
                    {
                        Console.WriteLine("No sample claimed a run");
                        return;
                    }

                    justSawStatusHeader = false;
                    lastCommandLine = null;
                } // Otherwise, it's an unintersting line.
            } // foreach input line

            //
            // Now generate output.
            //

            var outputFile = ASETools.CreateStreamWriterWithRetry(@"f:\temp\EffectOfDisablingOptimizations.txt");
            if (null == outputFile)
            {
                return;
            }

            outputFile.Write("Optimization Disabled");
            foreach (var sample in allSamples)
            {
                outputFile.Write("\t" + sample.shortName);   // this is normalized, but leave it without to be the legend key in the graph in excel
            }

            foreach (var sample in allSamples)
            {
                outputFile.Write("\t" + sample.shortName + " runtime");
            }

            foreach (var sample in allSamples)
            {
                outputFile.Write("\t" + sample.shortName + " n runs");
            }

            foreach (var sample in allSamples)
            {
                outputFile.Write("\t" + sample.shortName + " normalized error bar");
            }

            outputFile.WriteLine();

            foreach (var optimizationSet in allOptimizationSets)
            {
                var sortKey = optimizationSet.getSortKey();

                outputFile.Write(ASETools.ConvertToExcelString(optimizationSet.getFlagsAsString()));

                foreach (var sample in allSamples)
                {
                    outputFile.Write("\t" + sample.getNormalized(sortKey));
                }

                foreach (var sample in allSamples)
                {
                    outputFile.Write("\t" + sample.getRunTime(sortKey));
                }

                foreach (var sample in allSamples)
                {
                    outputFile.Write("\t" + sample.getNRuns(sortKey));
                }

                foreach (var sample in allSamples)
                {
                    outputFile.Write("\t" + sample.getNormalizedErrorBar(sortKey));
                }

                outputFile.WriteLine();
            } // foreach optimization set

            outputFile.WriteLine("**done**");
            outputFile.Close();

        } // Main
    } // Program
} // namespace
