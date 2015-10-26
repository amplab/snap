using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace AllExpression
{
    class Program
    {
        class Experiment
        {
            public static Experiment parse(string line)
            {
                Experiment experiment = new Experiment();

                var fields = line.Split('\t');

                experiment.disease_abbr = fields[0];
                experiment.reference = fields[1];
                if (experiment.disease_abbr == "ov" && experiment.reference == "ncbi36_bccagsc_variant")
                {
                    experiment.disease_abbr = "ov_hg18";    // Special case, since this is aligned against two different references
                }
                experiment.participantID = fields[2];
                experiment.TumorRNAAnalysis = fields[3];
                experiment.TumorDNAAnalysis = fields[4];
                experiment.NormalDNAAnalysis = fields[5];
                experiment.NormalRNAAnalysis = fields[6];
                experiment.TumorRNAPathname = fields[7];
                experiment.TumorDNAPathname = fields[8];
                experiment.NormalDNAPathname = fields[9];
                experiment.NormalRNAPathname = fields[10];
                experiment.VCFPathname = fields[11];
                experiment.gender = fields[12];
                experiment.daysToBirth = fields[13];
                experiment.daysToDeath = fields[14];
                experiment.OrigTumorDNAAliquotID = fields[15];
                experiment.TumorAllcountFile = fields[16];
                experiment.NormalAllcountFile = fields[17];

                return experiment;
            }

            public string disease_abbr;
            public string reference;
            public string participantID;
            public string TumorRNAAnalysis;
            public string TumorDNAAnalysis;
            public string NormalDNAAnalysis;
            public string NormalRNAAnalysis;
            public string TumorRNAPathname;
            public string TumorDNAPathname;
            public string NormalDNAPathname;
            public string NormalRNAPathname;
            public string VCFPathname;
            public string gender;
            public string daysToBirth;
            public string daysToDeath;
            public string OrigTumorDNAAliquotID;
            public string TumorAllcountFile;
            public string NormalAllcountFile;
        }

        class ReadCounts
        {
            public static ReadCounts parse(string line) 
            {
                ReadCounts readCounts = new ReadCounts();

                var fields = line.Split('\t');
                readCounts.analysisID = fields[0].ToLower();

                var subfields = fields[1].Split(' ');
                readCounts.hqNuclearReads = Convert.ToInt64(subfields[0]);
                readCounts.lqReads = Convert.ToInt64(subfields[6]);
                readCounts.mitochondrialReads = Convert.ToInt64(subfields[10]);
                readCounts.totalReads = Convert.ToInt64(subfields[13]);

                return readCounts;
            }

            public string analysisID;

            public long hqNuclearReads;
            public long lqReads;
            public long mitochondrialReads;
            public long totalReads;
        }

        class NMeanAndStdDev
        {
            public static NMeanAndStdDev parse(string line)
            {
                string[] fields = line.Split('\t');

                NMeanAndStdDev nMeanAndStdDev = new NMeanAndStdDev();

                nMeanAndStdDev.disease_abbr = fields[0].Substring(11);
                nMeanAndStdDev.n = Convert.ToInt64(fields[2]);
                nMeanAndStdDev.mean = Convert.ToDouble(fields[3]);
                nMeanAndStdDev.stdDev = Convert.ToDouble(fields[4]);

                return nMeanAndStdDev;
            }

            public string disease_abbr;
            public long n;
            public double mean;
            public double stdDev;
        }
        static void Main(string[] args)
        {
            var experiments = new List<Experiment>();

            StreamReader experimentReader = new StreamReader(@"f:\temp\expression\experiments.txt");

            experimentReader.ReadLine();    // Skip header
            string line;
            while (null != (line = experimentReader.ReadLine()))
            {
                Experiment experiment = Experiment.parse(line);
                experiments.Add(experiment);
            }

            experimentReader.Close();

            //
            // Get the read counts for all tumor RNA samples.
            //
            var readCountsByAnalysisID = new Dictionary<string, ReadCounts>();

            StreamReader readCountReader = new StreamReader(@"f:\sequence\reads\tcga\read_counts.txt");

            while (null != (line = readCountReader.ReadLine()))
            {
                ReadCounts readCounts = ReadCounts.parse(line);
                readCountsByAnalysisID.Add(readCounts.analysisID, readCounts);
            }
            readCountReader.Close();

            //
            // And get the per-disease distributions for the tp53 locus.
            //
            var statsByDisease = new Dictionary<string, NMeanAndStdDev>();

            StreamReader statsReader = new StreamReader(@"f:\sequence\reads\tcga\p53_locus_per_disease_counts.txt");
            while (null != (line = statsReader.ReadLine()))
            {
                var stats = NMeanAndStdDev.parse(line);
                statsByDisease.Add(stats.disease_abbr, stats);
            }
            statsReader.Close();

            //
            // Get the per-experiment counts of reads matching the P53 locus.
            //
            var p53LocusReadCountsByAnalysisID = new Dictionary<string, double>();

            StreamReader p53Reader = new StreamReader(@"f:\sequence\reads\tcga\p53-counts.txt");
            while (null != (line = p53Reader.ReadLine())) {
                string[] fields = line.Split('\t');
                string [] subfields = fields[3].Split('.');

                p53LocusReadCountsByAnalysisID.Add(subfields[0], Convert.ToDouble(fields[0]));
            }

            p53Reader.Close();


            StreamWriter output = new StreamWriter(@"f:\temp\expression\p53_rates.txt");
            output.WriteLine("TumorRNAAnalysisID\tdisease_abbr\tfractionOfMean");

            foreach (var experiment in experiments)
            {
                var analysisID = experiment.TumorRNAAnalysis;
                double mappedReads = p53LocusReadCountsByAnalysisID[analysisID];
                double totalReads = (double)readCountsByAnalysisID[experiment.TumorRNAAnalysis].hqNuclearReads;
                double meanReads = statsByDisease[experiment.disease_abbr].mean;

                output.WriteLine(analysisID + "\t" + experiment.disease_abbr + "\t" + (mappedReads / totalReads) / meanReads);
            }

            output.Close();

        }
    }
}
