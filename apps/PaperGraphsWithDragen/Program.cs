using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;
using System.Diagnostics;

namespace PaperGraphsWithDragen
{
    class Program
    {
        enum Aligner { SNAP, Bowtie, BWA, DRAGEN, GEM, Novoalign};

        static Dictionary<Aligner, string> alignerDesignation = new Dictionary<Aligner, string>();  // This is the string used in the filenames, i.e. "bowtie_s1000".
        static Dictionary<Aligner, string> alignerName = new Dictionary<Aligner, string>(); // The printed name of the aligner for the output, i.e., "bowtie".

        static void SetAlignerStrings()
        {
            alignerDesignation.Add(Aligner.SNAP, "snap");
            alignerDesignation.Add(Aligner.Bowtie, "bowtie_s1000");
            alignerDesignation.Add(Aligner.BWA, "bwa");
            alignerDesignation.Add(Aligner.DRAGEN, "dragen3.8");
            alignerDesignation.Add(Aligner.GEM, "gem");
            alignerDesignation.Add(Aligner.Novoalign, "novoalign");

            alignerName.Add(Aligner.SNAP, "SNAP");
            alignerName.Add(Aligner.Bowtie, "Bowtie");
            alignerName.Add(Aligner.BWA, "BWA");
            alignerName.Add(Aligner.DRAGEN, "DRAGEN");
            alignerName.Add(Aligner.GEM, "GEM");
            alignerName.Add(Aligner.Novoalign, "Novoalign");

        } // SetAlignerStrings

        class Sample
        {
            public readonly string name;
            public Dictionary<Aligner, ASETools.ConcordanceResults> concordance = new Dictionary<Aligner, ASETools.ConcordanceResults>();
            public Dictionary<Aligner, ASETools.RunTiming> timings = new Dictionary<Aligner, ASETools.RunTiming>();

            public Sample(string name_)
            {
                name = name_;
            } // ctor


        } // Sample
        static void Main(string[] args)
        {
            SetAlignerStrings();

            string concordanceDirectory = @"d:\temp\concordance-dragen3.8VC\";
            string timingDirectory = @"d:\temp\timings\";

            var samples = new Dictionary<string, Sample>();
            for (int i = 1; i <= 7; i++)
            {
                samples.Add("hg00" + i, new Sample("hg00" + i));
            }

            samples.Add("ERR194146", new Sample("ERR194146"));
            samples.Add("ERR194147", new Sample("ERR194147"));

            string oneTabPerAligner = "";
            foreach (var aligner in ASETools.EnumUtil.GetValues<Aligner>()) oneTabPerAligner += "\t";

            string tabSeparatedAlignerNamesWithLeadingTab = "";
            foreach (var aligner in ASETools.EnumUtil.GetValues<Aligner>())
            {
                tabSeparatedAlignerNamesWithLeadingTab += "\t" + alignerName[aligner];
            }

            foreach (var sample in samples.Select(_ => _.Value))
            {
                foreach (var aligner in ASETools.EnumUtil.GetValues<Aligner>())
                {
                    var concordanceFilename = concordanceDirectory + sample.name + "." + alignerDesignation[aligner] + ".dragen3.8VC.concordance.tar";
                    if (File.Exists(concordanceFilename))
                    {
                        sample.concordance.Add(aligner, new ASETools.ConcordanceResults(concordanceFilename));
                    }

                    var timingFilename = timingDirectory + sample.name + "." + (aligner == Aligner.SNAP ? "snap" : alignerDesignation[aligner]) + "_timings.txt";
                    if (File.Exists(timingFilename))
                    {
                        if (aligner == Aligner.SNAP)
                        {
                            sample.timings.Add(aligner, ASETools.RunTiming.LoadFromSNAPFile(timingFilename));
                        }
                        else
                        {
                            sample.timings.Add(aligner, ASETools.RunTiming.LoadFromLinuxFile(timingFilename));
                        }
                    } else
                    {
                        var splitTimingsFilenames = new List<string>();
                        for (int i = 0; i < 9; i++)
                        {
                            splitTimingsFilenames.Add(timingDirectory + sample.name + "." + i + "." + alignerDesignation[aligner] + "_timings.txt");
                        }

                        var secondStageFilename = timingDirectory + sample.name + "." + alignerDesignation[aligner] + "_second_stage_timings.txt";

                        if (File.Exists(secondStageFilename) && splitTimingsFilenames.TrueForAll(_ => File.Exists(_)))
                        {
                            sample.timings.Add(aligner, ASETools.RunTiming.LoadFromSplitFiles(splitTimingsFilenames, secondStageFilename));
                        }
                    } // case for split timings
                } // aligner
            } // sample

            var outputFile = ASETools.CreateStreamWriterWithRetry(@"d:\temp\Dragen3.8VCConcordanceAndTimings.txt");

            outputFile.WriteLine("\tMean F1" + oneTabPerAligner + "\tIndel F1" +  oneTabPerAligner + "\tSNV F1" + oneTabPerAligner+ "\tIndel Recall" + oneTabPerAligner + "\tIndel Precision" + oneTabPerAligner + "\tSNV Recall" + oneTabPerAligner + "\tSNV Precision");
            for (int i = 0; i < 7; i++)
            {
                outputFile.Write(tabSeparatedAlignerNamesWithLeadingTab + "\t");  // The extra tab at the end is for the next header
            }
            outputFile.WriteLine();

            foreach (var sample in samples.Select(_ => _.Value))
            {
                outputFile.Write(sample.name);

                foreach (var aligner in ASETools.EnumUtil.GetValues<Aligner>())
                {
                    outputFile.Write("\t");
                    if (sample.concordance.ContainsKey(aligner))
                    {
                        outputFile.Write(sample.concordance[aligner].Mean_F1_Score());
                    }
                } // aligner

                outputFile.Write("\t" + sample.name);

                foreach (var aligner in ASETools.EnumUtil.GetValues<Aligner>())
                {
                    outputFile.Write("\t");
                    if (sample.concordance.ContainsKey(aligner))
                    {
                        outputFile.Write(sample.concordance[aligner].results[ASETools.VariantType.Indel].F1_score);
                    }
                } // aligner

                outputFile.Write("\t" + sample.name);

                foreach (var aligner in ASETools.EnumUtil.GetValues<Aligner>())
                {
                    outputFile.Write("\t");
                    if (sample.concordance.ContainsKey(aligner))
                    {
                        outputFile.Write(sample.concordance[aligner].results[ASETools.VariantType.SNV].F1_score);
                    }
                } // aligner

                outputFile.Write("\t" + sample.name);

                foreach (var aligner in ASETools.EnumUtil.GetValues<Aligner>())
                {
                    outputFile.Write("\t");
                    if (sample.concordance.ContainsKey(aligner))
                    {
                        outputFile.Write(sample.concordance[aligner].results[ASETools.VariantType.Indel].recall);
                    }
                } // aligner

                outputFile.Write("\t" + sample.name);

                foreach (var aligner in ASETools.EnumUtil.GetValues<Aligner>())
                {
                    outputFile.Write("\t");
                    if (sample.concordance.ContainsKey(aligner))
                    {
                        outputFile.Write(sample.concordance[aligner].results[ASETools.VariantType.Indel].precision);
                    }
                } // aligner

                outputFile.Write("\t" + sample.name);

                foreach (var aligner in ASETools.EnumUtil.GetValues<Aligner>())
                {
                    outputFile.Write("\t");
                    if (sample.concordance.ContainsKey(aligner))
                    {
                        outputFile.Write(sample.concordance[aligner].results[ASETools.VariantType.SNV].recall);
                    }
                } // aligner

                outputFile.Write("\t" + sample.name);

                foreach (var aligner in ASETools.EnumUtil.GetValues<Aligner>())
                {
                    outputFile.Write("\t");
                    if (sample.concordance.ContainsKey(aligner))
                    {
                        outputFile.Write(sample.concordance[aligner].results[ASETools.VariantType.SNV].precision);
                    }
                } // aligner

                outputFile.WriteLine();
            } // sample

            outputFile.WriteLine();

            outputFile.WriteLine("\tLoad & Align" + oneTabPerAligner + "\tLoad & Align Error Bar" + oneTabPerAligner + "\tWhole pipeline");

            foreach (var sample in samples.Select(_ => _.Value))
            {
                outputFile.Write(sample.name);

                foreach (var aligner in ASETools.EnumUtil.GetValues<Aligner>())
                {

                }
            }

            outputFile.Close();

        } // Main
    } // Program
} // namespace
