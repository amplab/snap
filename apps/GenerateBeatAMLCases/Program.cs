using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Diagnostics;
using ASELib;

namespace GenerateBeatAMLCases
{
    class Program
    {

        class Patient
        {
            public readonly int patientId;
            public Dictionary<bool, ASETools.BeatAMLSamplesSummaryLine> samples = new Dictionary<bool, ASETools.BeatAMLSamplesSummaryLine>();   // true->tumor, false->normal

            public Patient(int patientId_)
            {
                patientId = patientId_;
            }
        }

        static ASETools.Configuration configuration;

        static long GetFileSize(string sampleId, bool rna)
        {
            var fileNames = Directory.EnumerateFiles(configuration.synapseDirectory + (rna ? @"rnaseq\bam\": @"seqcap\bam\"), (rna ? (@"sample_" + sampleId + "*.bam") : ("Sample_" + sampleId + "*.realign.bam"))).ToList();

            if (fileNames.Count() > 1)
            {
                Console.WriteLine("Got wrong count of files for sample id " + sampleId + ": " + fileNames);
                return 0;
            } else if (fileNames.Count() == 0)
            {
                Console.WriteLine("Missing file for sample ID " + sampleId + ", rna: " + rna);
                return 0;
            }

            return new FileInfo(fileNames[0]).Length;
        }

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            configuration = ASETools.Configuration.loadFromFile(args);
            if (null == configuration)
            {
                Console.WriteLine("Unable to load configuration.");
                return;
            }

            if (!configuration.isBeatAML)
            {
                Console.WriteLine("Configuration is not for BeatAML.  Did you forget to add the -configuration option?");
                return;
            }

            var samplesSummaryLines = ASETools.BeatAMLSamplesSummaryLine.ReadFile(configuration.samplesSummaryPathname);
            if (null == samplesSummaryLines)
            {

                Console.WriteLine("Unable to read samples summary lines.");
            }

            var patients = new Dictionary<int, Patient>();

            foreach (var sample in samplesSummaryLines)
            {
                var patientId = sample.patient_id;
                bool tumor = !sample.specimen_type.Contains("Skin Biopsy"); // Using .Contains because of "Skin Biopsy Diagnosis" and "Skin Biopsy Recovery", each of which occur once.

                if (!patients.ContainsKey(patientId))
                {
                    patients.Add(patientId, new Patient(patientId));
                }

                var patient = patients[patientId];

                if (!patient.samples.ContainsKey(tumor))
                {
                    patient.samples.Add(tumor, sample);
                } else
                {
                    if (sample.specimen_has_rna_seq && !patient.samples[tumor].specimen_has_rna_seq)
                    {
                        patient.samples[tumor] = sample;
                    } else
                    {
                        // Should probably decide which one to take here, but for now we'll just keep the first one.
                    }
                }
            } // foreach sample

            var patientsPutativelyWithSufficientData = patients.Select(x => x.Value).Where(x => x.samples.ContainsKey(true) && x.samples.ContainsKey(false) && x.samples[true].specimen_has_rna_seq);

            Console.WriteLine("Total of " + patients.Count() + " patients, of which " + patientsPutativelyWithSufficientData.Count() + " putatively have sufficient data.");

            int nNoMaf = 0;
            int nNoTumorDNA = 0;
            int nNoTumorRNA = 0;
            int nNoNormalDNA = 0;
            int nNoNormalRNA = 0;

            var goodCases = new Dictionary<string,ASETools.Case>();

            foreach (var patient in patientsPutativelyWithSufficientData)
            {
                var mafFilename = configuration.synapseDirectory + @"seqcap\maf\pair_" + patient.samples[true].seq_id + "_aml_" + patient.samples[false].seq_id + "_skin.mutect.somatic.filter.maf";
                if (!File.Exists(mafFilename))
                {
                    nNoMaf++;
                    continue;
                }

                var case_ = new ASETools.Case();

                case_.case_id = ASETools.PadStringToGuidLength("Case-" + patient.patientId);

                case_.maf_filename = mafFilename;   // In GDC, there is one MAF per disease type.  Here, it's one per patient.  The pipeline does a maf extraction step, which for BeatAML will just copy and reformat the MAF into a data directory.
                case_.maf_file_id = patient.patientId + "-maf";

                var tumorDNAFiles = Directory.EnumerateFiles(configuration.synapseDirectory + @"seqcap\bam\", "Sample_" + patient.samples[true].seq_id + "*.realign.bam").ToList();
                if (tumorDNAFiles.Count() == 0)
                {
                    nNoTumorDNA++;
                    continue;
                }


                var tumorRNAFiles = Directory.EnumerateFiles(configuration.synapseDirectory + @"rnaseq\bam\", @"Sample_" + patient.samples[true].seq_id + "*.bam").ToList();
                if (tumorRNAFiles.Count() == 0)
                {
                    nNoTumorRNA++;
                    continue;
                }

                var normalDNAFiles = Directory.EnumerateFiles(configuration.synapseDirectory + @"seqcap\bam\", "Sample_" + patient.samples[false].seq_id + "*.realign.bam").ToList();
                if (normalDNAFiles.Count() == 0)
                {
                    nNoNormalDNA++;
                    continue;
                }


                var normalRNAFiles = Directory.EnumerateFiles(configuration.synapseDirectory + @"rnaseq\bam\", @"sample_" + patient.samples[false].seq_id + "*.bam").ToList();
                if (normalRNAFiles.Count() == 0)
                {
                    nNoNormalRNA++;
                    // Don't skip this, we don't need it.
                }

                //
                // There are no GUIDs for these, so generate them.
                //
                case_.tumor_dna_file_id = ASETools.PadStringToGuidLength(patient.samples[true].seq_id + "-dna");
                case_.normal_dna_file_id = ASETools.PadStringToGuidLength(patient.samples[false].seq_id + "-dna");
                case_.tumor_rna_file_id = ASETools.PadStringToGuidLength(patient.samples[true].seq_id + "-rna");
                if (normalRNAFiles.Count() > 0)
                {
                    case_.normal_rna_file_id = ASETools.PadStringToGuidLength(patient.samples[false].seq_id + "-rna");
                    case_.normal_rna_size = new FileInfo(normalRNAFiles[0]).Length;
                }

                case_.tumor_dna_size = new FileInfo(tumorDNAFiles[0]).Length;
                case_.normal_dna_size = new FileInfo(normalDNAFiles[0]).Length;
                case_.tumor_rna_size = new FileInfo(tumorRNAFiles[0]).Length;

                case_.project_id = "BeatAML-LAML";

                case_.sample_ids = new List<string>();
                case_.sample_ids.Add(patient.samples[true].seq_id);
                case_.sample_ids.Add(patient.samples[false].seq_id);
                case_.maf_filename = mafFilename;

                goodCases.Add(case_.case_id,case_);
            }

            Console.WriteLine(goodCases.Count() + " total cases actually have data.  " + nNoMaf + " don't have a MAF, " + nNoTumorDNA + " don't have tumor DNA, " + nNoTumorRNA + " don't have tumor RNA and " + nNoNormalDNA + " don't have normal DNA.");

            ASETools.Case.SaveToFile(goodCases, configuration.casesFilePathname);
            Console.WriteLine(ASETools.ElapsedTimeInSeconds(timer));

        }
    }
}
