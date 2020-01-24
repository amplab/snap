using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.InteropServices;
using System.Text;
using System.Threading.Tasks;
using ASELib;

namespace AnalyzeMiRNAExpression
{
    class Program
    {
        static ASETools.CommonData commonData;
        static List<string> miRNANames;
        static List<string> geneNames;
        static PerThreadState finalResult;  // Don't initialize this here; we need to set up some globals that are used by the constructor first.

        static void Main(string[] args)
        {
            commonData = ASETools.CommonData.LoadCommonData(args);

            if (null == commonData)
            {
                return;
            }

            if (commonData.listOfCases.Any(_ => _.tumor_miRNA_expression_quantification_filename == "" && _.tumor_miRNA_expression_quantification_file_id != "" || _.annotated_selected_variants_filename == ""))
            {
                Console.WriteLine("At least one case that has tumor miRNA quantification is missing an input.  Download it and try again.");
                //BJB debug return;
            }

            var casesToRun = commonData.listOfCases.Where(_ => _.tumor_miRNA_expression_quantification_filename != "" && _.annotated_selected_variants_filename != "").ToList();
            if (casesToRun.Count() == 0)
            {
                Console.WriteLine("No cases have data (?)");
                return;
            }

            //
            // We need to initialize the list of miRNANames.  Read in a single case's miRNA quantification.
            //
            var miRNAExpressionQuantification = ASETools.miRNAExpressionQuantification.LoadFromFile(casesToRun[0].tumor_miRNA_expression_quantification_filename);
            miRNANames = miRNAExpressionQuantification.Select(_ => _.Key).ToList();
            geneNames = ASETools.SelectedGene.LoadFromFile(commonData.configuration.selectedGenesFilename).Select(_ => _.Hugo_Symbol).ToList();

            finalResult = new PerThreadState(); // Don't initialize this before miRNANames and geneNames.

            int nPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "cases", casesToRun.Count(), out nPerDot);

            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, PerThreadState>(casesToRun, HandleOneCase, FinishUp, null, nPerDot);

            threading.run();
            Console.WriteLine(ASETools.ElapsedTimeInSeconds(commonData.timer));





        } // Main


        class PerThreadState
        {
            public Dictionary<string, Dictionary<string, Dictionary<int, List<double>>>> miRNAExpressionByGeneAndMutationCount = new Dictionary<string, Dictionary<string, Dictionary<int, List<double>>>>();    // miRNA -> hugo_symbol -> 0,1,2 -> list of reads/million

            public PerThreadState()
            {
                foreach (var miRNA in miRNANames)
                {
                    miRNAExpressionByGeneAndMutationCount.Add(miRNA, new Dictionary<string, Dictionary<int, List<double>>>());
                    foreach (var geneName in geneNames)
                    {
                        miRNAExpressionByGeneAndMutationCount[miRNA].Add(geneName, new Dictionary<int, List<double>>());
                        foreach (var nMutations in ASETools.ZeroOneTwo)
                        {
                            miRNAExpressionByGeneAndMutationCount[miRNA][geneName].Add(nMutations, new List<double>());
                        } // mutation count
                    } // gene
                } // miRNA
            } // ctor

            public void merge(PerThreadState peer)
            {
                foreach (var miRNA in miRNANames)
                {
                    foreach (var geneName in geneNames)
                    {
                        foreach (var nMutations in ASETools.ZeroOneTwo)
                        {
                            miRNAExpressionByGeneAndMutationCount[miRNA][geneName][nMutations].AddRange(peer.miRNAExpressionByGeneAndMutationCount[miRNA][geneName][nMutations]);
                        } // mutation count
                    } // gene
                } // miRNA
            }// merge
        } // PerThreadState

        static void FinishUp(PerThreadState state)
        {
            lock (finalResult)
            {
                finalResult.merge(state);
            } // lock
        } // FinishUp

        static void HandleOneCase(ASETools.Case case_, PerThreadState state)
        {
            var miRNAExpressionQuantification = ASETools.miRNAExpressionQuantification.LoadFromFile(case_.tumor_miRNA_expression_quantification_filename);
            var asv = ASETools.AnnotatedVariant.readFile(case_.annotated_selected_variants_filename);

            foreach (var geneName in geneNames)
            {
                var mutationIndex = ASETools.ZeroOneMany(asv.Where(_ => _.Hugo_symbol == geneName && _.somaticMutation && !_.isSilent()).Count());

                foreach (var miRNAName in miRNAExpressionQuantification.Select(_ => _.Key).ToList())
                {
                    state.miRNAExpressionByGeneAndMutationCount[miRNAName][geneName][mutationIndex].Add(miRNAExpressionQuantification[miRNAName].reads_mer_milion_miRNA_mapped);
                } // miRNA
            } // gene
        } // HandleOneCase


    }
}
