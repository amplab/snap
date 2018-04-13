using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;

namespace ASEStats
{
    class Program
    {

        static void DescribeGeneASE(Dictionary<string, ASETools.ASEMapPerGeneLine> perGeneASE, string hugo_symbol)
        {
            if (!perGeneASE.ContainsKey(hugo_symbol))
            {
                Console.WriteLine("Gene " + hugo_symbol + " is not in the per-gene ASE list.");
                return;
            }

            var x = perGeneASE[hugo_symbol];
            Console.WriteLine(hugo_symbol + " has normal ASE of " + x.sampleData[false].meanASE + " with " + x.sampleData[false].nSamples + " samples.");
        }
        static void Main(string[] args)
        {
            var configuration = ASETools.Configuration.loadFromFile(args);

            if (null == configuration)
            {
                Console.WriteLine("Giving up because we were unable to load configuration.");
                return;
            }

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname).Select(x => x.Value).ToList();

            if (null == cases)
            {
                Console.WriteLine("Unable to load cases.");
            }

            Console.WriteLine("Total of " + cases.Count() + " cases.");
            Console.WriteLine(cases.Where(x => x.normal_rna_file_id != "").Count() + " have normal RNA");
            Console.WriteLine(cases.Where(x => x.normal_methylation_file_id != "").Count() + " have methylation");

            var byDisease = cases.GroupByToDict(x => x.disease());

            Console.WriteLine();
            Console.WriteLine("count\tnormal\tcopy\tdisease");
            
            foreach (var disease in byDisease)
            {
                Console.WriteLine(disease.Value.Count() + "\t" + disease.Value.Where(x => x.normal_rna_file_id != "").Count() + "\t" +  disease.Value.Where(x => x.tumor_copy_number_filename != "").Count() + "\t" + disease.Key);
            }

            Console.WriteLine("----    ----    ----");
            Console.WriteLine(byDisease.Select(x => x.Value.Count()).Sum() + "\t" + byDisease.Select(x => x.Value.Where(y => y.normal_rna_file_id != "").Count()).Sum() + "\t" + byDisease.Select(x => x.Value.Where(y => y.tumor_copy_number_filename != "").Count()).Sum());

            Console.WriteLine();

            var perGeneASEMap = ASETools.ASEMapPerGeneLine.ReadFromFileToDictionary(configuration.finalResultsDirectory + ASETools.PerGeneASEMapFilename);

            var tooHighASEGenes = perGeneASEMap.Select(x => x.Value).Where(x => x.sampleData[false].meanASE >= configuration.ASEInNormalAtWhichToExcludeGenes).Select(x => x.hugo_symbol).ToList();
            tooHighASEGenes.Sort();

            Console.Write("The following " + tooHighASEGenes.Count() + " genes are excluded as ASE measurement sites due to too high normal ASE: ");
            tooHighASEGenes.ForEach(x => Console.Write(x + ", "));

            Console.WriteLine();
            Console.WriteLine();

            var highishASEGenes = perGeneASEMap.Select(x => x.Value).Where(x => x.sampleData[false].meanASE >= 0.33 && x.sampleData[false].meanASE < configuration.ASEInNormalAtWhichToExcludeGenes).Select(x => x.hugo_symbol).ToList();
            highishASEGenes.Sort();

            Console.Write("The following additional " + highishASEGenes.Count() + " genes would be excluded ASE measurement sites due to too high normal ASE if we used a cutoff of 0.33: ");
            highishASEGenes.ForEach(x => Console.Write(x + ", "));

            Console.WriteLine();

            DescribeGeneASE(perGeneASEMap, "SNRPN");
            DescribeGeneASE(perGeneASEMap, "NDN");
            DescribeGeneASE(perGeneASEMap, "UBE3A");
            DescribeGeneASE(perGeneASEMap, "DIRAS3");
            DescribeGeneASE(perGeneASEMap, "CLPG");

        } // Main
    }
}
