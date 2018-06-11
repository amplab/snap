using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;


namespace QNDGetP53GeneHancers
{
    class Program
    {
        static void Main(string[] args)
        {
            var configuration = ASETools.Configuration.loadFromFile(args);
            if (null == configuration)
            {
                Console.WriteLine("Unable to load configuration.");
                return;
            }

            var geneHancersByGene = ASETools.GeneHancerLine.ReadFromFileToDict(configuration.geneHancerFilename);

            var outputFile = ASETools.CreateStreamWriterWithRetry("tp53_geneHancers.txt");
            outputFile.WriteLine("chrom\tstart\tend\tgeneHancer ID\tfeature\toverall geneHancer score\tTP53 specific score");

            foreach (var geneHancer in geneHancersByGene["TP53"])
            {
                outputFile.WriteLine(geneHancer.chromosome + "\t" + geneHancer.start + "\t" + geneHancer.end + "\t" + geneHancer.geneHancerId + "\t" + geneHancer.feature + "\t" + geneHancer.score + "\t" + geneHancer.geneScores["TP53"]);
            }

            outputFile.Close();
        }
    }
}
