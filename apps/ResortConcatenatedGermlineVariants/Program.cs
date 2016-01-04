using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace ResortConcatenatedGermlineVariants
{
    class Program
    {
        static void Main(string[] args)
        {
            if (args.Count() < 2)
            {
                Console.WriteLine("usage: ResortConcatenatedGermlineVariants inputFile outputFile {list of genes to extract to their own files}");
                return;
            }

            StreamReader inputFile = new StreamReader(args[0]);
            StreamWriter outputFile = new StreamWriter(args[1]);

            var variantsByGene = new Dictionary<string, List<string>>();

            int totalVariants = 0;

            string line;
            while (null != (line = inputFile.ReadLine()))
            {
                var fields = line.Split('\t');
                if (fields.Count() != 2)
                {
                    Console.WriteLine("Expected header line, got " + line);
                    return;
                }

                string gene = fields[0];
                int nVariants = Convert.ToInt32(fields[1]);

                if (!variantsByGene.ContainsKey(gene))
                {
                    variantsByGene.Add(gene, new List<string>());
                }

                for (int i = 0; i < nVariants; i++)
                {
                    variantsByGene[gene].Add(inputFile.ReadLine());
                }

                totalVariants++;
            }

            foreach (var entry in variantsByGene)
            {
                outputFile.WriteLine("" + entry.Key + "\t" + entry.Value.Count());
                foreach (var variant in entry.Value)
                {
                    outputFile.WriteLine(variant);
                }
            }

            outputFile.Close();

            for (int i = 2; i < args.Count(); i++)
            {
                string gene = args[i];
                if (!variantsByGene.ContainsKey(gene))
                {
                    Console.WriteLine("Found no recorded variants for gene " + gene);
                }
                else
                {
                    outputFile = new StreamWriter("germline_variants_for_" + gene + ".txt");
                    outputFile.WriteLine(gene + "\t" + variantsByGene[gene].Count());

                    foreach (var variant in variantsByGene[gene])
                    {
                        outputFile.WriteLine(variant);
                    }

                    outputFile.Close();
                }
            }

            Console.WriteLine("" + variantsByGene.Count() + " genes containing " + totalVariants + " variants");
        }
    }
}
