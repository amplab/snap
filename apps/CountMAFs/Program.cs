using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Threading;
using ASELib;

namespace CountMAFs
{
    class Program
    {
        static void Main(string[] args)
        {
            var configuration = ASETools.Configuration.loadFromFile(args);

            if (configuration.casesFilePathname == "")
            {
                Console.WriteLine("Must generate cases first.");
                return;
            }

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

            var casesToProcess = cases.Where(x => x.Value.maf_filename != "").Select(x => x.Value).ToList();

            Console.WriteLine("Processing " + casesToProcess.Count() + " of " + cases.Count() + " cases.");

            var mutationsByType = new Dictionary<string, int>();
            var mutationsByDisease = new Dictionary<string, int>();
            var casesByDisease = new Dictionary<string, int>();

            foreach (var case_ in casesToProcess)
            {
                if (!casesByDisease.ContainsKey(case_.disease()))
                {
                    casesByDisease.Add(case_.disease(), 0);
                }

                casesByDisease[case_.disease()]++;
            }

            var threads = new List<Thread>();
            for (int i = 0; i < Environment.ProcessorCount; i++)
            {
                threads.Add(new Thread(() => WorkerThread(casesToProcess, mutationsByType, mutationsByDisease)));
            }

            threads.ForEach(t => t.Start());
            threads.ForEach(t => t.Join());

            var nMutations = mutationsByDisease.Select(x => x.Value).Sum();

            Console.WriteLine("Disease\tnMutations\tfrac");
            foreach (var disease in mutationsByDisease)
            {
                Console.WriteLine(disease.Key + "\t" + disease.Value + "\t" + (double)disease.Value / nMutations);
            }

            Console.WriteLine();
            Console.WriteLine("Variant\tnMutations\tfrac");
            foreach (var variant in mutationsByType)
            {
                Console.WriteLine(variant.Key + "\t" + variant.Value + "\t" + (double)variant.Value / nMutations);
            }

            Console.WriteLine();

            Console.WriteLine("Total of " + nMutations + " mutations.");
        }

        static void WorkerThread(List<ASETools.Case> casesToProcess, Dictionary<string, int> mutationsByType, Dictionary<string, int> mutationsByDisease)
        {
            var localMutationsByType = new Dictionary<string, int>();
            var localMutationsByDisease = new Dictionary<string, int>();

            while (true)
            {
                ASETools.Case case_;
                lock (casesToProcess)
                {
                    if (casesToProcess.Count() == 0)
                    {
                        foreach (var mutationByType in localMutationsByType)
                        {
                            if (!mutationsByType.ContainsKey(mutationByType.Key)) {
                                mutationsByType.Add(mutationByType.Key, 0);
                            }
                            mutationsByType[mutationByType.Key] += mutationByType.Value;
                        }

                        foreach (var mutationByDisease in localMutationsByDisease)
                        {
                            if (!mutationsByDisease.ContainsKey(mutationByDisease.Key))
                            {
                                mutationsByDisease.Add(mutationByDisease.Key, 0);
                            }

                            mutationsByDisease[mutationByDisease.Key] += mutationByDisease.Value;
                        }

                        return;
                    }

                    case_ = casesToProcess[0];
                    casesToProcess.RemoveAt(0);
                }

                var mafLines = ASETools.MAFLine.ReadFile(case_.extracted_maf_lines_filename, case_.maf_file_id, false);
                var disease = case_.disease();

                if (!localMutationsByDisease.ContainsKey(disease))
                {
                    localMutationsByDisease.Add(disease, 0);
                }
                localMutationsByDisease[disease] += mafLines.Count();

                foreach (var mafLine in mafLines)
                {
                    if (!localMutationsByType.ContainsKey(mafLine.Variant_Classification))
                    {
                        localMutationsByType.Add(mafLine.Variant_Classification, 0);
                    }

                    localMutationsByType[mafLine.Variant_Classification]++;
                }
            }
        }
    }
}
