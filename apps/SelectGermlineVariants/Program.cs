using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ASELib;
using System.IO;
using System.Threading;
using System.Diagnostics;


namespace SelectGermlineVariants
{
    class Program
    {

        const long granularity = 1; // How often to select a variant if one is available. (This used to be 1000, but we're skipping the granularity here at the tentative stage, and doing it when we go from tentative to final, since we don't know what will be eliminated later and don't want to lose a candidate by choosing a bad one from a grain that otherwise has a good one.)
        const long isolationDistance = 150; // How many bases around a selected variant must match the germline exclusively.

        class CandidateVariant : IEquatable<CandidateVariant>, IComparable<CandidateVariant>
        {
            public CandidateVariant(string line_, long locus_, double odds_, string chromosome_, double alleleFrequency_, double alleleBlanace_, int alleleCount_, int alleleNumber_, string cigar_, bool goodcandidate_)
            {
                line = line_;
                locus = locus_;
                odds = odds_;
                chromosome = chromosome_;
                alleleFrequency = alleleFrequency_;
                alleleBalance = alleleBlanace_;
                alleleCount = alleleCount_;
                alleleNumber = alleleNumber_;
                cigar = cigar_;
                goodCandidate = goodcandidate_;
            }

            public bool Equals(CandidateVariant peer)
            {
                return peer.line == line;   // The rest is all derived from line, so equal lines are equivalent to equal objects.
            }

            public int CompareTo(CandidateVariant peer)
            {
                if (chromosome != peer.chromosome) return chromosome.CompareTo(peer.chromosome);
                return locus.CompareTo(peer.locus);
            }

            public readonly string chromosome;
            public readonly string line;
            public readonly long locus;
            public readonly double odds;
            public readonly double alleleFrequency;
            public readonly double alleleBalance;
            public readonly int alleleCount;
            public readonly int alleleNumber;
            public readonly string cigar;
            public readonly bool goodCandidate;
            public int nRNAReads = -1;
            public int nDNAReads = -1;
        }

        static void EmitBestCandidate(StreamWriter outputFile, string chromosome, List<CandidateVariant> liveCandidates, Dictionary<int, int> variantCountByOdds)
        {
            if (liveCandidates.Count() == 0)
            {
                return;
            }

            CandidateVariant bestCandiate = liveCandidates[0];

            for (int i = 0; i < liveCandidates.Count(); i++) // Yes, I know there's some clever c# way to do this with linq, but it would take me longer to look up than to write this loop
            {
                if (liveCandidates[i].odds > bestCandiate.odds)
                {
                    bestCandiate = liveCandidates[i];
                }
            }

            if (variantCountByOdds != null)
            {
                if (!variantCountByOdds.ContainsKey((int)bestCandiate.odds))
                {
                    variantCountByOdds.Add((int)bestCandiate.odds, 0);
                }

                variantCountByOdds[(int)bestCandiate.odds]++;
            }
            outputFile.WriteLine(chromosome + "\t" + bestCandiate.locus + "\t" + bestCandiate.line);
        }

        static void EliminateTooNearCandidates(List<CandidateVariant> liveCandidates, long locus)
        {
            var candidatesToEliminate = new List<CandidateVariant>();
            foreach (var candidate in liveCandidates)
            {
                if (candidate.locus + isolationDistance >= locus)
                {
                    candidatesToEliminate.Add(candidate);
                }
            }

            foreach (var candidate in candidatesToEliminate)
            {
                liveCandidates.Remove(candidate);
            }
        }


        static void markReadCount(Dictionary<string, Dictionary<long, CandidateVariant>> viableCandidates, bool dna, string contigName, long locus, int mappedReadCount)
        {
            contigName = contigName.ToLower();

            if (!viableCandidates.ContainsKey(contigName))
            {
                //
                // Probably a minor contig that doesn't have a high quality SNV called for it.  Just ignore it.
                //
                return;
            }

            if (viableCandidates[contigName].ContainsKey(locus))
            {
                if (dna && viableCandidates[contigName][locus].nDNAReads != -1 || !dna && viableCandidates[contigName][locus].nRNAReads != -1)
                {
                    Console.WriteLine("Got read count more than once for same variant " + contigName + ":" + locus);
                    throw new FormatException();
                }

                if (dna)
                {
                    viableCandidates[contigName][locus].nDNAReads = mappedReadCount;
                }
                else
                {
                    viableCandidates[contigName][locus].nRNAReads = mappedReadCount;
                }
            }
        }

        static void HandleOneCase(ASETools.Case case_, int state)
        { 
            StreamReader vcfFile = null;
            vcfFile = ASETools.CreateStreamReaderWithRetry(case_.vcf_filename);

            string line;
            while (null != (line = vcfFile.ReadLine()) && line.Count() != 0 && line[0] == '#') {
                // Skip the header lines
            }

            if (null == line || line.Count() == 0) {
                Console.WriteLine("Corrupt vcf: missing body: " + case_.vcf_filename);
                return;
            }

            Dictionary<int, int> variantCountByOdds = computeDistribution ? new Dictionary<int, int>() : null;
            var mafLines = ASETools.AVLTree<ASETools.MAFLine>.CreateFromList(ASETools.MAFLine.ReadFileOnlyOnePerLocus(case_.all_maf_lines_filename, case_.case_id, false));

            //
            // First, read in all of the variants, saving those that we can't immediately exclude because of one reason or another.
            //

            var liveCandidates = new List<CandidateVariant>();
            var previousGrainsCandidates = new List<CandidateVariant>();
            var savedGrains = new List<List<CandidateVariant>>();               // Grains and all the candidate variants in them


            //vcfuniq doesn't seen to actually work with the current version of the tools.  So, read in all the lines first, then sort and uniq them, then process them.
            var candidateVariants = new List<CandidateVariant>();

            while (null != (line = vcfFile.ReadLine()))
            {
                var fields = line.Split('\t');

                if (fields.Count() != 10)
                {
                    Console.WriteLine("Wrong number of fields (" + fields.Count() + " != 10) in vcf line: '" + line + "' in file " + case_.vcf_filename + ".  Ignoring file.");
                    return;
                }

                var infoFields = fields[7].Split(';');

                var info = new Dictionary<string, string>();
                foreach (var infoField in infoFields)
                {
                    var keyValue = infoField.Split('=');
                    if (keyValue.Count() != 2)
                    {
                        Console.WriteLine("Unable to parse info field '" + infoField + " in file " + case_.vcf_filename);
                        return;
                    }

                    info.Add(keyValue[0], keyValue[1]);
                }

                if (!info.ContainsKey("AN") || !info.ContainsKey("AC") || !info.ContainsKey("CIGAR") || !info.ContainsKey("DP") || !info.ContainsKey("AF") || !info.ContainsKey("AB") || !info.ContainsKey("ODDS"))
                {
                    Console.WriteLine("vcf line '" + line + " doesn't contain one or more required info fields.  Skipping file " + case_.vcf_filename);
                    return;
                }

                bool goodCandidate = true;
                double alleleFrequency = 0;
                double alleleBalance = 0;
                int alleleCount = 0;
                int alleleNumber = 0;
                double odds = 0;
                string cigar = info["CIGAR"];
                int locus = 0;
                double qual;

                try
                {
                    if (info["AF"].Contains(',')) // This happens with multiple alleles (typically, heterozygous with both different from the reference).  Don't try to parse them, just don't select it
                    {
                        goodCandidate = false;
                    }
                    else
                    {
                        alleleFrequency = Convert.ToDouble(info["AF"]);
                        alleleBalance = Convert.ToDouble(info["AB"]);
                        alleleCount = Convert.ToInt32(info["AC"]);
                        alleleNumber = Convert.ToInt32(info["AN"]);
                        odds = Convert.ToDouble(info["ODDS"]);
                    }
                    locus = Convert.ToInt32(fields[1]);
                    qual = Convert.ToDouble(fields[5]);
                }
                catch (FormatException)
                {
                    Console.WriteLine("Error parsing info fields in line " + line + " of file " + case_.vcf_filename + ".  Skipping file.");
                    return;
                }

                if (qual < 1)
                {
                    //
                    // Just ignore this.  Keeping it will result in the rejection of nearby decent variants.  These things typically are caused just by noise in the reads, low base-call quality bases, or
                    // other garbage.
                    //
                    continue;
                }

                string chromosome = fields[0];

                //
                // Eliminate any candidates that are too near to a somatic mutation (any somatic mutation, not just the ones we
                // select as interesting).
                //
                var key = new ASETools.MAFLine(chromosome, (int)locus);
                ASETools.MAFLine nextLowestMAFLine = null, nextHighestMAFLine = null;
                if (mafLines.FindFirstLessThanOrEqualTo(key, out nextLowestMAFLine) &&
                    nextLowestMAFLine.Chromosome == chromosome && nextLowestMAFLine.End_Positon + isolationDistance > locus)
                {
                    goodCandidate = false;
                }
                else if (mafLines.FindFirstGreaterThanOrEqualTo(key, out nextHighestMAFLine) &&
                    nextHighestMAFLine.Chromosome == chromosome && locus + isolationDistance > nextHighestMAFLine.Start_Position)
                {
                    goodCandidate = false;
                }
                /*
                    else if (repetitiveRegionMap.isCloseToRepetitiveRegion(ASETools.ChromosomeNameToIndex(fields[0]), (int)locus, (int)isolationDistance))
                {
                    goodCandidate = false;
                }
                */

                //
                // Eliminate any variant in a gene that also contains a somatic mutation that might cause nonsense mediated decay.
                // We're a little over cautious here, but since we're just eliminating potential germline variant ASE measurement
                // sites it's not a problem.
                //
                var genesMappedHere = geneMap.getGenesMappedTo(chromosome, locus).ToList();
                if (nextLowestMAFLine != null && genesMappedHere.Count() != 0)
                {
                    int minRange = genesMappedHere.Select(x => x.minLocus).Min();
                    for (var mafLine = nextLowestMAFLine;
                        mafLine != null &&
                        mafLine.Chromosome == chromosome && mafLine.Start_Position >= minRange;
                        mafLine = mafLines.FindFirstLessThan(mafLine))
                    {
                        if (mafLine.IsNonsenseMediatedDecayCausing())
                        {
                            goodCandidate = false;
                            break;
                        }
                    }
                }

                if (goodCandidate && nextHighestMAFLine != null && genesMappedHere.Count() != 0)
                {
                    int maxRange = genesMappedHere.Select(x => x.maxLocus).Max();
                    for (var mafLine = nextHighestMAFLine;
                        mafLine != null &&
                        mafLine.Chromosome == chromosome && mafLine.End_Positon <= maxRange;
                        mafLine = mafLines.FindFirstGreaterThan(mafLine))
                    {
                        if (mafLine.IsNonsenseMediatedDecayCausing())
                        {
                            goodCandidate = false;
                            break;
                        }
                    }
                }


                candidateVariants.Add(new CandidateVariant(line, locus, odds, chromosome.ToLower(), alleleFrequency, alleleBalance, alleleCount, alleleNumber, cigar, goodCandidate));

            } // while we have a VCF line.

            candidateVariants.Sort();
            var lociToDelete = Enumerable.Range(1, candidateVariants.Count() - 1).Where(_ => candidateVariants[_].chromosome == candidateVariants[_ - 1].chromosome && candidateVariants[_].locus == candidateVariants[_ - 1].locus).ToList();

            for (int i = lociToDelete.Count() - 1; i >= 0; i--) // Run through the list backwards so that when we delete an index, it doesn't affect the positions of the other ones we need to delete.
            {
                candidateVariants.RemoveAt(lociToDelete[i]);
            }

            string currentChromosome = "";
            long lastLocus = -isolationDistance - 1;

            foreach (var candidateVariant in  candidateVariants)
            { 
                //
                // Figure out if we've moved into another grain, in which case we save the candidates we have from previous grains
                //
                if (currentChromosome != candidateVariant.chromosome || candidateVariant.locus / granularity != lastLocus / granularity)
                {
                    savedGrains.Add(previousGrainsCandidates);
                    previousGrainsCandidates = liveCandidates;
                    liveCandidates = new List<CandidateVariant>();

                    if (currentChromosome != candidateVariant.chromosome)
                    {
                        //
                        // Starting a new chromosome, so we don't need to hang on to an old
                        // chromosome's candidates to make sure that we don't have any variants too close to
                        // the end of the old grain.
                        //
                        savedGrains.Add(previousGrainsCandidates);
                        previousGrainsCandidates = new List<CandidateVariant>();
                        currentChromosome = candidateVariant.chromosome;
                        lastLocus = -isolationDistance - 1;
                    }
                }

                EliminateTooNearCandidates(liveCandidates, candidateVariant.locus);
                EliminateTooNearCandidates(previousGrainsCandidates, candidateVariant.locus);

                var goodCandidate = candidateVariant.goodCandidate && candidateVariant.alleleFrequency == 0.5 && candidateVariant.alleleBalance > 0.4 && candidateVariant.alleleBalance < 0.6 && candidateVariant.alleleCount == 1 && 
                    candidateVariant.alleleNumber == 2 && candidateVariant.cigar == "1X" && (candidateVariant.odds > 20 || computeDistribution) && lastLocus + isolationDistance < candidateVariant.locus;

                if (goodCandidate)
                {
                    liveCandidates.Add(candidateVariant);
                }

                lastLocus = candidateVariant.locus;
                currentChromosome = candidateVariant.chromosome;
            } // candidate variants

            //
            // We now have a list of grains.  Make a map of the candidates in those grains that we can use to add in the DNA/RNA read count.
            //

            var viableCandidates = new Dictionary<string, Dictionary<long, CandidateVariant>>(); // Maps chromosome -> (locus -> candidate)

            foreach (var grain in savedGrains)
            {
                foreach (var candidateVariant in grain)
                {
                    if (!viableCandidates.ContainsKey(candidateVariant.chromosome))
                    {
                        viableCandidates.Add(candidateVariant.chromosome, new Dictionary<long, CandidateVariant>());
                    }

                    viableCandidates[candidateVariant.chromosome].Add(candidateVariant.locus, candidateVariant);
                }
            }

            //
            // Now read in the allcount files and use them to annotate the candidates.
            //

            ASETools.AllcountReader.ProcessBase processRNABase = (contigName, locus, mappedReadCount) => markReadCount(viableCandidates, false, contigName, locus, mappedReadCount);
            ASETools.AllcountReader.ProcessBase processDNABase = (contigName, locus, mappedReadCount) => markReadCount(viableCandidates, true, contigName, locus, mappedReadCount);
            var rnaAllcountReader = new ASETools.AllcountReader(case_.tumor_rna_allcount_filename);
            var dnaAllcountReader = new ASETools.AllcountReader(case_.tumor_dna_allcount_filename);
            long mappedHQNUclearReads;
            int numContigs;
            if (!dnaAllcountReader.openFile(out mappedHQNUclearReads, out numContigs))
            {
                Console.WriteLine("Couldn't open or bad header format in " + case_.tumor_dna_allcount_filename);
                return;
            }

            if (!dnaAllcountReader.ReadAllcountFile(processDNABase))
            {
                Console.WriteLine("Bad internal format or truncation in " + case_.tumor_dna_allcount_filename);
                return;
            }

            if (!rnaAllcountReader.openFile(out mappedHQNUclearReads, out numContigs)) {
                Console.WriteLine("Couldn't open or bad header format in " + case_.tumor_rna_allcount_filename);
                return;
            }
                    
            if (!rnaAllcountReader.ReadAllcountFile(processRNABase)) {
                Console.WriteLine("Bad internal format or truncation in " + case_.tumor_rna_allcount_filename);
                return;
            }

            //
            // Now run through the grains, select only the variants that have enough reads, and emit the best one for each grain.
            //
            var outputFilename = case_.vcf_filename.Substring(0, case_.vcf_filename.LastIndexOf('.')) + ASETools.tentativeSelectedVariantsExtension;
            var outputFile = computeDistribution ? StreamWriter.Null :  ASETools.CreateStreamWriterWithRetry(outputFilename);    // Don't actually write the output if we're computing distribution.
            outputFile.WriteLine("SelectGermlineVariants v1.1 for input file " + case_.vcf_filename);      // v1.0 didn't take into account the read counts when selecting variants.

            foreach (var grain in savedGrains)
            {
                var remainingCandidates = new List<CandidateVariant>();

                foreach (var candidate in grain)
                {
                    if (candidate.nDNAReads >= 10 && candidate.nRNAReads >= 10 && candidate.nDNAReads < 10000000 && candidate.nRNAReads < 1000000)  // Get rid of the extraorindarily high read counts because they can make extracted files > 2GB, which breaks the code
                    {
                        remainingCandidates.Add(candidate);
                    }
                }

                if (remainingCandidates.Count() > 0)
                {
                    EmitBestCandidate(outputFile, remainingCandidates[0].chromosome, remainingCandidates, variantCountByOdds);
                }
            }

            outputFile.WriteLine("**done**");
            outputFile.Close();


            if (computeDistribution)
            {
                lock (configuration)    // Overloading this lock here, but it's harmless
                {
                    var writer = ASETools.CreateAppendingStreamWriterWithRetry(configuration.finalResultsDirectory + "SelectedVariantDistribution.txt");
                    if (writer != null)
                    {
                        writer.Write(case_.case_id);

                        var maxOdds = variantCountByOdds.Select(x => x.Key).Max();
                        for (int i = 0; i < 20; i++)
                        {
                            writer.Write("\t" + (variantCountByOdds.ContainsKey(i) ? variantCountByOdds[i] : 0));
                        }
                        writer.WriteLine("\t" + variantCountByOdds.Where(x => x.Key >= 20).Select(x => x.Value).Sum());
                        writer.Close();
                    }
                }
            }
            vcfFile.Close();
        } // HandleOneCase


        static bool computeDistribution;
        static ASETools.Configuration configuration;
        //static ASETools.ASERepetitiveRegionMap repetitiveRegionMap;
        static ASETools.GeneLocationsByNameAndChromosome geneLocationInformation;
        static ASETools.GeneMap geneMap;

        static void Main(string[] args)
        {
            var timer = new Stopwatch();
            timer.Start();

            configuration = ASETools.Configuration.loadFromFile(args);

            if (null == configuration)
            {
                Console.WriteLine("Giving up because we were unable to load configuration.");
                return;
            }

            geneLocationInformation = new ASETools.GeneLocationsByNameAndChromosome(ASETools.readKnownGeneFile(configuration.geneLocationInformationFilename));
            geneMap = new ASETools.GeneMap(geneLocationInformation.genesByName);

            var caseIds = configuration.commandLineArgs.Where(x => x != "-d").ToList();

            if (caseIds.Count() == 0)
            {
                Console.WriteLine("usage: SelectGermlineVariants {-configuration configurationFile} {-d} <patrticipantIds>");
                Console.WriteLine("-d says to compute the distribution of how many variants we would select by varying the required odds from FreeBayes.");
                return;
            }

            computeDistribution = configuration.commandLineArgs.Contains("-d");

            var cases = ASETools.Case.LoadCases(configuration.casesFilePathname);

            int nValidCases = 0;

            var casesToProcess = new List<ASETools.Case>();

            foreach (var caseId in caseIds)
            {
                if (!cases.ContainsKey(caseId))
                {
                    Console.WriteLine(caseId + " does not appear to be a case ID; ignoring.");
                    continue;
                }

                if (cases[caseId].vcf_filename == "" || cases[caseId].tumor_rna_allcount_filename == "" || cases[caseId].all_maf_lines_filename == "")
                {
                    Console.WriteLine(caseId + " doesn't appear to have a complete set of vcf, allcount and all MAF lines files yet.  Ignoring.");
                    continue;
                }

                casesToProcess.Add(cases[caseId]);
                nValidCases++;
            }


            int nCasesPerDot;
            ASETools.PrintMessageAndNumberBar("Processing", "cases", casesToProcess.Count(), out nCasesPerDot);
            var threading = new ASETools.WorkerThreadHelper<ASETools.Case, int>(casesToProcess, HandleOneCase, null, null, nCasesPerDot);
            threading.run();


            timer.Stop();
            Console.WriteLine("Processed " + nValidCases + " cases in " + ASETools.ElapsedTimeInSeconds(timer));
        }
    }
}
