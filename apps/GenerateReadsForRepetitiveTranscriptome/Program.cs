using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;
using System.IO;
using ASELib;

namespace GenerateReadsForRepetitiveTranscriptome
{
    class Program
    {

        static void Main(string[] args)
        {
            Stopwatch timer = new Stopwatch();
            timer.Start();

            var configuration = ASETools.Configuration.loadFromFile(args);

            if (null == configuration)
            {
                return;
            }

            if (configuration.commandLineArgs.Count() != 3 && configuration.commandLineArgs.Count() != 4)
            {
                Console.WriteLine("usage: GenerateReadsForRepetitiveTranscriptome snapIndex generatedReadsFilename chromosomeNumber {generatedFastaFilename}");
                return;
            }

            string snapIndexName = configuration.commandLineArgs[0];
            string generatedReadsFilename = configuration.commandLineArgs[1];
            int chromosomeNumberToProcess;
            try
            {
                chromosomeNumberToProcess = Convert.ToInt32(configuration.commandLineArgs[2]);
            } catch
            {
                Console.WriteLine("Chromosome number must be a number,");
                return;
            }

            string generatedFastaName;
            if (configuration.commandLineArgs.Count() == 4)
            {
                generatedFastaName = configuration.commandLineArgs[3];
            } else
            {
                generatedFastaName = null;
            }


            Console.Write("Loading gene information...");

            var geneLocationInformation = new ASETools.GeneLocationsByNameAndChromosome(ASETools.readKnownGeneFile(configuration.geneLocationInformationFilename));
            var geneMap = new ASETools.GeneMap(geneLocationInformation.genesByName);

            var transcriptomeContigs = new Dictionary<string, Dictionary<ASETools.Exon, List<List<ASETools.Exon>>>>();  // chromosome->starting exon->{list of exons in contig}
            var includedExons = new Dictionary<string, Dictionary<ASETools.Exon, bool>>();  // chromosome->exons that are included anywhere in a contig.  Used to eliminate singleton exons that exist in a contig elsewhere.


            //
            // We're looking for places where exons abut different exons within the read length limit.
            //

            int nContigs = 0;
            long nBases = 0;
            int nSingletons = 0;
            int nMoreThanTwo = 0;

            Console.WriteLine(ASETools.ElapsedTimeInSeconds(timer));
            Console.Write("Processing transcriptome...");

            foreach (var gene in geneLocationInformation.genesByName.Select(_ => _.Value).ToList())
            {
                if (gene.inconsistent)
                {
                    continue;
                }
#if false
Console.WriteLine(gene.hugoSymbol + " has " + gene.isoforms.Count() + " isoforms.");

var exonTable = new Dictionary<ASETools.Exon, char>();
char nextExon = 'A';
                
for (int i = 0; i < gene.isoforms.Count(); i++)
{
    var isoform = gene.isoforms[i];
    Console.WriteLine("isoform " + i + " has " + isoform.exons.Count() + " exons :");
    for (int j = 0; j < isoform.exons.Count(); j++)
    {
        var exon = isoform.exons[j];
        if (!exonTable.ContainsKey(exon))
        {
            exonTable.Add(exon, nextExon);
            nextExon++;
        }
        Console.WriteLine("\tExon " + j + "(" + exonTable[exon] + "): start: " + exon.start + ", end: " + exon.end + ", length: " + exon.length());
    }
}
#endif

                var currentChromosome = gene.chromosome;
                if (!transcriptomeContigs.ContainsKey(currentChromosome))
                {
                    transcriptomeContigs.Add(currentChromosome, new Dictionary<ASETools.Exon, List<List<ASETools.Exon>>>());
                    includedExons.Add(currentChromosome, new Dictionary<ASETools.Exon, bool>());
                }

                if (gene.isoforms.Any(_ => _.exons.Any(e => e.start <= 1 || e.end <= 1))) {
                    Console.WriteLine("Bogus exon for gene " + gene.hugoSymbol);
                    continue;
                }

                var contigsForThisChromosome = transcriptomeContigs[currentChromosome];

                foreach (var isoform in gene.isoforms)
                {
                    for (int i = 0; i < isoform.exons.Count(); i++)
                    {
                        var exon = isoform.exons[i];

                        var newContig = new List<ASETools.Exon>();
                        newContig.Add(exon);

                        var lengthAdded = 0;
                        for (int j = i + 1; j < isoform.exons.Count() && lengthAdded < configuration.readLengthForRepetitiveRegionDetection + 5 /* for indels */; j++)
                        {
                            newContig.Add(isoform.exons[j]);
                            lengthAdded += isoform.exons[j].length();
                        }

                        if (newContig.Count() == 1 && includedExons[gene.chromosome].ContainsKey(exon))
                        {
                            //
                            // This exon is already there in an existing contig, no need to add it as a singleton.
                            //
                            continue;
                        }

                        if (newContig.Select(_ => _.length()).Sum() < configuration.readLengthForRepetitiveRegionDetection)
                        {
                            //
                            // Too short to map.  Skip it.
                            //
                            continue;
                        }

                        if (!contigsForThisChromosome.ContainsKey(exon))
                        {
                            contigsForThisChromosome.Add(exon, new List<List<ASETools.Exon>>());
                            contigsForThisChromosome[exon].Add(newContig);
                            nContigs++;
                            nBases += newContig.Select(_ => _.end - _.start).Sum();

                            if (newContig.Count() == 1)
                            {
                                nSingletons++;
                            } else if (newContig.Count() > 2)
                            {
                                nMoreThanTwo++;
                            }

//Console.WriteLine("Adding " + new string(newContig.Select(_ => exonTable[_]).ToArray()));

                            addExonsToList(gene.chromosome, newContig, includedExons);
                        } else
                        {
                            //
                            // Figure out if this this is already there.
                            //
                            bool found = false;
                            foreach (var candidate in contigsForThisChromosome[exon])
                            {
                                if (candidate.Count() != newContig.Count())
                                {
                                    continue;
                                }

                                bool allMatch = true;
                                for (int j = 0; j < candidate.Count(); j++)
                                {
                                    if (!candidate[j].Equals(newContig[j]))
                                    {
                                        allMatch = false;
                                        break;
                                    }
                                }

                                if (allMatch)
                                {
                                    found = true;
                                    break;
                                }
                            } // foreach candidate

                            if (!found)
                            {
                                contigsForThisChromosome[exon].Add(newContig);
                                nContigs++;
                                nBases += newContig.Select(_ => _.end - _.start).Sum();
                                addExonsToList(gene.chromosome, newContig, includedExons);

//Console.WriteLine("Adding " + new string(newContig.Select(_ => exonTable[_]).ToArray()));

                                if (newContig.Count() == 1)
                                {
                                    nSingletons++;
                                }
                                else if (newContig.Count() > 2)
                                {
                                    nMoreThanTwo++;
                                }
                            }
                        }
                    } // exons
                } // isoforms
            } // genes

            Console.WriteLine(ASETools.ElapsedTimeInSeconds(timer) + " (total elapsed)");
            Console.WriteLine("Total of " + nContigs + " contigs and " + nBases + " bases.  " + nSingletons + " singletons and " + nMoreThanTwo + " contigs with > 2 exons");

            Console.Write("Loading genome...");
            var genome = new ASETools.Genome();

            if (!genome.load(snapIndexName))
            {
                Console.WriteLine();
                Console.WriteLine("Unable to load genome from directory " + args[0]);
                return;
            }
            Console.WriteLine(ASETools.ElapsedTimeInSeconds(timer) + " (total elapsed)");
            Console.Write("Generating reads" + ((generatedFastaName == null) ? "" : " and fasta") + "...");

            var outputFile = ASETools.CreateStreamWriterWithRetry(generatedReadsFilename);
            if (null == outputFile)
            {
                Console.WriteLine("Unable to open output file " + generatedReadsFilename);
                return;
            }


            outputFile.WriteLine("@HD\tVN:1.4");
            outputFile.Write("@PG\tID:GRRT.1\tPG:GenerateReadsForRepetitiveTranscriptome\tCL:");

            StreamWriter outputFasta;

            if (generatedFastaName == null)
            {
                outputFasta = StreamWriter.Null;
            } else {
                outputFasta = ASETools.CreateStreamWriterWithRetry(generatedFastaName);
                if (null == outputFasta)
                {
                    Console.WriteLine("Unable to open " + outputFasta);
                    return;
                }
            }

//outputFile = StreamWriter.Null;
//outputFasta = StreamWriter.Null;

            string qual = "";
            for (int i = 0; i < configuration.readLengthForRepetitiveRegionDetection; i++)
            {
                qual += "2";
            }

            foreach (var chromosomeEntry in transcriptomeContigs)
            {
                var chromosome = chromosomeEntry.Key;
                var mapOfLists = chromosomeEntry.Value;

                if (generatedFastaName == null && chromosome != "chr" + chromosomeNumberToProcess)
                {
                    continue;
                }

                foreach (var listOfLists in mapOfLists.Select(_ => _.Value)) {
                    foreach (var exonList in listOfLists)
                    {
                        int totalLength = exonList.Select(_ => _.length()).Sum();

                        // Write the new transcript contig into the fasta we're generating.
                        string contigName = "transcript_"+ chromosome;
                        foreach (var exon in exonList)
                        {
                            contigName += "_" + exon.start + "_" + exon.end;
                        }

                        outputFasta.WriteLine(">" + contigName);

                        const int dataLineMaxLength = 80;
                        string dataLine = "";

                        foreach (var exon in exonList)
                        {
                            for (int offsetInExon = 0; offsetInExon < exon.length(); offsetInExon++)
                            {
                                dataLine += char.ToUpper(genome.getBase(chromosome, exon.start + offsetInExon));
                                if (dataLine.Length >= dataLineMaxLength)
                                {
                                    outputFasta.WriteLine(dataLine);
                                    dataLine = "";
                                }
                            } // foreach base in the exon
                        } // foreach exon in the contig

                        if (dataLine != "")
                        {
                            outputFasta.WriteLine(dataLine);
                        }

                        if (chromosome != "chr" + chromosomeNumberToProcess)
                        {
                            continue;
                        }

                        // generate reads for the transcript contig
                        for (int offsetInContig = 0; offsetInContig < totalLength - configuration.readLengthForRepetitiveRegionDetection; offsetInContig++)
                        {
                            string readId = "trans." + chromosome + ".";
                            int offsetAccountedFor = 0;
                            int exonNumber = 0;
                            int offsetInExon = 0;
                            while (offsetAccountedFor < offsetInContig)
                            {
                                if (offsetAccountedFor + exonList[exonNumber].length() < offsetInContig - offsetAccountedFor)
                                {
                                    offsetAccountedFor += exonList[exonNumber].length();
                                    exonNumber++;
                                } else
                                {
                                    offsetInExon = offsetInContig - offsetAccountedFor;
                                    offsetAccountedFor = offsetInContig;
                                }
                            } // While we're still consuming exons prior to the start

                            readId += (exonList[exonNumber].start + offsetInExon);

                            string readData = "";
                            int readLength = 0;

                            while (readLength < configuration.readLengthForRepetitiveRegionDetection)
                            {
                                if (offsetInExon >= exonList[exonNumber].length())
                                {
                                    readId += "." + exonList[exonNumber].end + "." + exonList[exonNumber + 1].start;
                                    exonNumber++;
                                    offsetInExon = 0;
                                }

                                readData += char.ToUpper(genome.getBase(chromosome, exonList[exonNumber].start + offsetInExon));
                                readLength++;
                                offsetInExon++;
                            }

                            readId += "." + (exonList[exonNumber].start + offsetInExon);
                            outputFile.WriteLine(readId + "\t2\t" + contigName + "\t" + (offsetInContig + 1) + "\t70\t" + readLength + "M\t*\t*\t" + readLength + "\t" + readData + "\t" + qual);
                        } // for each offset in the contig that's small enough to generate a read
                    }
                }
            }

            outputFile.Close();
            outputFasta.Close();
            Console.WriteLine("Total elapsed time " + ASETools.ElapsedTimeInSeconds(timer));
        } // Main

        static void addExonsToList(string chromosome, List<ASETools.Exon> newContig, Dictionary<string, Dictionary<ASETools.Exon, bool>> includedExons)
        {
            foreach (var exon in newContig)
            {
                if (!includedExons[chromosome].ContainsKey(exon))
                {
                    includedExons[chromosome].Add(exon, true);
                }
            }
        }

    }
}
