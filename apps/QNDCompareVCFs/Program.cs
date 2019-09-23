using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Diagnostics;
using ASELib;

namespace QNDCompareVCFs
{
    class Program
    {

        class Location : IComparable<Location>
        {
            public readonly string chromosome;
            public readonly int pos;

            public int CompareTo(Location peer)
            {
                if (chromosome == peer.chromosome)
                {
                    return pos.CompareTo(peer.pos);
                }
                return chromosome.CompareTo(peer.chromosome);
            }

            public Location(string chromosome_, int locus_)
            {
                chromosome = chromosome_;
                pos = locus_;
            }

            static public bool operator<(Location first, Location second)
            {
                return first.CompareTo(second) < 0;
            }

            static public bool operator>(Location first, Location second)
            {
                return first.CompareTo(second) > 0;
            }

            static public bool operator==(Location first, Location second)
            {
                return first.CompareTo(second) == 0;
            }

            static public bool operator!=(Location first, Location second)
            {
                return first.CompareTo(second) != 0;
            }

            static public bool operator<=(Location first, Location second)
            {
                return first.CompareTo(second) <= 0;
            }

            static public bool operator>=(Location first, Location second)
            {
                return first.CompareTo(second) >= 0;
            }

            override public string ToString()
            {
                return chromosome + "\t" + pos;
            }

            public override bool Equals(object peer)
            {
                return this == (Location)peer;
            }

            public override int GetHashCode()
            {
                return ((object)(this)).GetHashCode();  // We use all the fields in ==, so it's OK to do this
            }
        }

        const int CHROM = 0;
        const int POS = 1;
        const int REF = 3;
        const int ALT = 4;
        const int QUAL = 5;

        static void Main(string[] args)
        {
            if (args.Count() != 2 && args.Count() != 3)
            {
                Console.WriteLine("Usage: QNDCompareVCFs vcf1 vcf2 {truthVCF}");
                return;
            }

            int nStreams = args.Count();

            var inputStreams = new StreamReader[nStreams];
            var nextLines = new string[nStreams];
            var lastLocation = new Location[nStreams];

            var lastQual = new double[nStreams];
            int nDifferentLengthIndels = 0;
            int nRefMismatch = 0;
            int nAltMismatch = 0;
            var nCalledInOnly = new int[nStreams];
            var nLines = 0;
            int nMissingTrueVariant = 0;

            for (int i = 0; i < nStreams; i++)
            {
                if (args[i].ToLower().EndsWith(".gz"))
                {
                    inputStreams[i] = ASETools.CreateCompressedStreamReaderWithRetry(args[i]);
                } else
                {
                    inputStreams[i] = ASETools.CreateStreamReaderWithRetry(args[i]);
                }

                if (inputStreams[i] == null)
                {
                    Console.WriteLine("Unable to open input file " + args[i]);
                    return;
                }

                //
                // Skip over the headers, which start with '#'
                //
                while ((nextLines[i] = inputStreams[i].ReadLine()) != null && nextLines[i].StartsWith("#"))
                {
                    // This loop body intentionally left blank
                }

                var fields = nextLines[i].Split('\t');
                lastLocation[i] = new Location(fields[CHROM], Convert.ToInt32(fields[POS]));
                lastQual[i] = Convert.ToDouble(fields[QUAL]);

                nCalledInOnly[i] = 0;
            } // for each input stream

            int min_qual = 30;  // Ignore calls with lower qual than this.

            Console.Write("Comparing " + args[0] + " and " + args[1]);
            if (nStreams == 3)
            {
                Console.Write(" with truth set " + args[2]);
            }

            Console.WriteLine(" and min qual " + min_qual);
            Console.WriteLine("Chromosome\tLocus\tmax qual\tDifference");

            while (nextLines.Any(_ => _ != null))
            {
                nLines++;
                int fileToRead = -1;

                var fields = new string[nStreams][];
                for (int i = 0; i < nStreams; i++)
                {
                    fields[i] = nextLines[i].Split('\t');
                    if (fields[i].Count() != 10)
                    {
                        Console.WriteLine(args[i] + " has an unparsable line: " + nextLines[i]);
                        return;
                    }

                    lastLocation[i] = new Location(fields[i][0], Convert.ToInt32(fields[i][1]));
                    lastQual[i] = Convert.ToDouble(fields[i][5]);

                }

                if (lastLocation.Any(_ => _ != lastLocation[0]))
                {
                    //
                    // A locus that doesn't occur in all inputs.
                    //

                    if (nStreams > 2 && (lastLocation[2] < lastLocation[0] || lastLocation[2] < lastLocation[1]))
                    {
                        //
                        // Called in the truth set but missing from both.
                        //
                        Console.Write(lastLocation[2] + "\tAppears in the truth set but not in ");
                        if (lastLocation[2] < lastLocation[0])
                        {
                            if (lastLocation[2] < lastLocation[1])
                            {
                                Console.WriteLine("either variant call.");
                            }
                            else
                            {
                                Console.WriteLine("0");
                            }
                        }
                        else
                        {
                            Console.WriteLine("1");
                        }
                        nMissingTrueVariant++;
                        fileToRead = 2;
                    } else
                    {
                        //
                        // Called in one loctaion but not the truth set.
                        //
                        int lower = lastLocation[0] < lastLocation[1] ? 0 : 1;
                        Console.WriteLine(lastLocation[lower] + "\tCalled only in " + lower);
                        fileToRead = lower;
                    }
                }
                else
                {
                    //
                    // Same location in all input files.
                    //

                    fileToRead = nStreams; // Means read all

                    if (lastQual.Max() >= min_qual)
                    {
                        if (fields.Any(_ => _[REF].Length != fields[0][REF].Length))
                        {
                            //
                            // Sometimes, freebayes will call indels of different lengths just because it keeps more of the reference in one or the other.
                            // Take the longest one and pad out the rest with it.
                            //
                            int longest = 0;
                            for (int i = 1; i < nStreams; i++)
                            {
                                if (fields[i][REF].Length > fields[longest][REF].Length)
                                {
                                    longest = i;
                                }
                            }

                            for (int i = 0; i < nStreams; i++)
                            {
                                if (fields[i][REF].Length < fields[longest][REF].Length)
                                {
                                    fields[i][ALT] = fields[i][ALT] + fields[longest][REF].Substring(fields[i][REF].Length);
                                    fields[i][REF] = fields[longest][REF];
                                }
                            }
                        } // If they didn't all have the same length REF

                        if (fields.Any(_ => _[REF] != fields[0][REF]))
                        {
                            Console.Write(lastLocation[0] + "\tREF doesn't match");
                            fields.ToList().ForEach(_ => Console.Write(" " + _[REF]));
                            Console.WriteLine();
                            nRefMismatch++;

                        }
                        else if (fields.Any(_ => _[ALT] != fields[0][ALT]))
                        {
                            Console.Write(lastLocation[0] + "\tALT mismatch");
                            fields.ToList().ForEach(_ => Console.Write(" " + _[ALT]));
                            Console.WriteLine();
                            nAltMismatch++;
                        }
                    } // Qual high enough

                } // saw the same locus

                for (int i = 0; i < nStreams; i++)
                {
                    if (i == fileToRead || fileToRead == nStreams)
                    {
                        nextLines[i] = inputStreams[i].ReadLine();
                    }
                } // for each file to read
            } // while we have data

            inputStreams.ToList().ForEach(_ => _.Close());

            Console.WriteLine("Saw " + nLines + " different called loci.  Of those, " + nDifferentLengthIndels + " had different length indels, " + nRefMismatch + " had ref mismatches, " + nAltMismatch + " had alt mismatches " + ", " + nCalledInOnly[0] + " were only in the first input and " + nCalledInOnly[1] + " were only in the second.  The remainder matched.");


        } // Main
    }
}
