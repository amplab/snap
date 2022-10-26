using System;
using System.Collections.Generic;
using System.Linq;
using System.Net.Sockets;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using ASELib;
using Microsoft.Win32;

namespace LongReadCoverage
{
    internal class Program
    {
        static int nBitsSetInULong(ulong value)
        {
            if (value == ulong.MaxValue) return nBitsPerULong;
            if (value == 0) return 0;

            int nBitsSet = 0;
            for (int i = 0; i < nBitsPerULong; i++)
            {
                if ((((ulong)1 << i) & value) != 0)
                {
                    nBitsSet++;
                }
            }

            return nBitsSet;
        } // nBitsSetInULong

        class Contig
        {
            public Contig(string name_, long size_)
            {
                name = name_;
                size = size_;
                nCoverageBitsLongs = (int)((size + 63) / 64);

                coverageBits = new ulong[nCoverageBitsLongs];

                reset();
            }

            public void reset()
            {
                for (var i = 0; i < nCoverageBitsLongs; i++)
                {
                    coverageBits[i] = 0;
                }
            } // reset

            public void markCovered(int readLength, int offsetInContig)
            {
                //
                // Set the bits in three parts: up to the first whole long, all of the whole longs, and the remainder.
                //
                if (offsetInContig % nBitsPerULong != 0)
                {
                    int partialLongIndex = offsetInContig / nBitsPerULong;
                    for (int whichBit = offsetInContig % nBitsPerULong; whichBit < nBitsPerULong; whichBit++)
                    {
                        coverageBits[partialLongIndex] |= (uint)1 << whichBit;
                    }
                }

                for (int longIndex = (offsetInContig + nBitsPerULong - 1) / nBitsPerULong; longIndex < (offsetInContig + readLength) / nBitsPerULong; longIndex++)
                {
                    coverageBits[longIndex] = ulong.MaxValue;
                }

                int finalLongIndex = (offsetInContig + readLength) / nBitsPerULong;
                for (int whichBit = 0; whichBit < (offsetInContig + readLength) % nBitsPerULong; whichBit++)
                {
                    coverageBits[finalLongIndex] |= (uint)1 << whichBit;
                }
            } // markCovered

            public long totalBasesCovered() // This returns a long because while we will never have more than maxInt, sums of multiple contigs will
            {
                return coverageBits.Select(_ => nBitsSetInULong(_)).Sum();
            }

            public bool isBaseCovered(long offset)
            {
                return (coverageBits[offset / nBitsPerULong] & ((ulong)1 << (int)(offset % nBitsPerULong))) != 0;
            }

            public List<long> chunkSizes()
            {
                List<long> sizes = new List<long>();

                long currentChunkStart = -1;
                bool inChunk = false;

                for (long i = 0; i < size; i++)
                {

                }
            }

            public readonly string name;
            public readonly long size;
            public long startingOffset = 0;
            int nCoverageBitsLongs;
            ulong[] coverageBits;
        } // Contig


        static List<Contig> contigs = new List<Contig>();
        static long totalSize;
        static Random random = new Random();

        const int nBitsPerULong = 64;

        static long getRandom(long maxValue)    // This is a bad implementation and doesn't give really uniform coverage.  Sue me.
        {
            if (maxValue <= int.MaxValue)
            {
                return random.Next((int)maxValue);
            }

            return (random.Next(int.MaxValue) + ((long)1 << 32) * random.Next(int.MaxValue)) % maxValue;
        } // getRandom

        static void randomlyChooseContigAndOffset(int readSize, out int contigNumber, out int offsetInContig)
        {
            while (true)
            {
                long rawOffset = getRandom(totalSize);

                if (rawOffset + readSize >= totalSize)
                {
                    continue;
                }

                long totalOffsetForContigs = 0;
                for (contigNumber = 0; contigNumber < contigs.Count(); contigNumber++)
                {
                    if (totalOffsetForContigs + contigs[contigNumber].size > rawOffset)
                    {
                        // It starts in this contig
                        if (totalOffsetForContigs + contigs[contigNumber].size <= rawOffset + readSize)
                        {
                            //
                            // But it doens't end in the contig.
                            //
                            break;
                        }
                        offsetInContig = (int)(rawOffset - totalOffsetForContigs);  // Since contigs can't be bigger than int.MaxValue
                        return;
                    }
                    totalOffsetForContigs += contigs[contigNumber].size;
                } // for each contig
            } // forever
        } // randomlyChooseContigAndOffset

        static double ComputeRawFractionCovered(int readLength, double coverage)
        {
            contigs.ForEach(_ => _.reset());

            long nReads = (long)((double)totalSize * coverage / readLength);

            for (long i = 0; i < nReads; i++)
            {
                int contigNumber, offsetInContig;

                randomlyChooseContigAndOffset(readLength, out contigNumber, out offsetInContig);
                contigs[contigNumber].markCovered(readLength, offsetInContig);
            }

            var totalBasesCovered = contigs.Select(_ => _.totalBasesCovered()).Sum();

            return (double)totalBasesCovered / totalSize;

        } // ComputeRawFractionCovered


        static void Main(string[] args)
        {
            //
            // Sizes from hg38
            //

            contigs.Add(new Contig("chr1", 248956422));
            contigs.Add(new Contig("chr2", 242193529));
            contigs.Add(new Contig("chr3", 198295559));
            contigs.Add(new Contig("chr4", 190214555));
            contigs.Add(new Contig("chr5", 181538259));
            contigs.Add(new Contig("chr6", 170805979));
            contigs.Add(new Contig("chr7", 159345973));
            contigs.Add(new Contig("chr8", 145138636));
            contigs.Add(new Contig("chr9", 138394717));
            contigs.Add(new Contig("chr10", 133797422));
            contigs.Add(new Contig("chr11", 135086622));
            contigs.Add(new Contig("chr12", 133275309));
            contigs.Add(new Contig("chr13", 114364328));
            contigs.Add(new Contig("chr14", 107043718));
            contigs.Add(new Contig("chr15", 101991189));
            contigs.Add(new Contig("chr16", 90338345));
            contigs.Add(new Contig("chr17", 83257441));
            contigs.Add(new Contig("chr18", 80373285));
            contigs.Add(new Contig("chr19", 58617616));
            contigs.Add(new Contig("chr20", 64444167));
            contigs.Add(new Contig("chr21", 46709983));
            contigs.Add(new Contig("chr22", 50818468));
            contigs.Add(new Contig("chrX", 156040895));
            contigs.Add(new Contig("chrY", 57227415));

            long startingOffset = 0;
            for (int i = 0; i < contigs.Count(); i++)
            {
                contigs[i].startingOffset = startingOffset;
                startingOffset += contigs[i].size;
            }

            totalSize = contigs.Select(_ => _.size).Sum();

            int [] readSizes = { 10000, 30000, 100000, 1000000 };

            foreach (var readSize in readSizes)
            {
                for (double coverage = 1; coverage < 10; coverage += 0.5)
                {
                    Console.WriteLine("Read size " + readSize + " total read size " + coverage + "x genome yields " + ComputeRawFractionCovered(readSize, coverage) + " coverage");
                }

                Console.WriteLine();

                for (double coverage = 1; coverage < 10; coverage += 0.5)
                {
                    Console.WriteLine("Read size " + readSize + " total read size " + coverage + "x genome yields " + ComputeRawFractionCovered(readSize - 1000, coverage * ((double)readSize - 1000) / readSize) + " coverage requiring 1kb of overlap");
                }

                Console.WriteLine();
            }

            contigs.Add(new Contig("chr1-2", 248956422));
            contigs.Add(new Contig("chr2-2", 242193529));
            contigs.Add(new Contig("chr3-2", 198295559));
            contigs.Add(new Contig("chr4-2", 190214555));
            contigs.Add(new Contig("chr5-2", 181538259));
            contigs.Add(new Contig("chr6-2", 170805979));
            contigs.Add(new Contig("chr7-2", 159345973));
            contigs.Add(new Contig("chr8-2", 145138636));
            contigs.Add(new Contig("chr9-2", 138394717));
            contigs.Add(new Contig("chr10-2", 133797422));
            contigs.Add(new Contig("chr11-2", 135086622));
            contigs.Add(new Contig("chr12-2", 133275309));
            contigs.Add(new Contig("chr13-2", 114364328));
            contigs.Add(new Contig("chr14-2", 107043718));
            contigs.Add(new Contig("chr15-2", 101991189));
            contigs.Add(new Contig("chr16-2", 90338345));
            contigs.Add(new Contig("chr17-2", 83257441));
            contigs.Add(new Contig("chr18-2", 80373285));
            contigs.Add(new Contig("chr19-2", 58617616));
            contigs.Add(new Contig("chr20-2", 64444167));
            contigs.Add(new Contig("chr21-2", 46709983));
            contigs.Add(new Contig("chr22-2", 50818468));
            contigs.Add(new Contig("chrX-2", 156040895));
            contigs.Add(new Contig("chrY-2", 57227415));

            startingOffset = 0;
            for (int i = 0; i < contigs.Count(); i++)
            {
                contigs[i].startingOffset = startingOffset;
                startingOffset += contigs[i].size;
            }

            totalSize = contigs.Select(_ => _.size).Sum();

            foreach (var readSize in readSizes)
            {
                for (double coverage = 1; coverage < 10; coverage += 0.5)
                {
                    Console.WriteLine("Read size " + readSize + " total read size " + coverage + "x (haploid) genome yields " + ComputeRawFractionCovered(readSize, coverage/2) + " coverage");
                }

                Console.WriteLine();

                for (double coverage = 1; coverage < 10; coverage += 0.5)
                {
                    Console.WriteLine("Read size " + readSize + " total read size " + coverage + "x (haploid) genome yields " + 
                        ComputeRawFractionCovered(readSize - 4000, coverage / 2 * ((double)readSize - 4000) / readSize) + " coverage requiring 4kb of overlap (this is what you'd need for phasing)");
                }

                Console.WriteLine();
            }

        } // Main
    } // Program
}
