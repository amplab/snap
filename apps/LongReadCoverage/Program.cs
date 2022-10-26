using System;
using System.Collections.Generic;
using System.Linq;
using System.Net.Sockets;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using ASELib;

namespace LongReadCoverage
{
    internal class Program
    {
        class Contig
        {
            public Contig(string name_, int size_)
            {
                name = name_;
                size = size_;
                nCoverageBitsLongs = (size + 63) / 64;

                coverageBits = new long[nCoverageBitsLongs];

                reset();
            }

            public void reset()
            {
                for (var i = 0; i < nCoverageBitsLongs; i++)
                {
                    coverageBits[i] = 0;
                }
            } // reset

            public readonly string name;
            public readonly int size;
            public long startingOffset = 0;
            int nCoverageBitsLongs;
            long[] coverageBits;
        } // Contig


        static List<Contig> contigs = new List<Contig>();
        static long totalSize;
        static Random random = new Random();

        static long getRandomLong()
        {
            long result = random.Next();
            result |= ((long)random.Next()) << 32;

            if (result < 0) result = -result;
            return result;
        } // getRandomLong

        static void randomlyChooseContigAndOffset(int readSize, out int contigNumber, out int offsetInContig)
        {
            while (true)
            {
                long rawOffset = random.NextInt64();
            }
        }

        static double ComputeRawFractionCovered(long readLength, double coverage)
        {
            contigs.ForEach(_ => _.reset());



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
            contigs.Add(new Contig("chrM", 16569));

            long startingOffset = 0;
            for (int i = 0; i < contigs.Count(); i++)
            {
                contigs[i].startingOffset = startingOffset;
                startingOffset += contigs[i].size;
            }

            totalSize = contigs.Select(_ => _.size).Sum();


        } // Main
    } // Program
}
