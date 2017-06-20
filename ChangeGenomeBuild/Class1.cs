using System;
using System.Collections.Generic;
using System.Linq;
using ASELib;
using System.Text;
using System.Threading.Tasks;
using System.IO;

namespace GenomeBuild
{

	public class Interval
	{
		public readonly string name;
		public readonly long start;
		public readonly long end;
		public readonly char strand;

		public readonly string targetName;
		public readonly long targetStart;
		public readonly long targetEnd;

		public Interval(string name_, long start_, long end_, char strand_ = '+', string targetName_ = "", long targetStart_ = 0, long targetEnd_ = 0)
		{
			if (start_ > end_)
			{
				throw new Exception("Invalid Interval");
			}

			name = name_;
			start = start_;
			end = end_;

			targetName = targetName_;
			targetStart = targetStart_;
			targetEnd = targetEnd_;
			strand = strand_;
		}

		// Gets intersection of intervals
		public Interval Intersection(Interval other)
		{
			if (this.name != other.name)
			{
				return null;
			} else if(this.start > other.end || this.end < other.start) // the intervals do not overlap
			{
				return null;
			}

			return new Interval(this.name, Math.Max(this.start, other.start), Math.Min(this.end, other.end), strand, targetName, targetStart, targetEnd);

		}
	}
    public class LiftOver
    {


		// dictionary of genome build to genome build mappins
		Dictionary<string, List<Interval>> genomeMap = new Dictionary<string, List<Interval>>();

		// size of source chromosomes
		Dictionary<string, long> src_chrSize = new Dictionary<string, long>();

		// size of target chromosomes
		Dictionary<string, long> target_chrSize = new Dictionary<string, long>();

		public void readChainFile(string filename)
		{

			if (!filename.EndsWith(".chain"))
			{
				throw new Exception("invalid chain file format");
			}

			StreamReader inputFile = ASETools.CreateStreamReaderWithRetry(filename);

			// information that changes for each new chromosome
			long srcFrom = -1;
			long targetFrom = -1;

			long score = -1;
			string srcName = "";
			long srcSize = -1;
			char srcStrand = '.';
			long srcStart = -1;
			long srcEnd = -1;

			string targetName = "";
			long targetSize = -1;
			char targetStrand = '.';
			long targetStart = -1;
			long targetEnd = -1;

			long size = -1;

			string line;
			while ((line = inputFile.ReadLine()) != null)
			{
				string[] split = line.Split(null);

				if (split[0] == "chain")
				{
					score = Convert.ToInt64(split[1]);
					srcName = split[2];
					srcSize = Convert.ToInt64(split[3]);
					srcStrand = split[4].ToCharArray()[0];
					srcStart = Convert.ToInt64(split[5]);
					srcEnd = Convert.ToInt64(split[6]);

					if (src_chrSize.ContainsKey(srcName))
					{
						src_chrSize[srcName] = srcSize;
					}
					else {
						src_chrSize.Add(srcName, srcSize);
					}

					targetName = split[7];
					targetSize = Convert.ToInt64(split[8]);
					targetStrand = split[9].ToCharArray()[0];
					targetStart = Convert.ToInt64(split[10]);
					targetEnd = Convert.ToInt64(split[11]);

					if (target_chrSize.ContainsKey(targetName))
					{
						target_chrSize[targetName] = targetSize;
					}
					else
					{
						target_chrSize.Add(targetName, targetSize);
					}


					srcFrom = srcStart;
					targetFrom = targetStart;

					if (!genomeMap.ContainsKey(srcName))
						genomeMap.Add(srcName, new List<Interval>());

				} else if (split.Length == 3 && split[0] != "chain")
				{
					size = Convert.ToInt64(split[0]);
					long sgap = Convert.ToInt64(split[1]);
					long tgap = Convert.ToInt64(split[2]);

					if (targetStrand == '+')
					{
						genomeMap[srcName].Add(new Interval(srcName, srcFrom, srcFrom + size, targetStrand, targetName, targetFrom, targetFrom + size));


					} else if (targetStrand == '-')
					{
						genomeMap[srcName].Add(new Interval(srcName, srcFrom, srcFrom + size, targetStrand, targetName, targetSize - (targetFrom + size), targetSize - targetFrom));
					}
					else
					{
						throw new Exception("invalid target strand " + targetStrand);
					}

					srcFrom += size + sgap;
					targetFrom += size + tgap;


				} else if (split.Length == 1 && line.Length > 0 && split[0] != "chain")
				{
					size = Convert.ToInt64(split[0]);

					if (targetStrand == '+')
					{
						genomeMap[srcName].Add(new Interval(srcName, srcFrom, srcFrom + size, targetStrand, targetName, targetFrom, targetFrom + size));


					}
					else if (targetStrand == '-')
					{
						genomeMap[srcName].Add(new Interval(srcName, srcFrom, srcFrom + size, targetStrand, targetName, targetSize - (targetFrom + size), targetSize - targetFrom));
					}
					else
					{
						throw new Exception("invalid target strand " + targetStrand);
					}
				}
				else
				{
					if (line.Length > 0)
					{
						// not a blank line between contigs
						throw new Exception("Invalid file format for " + filename + " for line " + line);
					}
				}
			}

			if (srcFrom + size != srcEnd || (targetFrom + size) != targetEnd)
			{
				throw new Exception("Alignment blocks do not match block size");
			}
		}

		// similar to that of https://github.com/huangzhibo/CrossMap
		public List<Interval> mapCoordinates(Interval interval)
		{

			List<Interval> matches = new List<Interval>();

			if (!genomeMap.ContainsKey(interval.name))
			{
				Console.WriteLine(interval.name + " not found in genomeMap. returning none");
				return matches;
			}

			// overlaps query TODO is this correct
			List<Interval> overlappingIntervals = genomeMap[interval.name].Where(r => interval.start < r.end && interval.end > r.start).ToList();

			foreach (var overlap in overlappingIntervals)
			{
				Interval intersection = interval.Intersection(overlap);

				long offset = Math.Abs(intersection.start - overlap.start);
				long size = Math.Abs(intersection.end - intersection.start);

				if (overlap.strand == '+')
				{
					var thisStart = overlap.targetStart + offset;
					if (interval.strand == '+')
					{
						matches.Add(new Interval(overlap.targetName, thisStart, thisStart + size, overlap.strand));
					}
					else
					{
						matches.Add(new Interval(overlap.targetName, thisStart, thisStart + size, '-'));
					}
				} else if (overlap.strand == '-')
				{
					var thisStart = overlap.targetEnd - offset - size;
					if (interval.strand == '+')
					{
						matches.Add(new Interval(overlap.targetName, thisStart, thisStart + size, overlap.strand));
					}
					else
					{
						matches.Add(new Interval(overlap.targetName, thisStart, thisStart + size, '+'));
					}
				}

			}

			return matches;
		}

    }
}
