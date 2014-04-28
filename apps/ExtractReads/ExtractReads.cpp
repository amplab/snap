/*++

Module Name:

    ComputeROC.cpp

Abstract:

   Take a SAM file with simulated reads and compute a ROC curve from it.

Authors:

    Bill Bolosky, December, 2012

Environment:
`
    User mode service.

Revision History:

    
--*/

#include "stdafx.h"
#include "SAM.h"
#include "BAM.h"
#include "Genome.h"
#include "Compat.h"
#include "Read.h"
#include "RangeSplitter.h"
#include "BigAlloc.h"

void usage()
{
    fprintf(stderr,"usage: ExtractReads genomeDirectory input.bam output.sam chromosome\n");
  	exit(1);
}


int main(int argc, char * argv[])
{
    BigAllocUseHugePages = false;

    if (argc != 5) usage();

    static const char *genomeSuffix = "Genome";
	size_t filenameLen = strlen(argv[1]) + 1 + strlen(genomeSuffix) + 1;
	char *fileName = new char[strlen(argv[1]) + 1 + strlen(genomeSuffix) + 1];
	snprintf(fileName,filenameLen,"%s%c%s",argv[1],PATH_SEP,genomeSuffix);
	const Genome *genome = Genome::loadFromFile(fileName, 0);
	if (NULL == genome) {
		fprintf(stderr,"Unable to load genome from file '%s'\n",fileName);
		return -1;
	}
	delete [] fileName;
	fileName = NULL;

    ReaderContext readerContext;
    readerContext.header = NULL;
    readerContext.genome = genome;
    readerContext.clipping = NoClipping;
    readerContext.paired = false;
    readerContext.ignoreSecondaryAlignments = true;
    readerContext.ignoreSupplementaryAlignments = true;
    readerContext.defaultReadGroup = "";

    ReadSupplierGenerator *readSupplierGenerator = BAMReader::createReadSupplierGenerator(fileName,1, readerContext);
    ReadSupplier *readSupplier = readSupplierGenerator->generateNewReadSupplier();

    AlignerOptions options("");
    options.outputFileTemplate = argv[4];
    options.sortOutput = false;
    
    const FileFormat* format = 
    FileFormat::SAM[0]->isFormatOf(argv[3]) ? FileFormat::SAM[false] :
    FileFormat::BAM[0]->isFormatOf(argv[3]) ? FileFormat::BAM[false] :
    NULL;

    if (NULL == format) {
        fprintf(stderr,"Can't determine format for output file (does it end in .sam or .bam?)\n");
        return 0;
    }

    ReadWriterSupplier *writerSupplier = format->getWriterSupplier(&options, readerContext.genome);
    ReadWriter* writer = writerSupplier->getWriter();
    writer->writeHeader(readerContext, options.sortOutput, argc, (const char **)argv, "", "");

    Read *read;
    AlignmentResult alignmentResult;
    unsigned genomeLocation;
    bool isRC;
    unsigned mapQ;
    _int64 totalReads = 0;
    _int64 emittedReads = 0;
    while (readSupplier->getNextRead()) {
        totalReads++;
        const Genome::Contig *contig = genome->getContigAtLocation(genomeLocation);
        if (NULL != contig && !strcmp(contig->name, argv[4])) {
            emittedReads++;
            writer->writeRead(read, alignmentResult, mapQ, genomeLocation, isRC ? RC : FORWARD);
        }
    }
    writer->close();

    printf("Processed %lld reads, of which %lld were emitted\n", totalReads, emittedReads);

	return 0;
}

