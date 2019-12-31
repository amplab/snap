// A c++ version of vcffirstheader, because the script from vcftools is much slower than the rest of the pipeline for no apparent reason.

#include "stdio.h"
#include "stdlib.h"
#include "string.h"

int main()
{
	int bufferSize = 1024;
	char *buffer = new char[bufferSize];
	int usedBuffer = 0;

	bool inHeader = true;

	while (NULL != fgets(buffer + usedBuffer, bufferSize - usedBuffer - 1, stdin)) 
	{
		if (buffer[0] == '\0')
		{
			fprintf(stderr, "fgets returned the empty string\n");
			exit(1);
		}

		usedBuffer = (int)strlen(buffer);

		if (buffer[usedBuffer - 1 != '\n'])
		{
			//
			// Overran the buffer.  Double it, copy in what we've already read and continue.
			//
			bufferSize *= 2;
			char *newBuffer = new char[bufferSize];
			memcpy(newBuffer, buffer, usedBuffer);
			delete[] buffer;
			buffer = newBuffer;
			continue;
		}


		if (buffer[0] == '#') 
		{
			if (!inHeader) 
			{
				usedBuffer = 0;
				continue;
			}

			if (buffer[1] != '#') 
			{
				//
				// The last line of the header has a single #.  Once we've seen it, we're done.  But emit the first copy.
				//
				inHeader = false;
			}



		} // while we have an input line
		fwrite(buffer, 1, usedBuffer, stdout);
		usedBuffer = 0;
	}

}
