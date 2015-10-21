/*

bedUniformLength.c
Written by: Craig Lowe

*/

#include "common.h"
#include "linefile.h"
#include "options.h"
#include "memalloc.h"
#include "bed.h"

/*---------------------------------------------------------------------------*/

static struct optionSpec optionSpecs[] =
/* command line option specifications */
{
	{"length", OPTION_INT},
	{NULL, 0}
};

int optLength = 1000;

/*---------------------------------------------------------------------------*/

void usage()
{
errAbort(
	"bedUniformLength - expand or shrink all regions in the bed file\n"
	"                   to be a uniform length.\n"
	"usage:\n"
	"   bedUniformLength in.bed noGap.bed randomlyPlaced.bed\n"
	"options:\n"
	"    -length    (1000)     new length of all the elements\n"
	"notes:\n"
	"    The noGap.bed file is a bed file of all the ungapped regions\n"
	"     in the genome so that regions are not expanded into gaps.\n"
	);
}

/*---------------------------------------------------------------------------*/

unsigned int maxUnsigned(unsigned int a, unsigned int b)
{
	if(a>=b){return(a);}
	else{return(b);}
}

unsigned int minUnsigned(unsigned int a, unsigned int b)
{
	if(a<=b){return(a);}
	else{return(b);}
}

boolean bedOverlap(struct bed *a, struct bed *b)
{
	return(sameString(a->chrom, b->chrom) && maxUnsigned(a->chromStart, b->chromStart) < minUnsigned(a->chromEnd, b->chromEnd));
}

void adjustBed(struct bed *a, struct bed *noGap, unsigned int newLength)
{
	unsigned int currLength = 0, newStart = 0;
	currLength = a->chromEnd - a->chromStart;
	
	if(noGap->chromEnd - noGap->chromStart < newLength){fprintf(stderr, "Warning: not enough room to expand region: %s %u %u\n", a->chrom, a->chromStart, a->chromEnd);}

	if(currLength > newLength)
	{
		newStart = maxUnsigned(a->chromStart + (currLength - newLength) / 2, noGap->chromStart);
	}
	else if(currLength < newLength)
	{
		
		newStart = maxUnsigned(a->chromStart - (newLength - currLength) / 2, noGap->chromStart);
	}

	if(newStart + newLength <= noGap->chromEnd)	
	{
		a->chromStart = newStart;
		a->chromEnd = a->chromStart + newLength;
	}
	else
	{
		a->chromEnd = noGap->chromEnd;
		a->chromStart = a->chromEnd - newLength;
	}
}

/*---------------------------------------------------------------------------*/

void bedUniformLength(char *bedFile, char *noGapBedFile, char *outFilename)
{
	struct bed *regionBedList = NULL, *noGapBedList = NULL, *currRegion = NULL, *currNoGap = NULL;

	noGapBedList = bedLoadNAll(noGapBedFile, 3);
	slSort(&noGapBedList, bedCmp);
	regionBedList = bedLoadNAll(bedFile, 4);
	slSort(&regionBedList, bedCmp);

	FILE *fout = mustOpen(outFilename, "w");
	for(currRegion=regionBedList, currNoGap=noGapBedList; currRegion!=NULL; currRegion=currRegion->next)
	{
		while(!bedOverlap(currRegion, currNoGap))
		{
			currNoGap=currNoGap->next;
			if(currNoGap == NULL){errAbort("No overlap found: %s\n", currRegion->chrom);}
		}
		adjustBed(currRegion, currNoGap, (unsigned int)optLength);
		bedTabOut(currRegion, fout);
	}
	carefulClose(&fout);
	if(currRegion != NULL){errAbort("Not all regions overlapped an ungapped region\n");}
}

/*---------------------------------------------------------------------------*/


int main(int argc, char *argv[])
{
	optionInit(&argc, argv, optionSpecs);
	if (argc != 4){usage();}

	optLength = optionInt("length", optLength);
	if(optLength < 1){errAbort("Error: length can not be less than 1\n");}

	bedUniformLength(argv[1],argv[2],argv[3]);
	return(0);
}

