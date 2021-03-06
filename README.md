bedUniformLength
================

A genomics tool to make all regions in a bed file the same length.  The noGap.bed file is used to ensure
that regions are not extended into assembly gaps.  I often use this tool when you have a set of regions,
often of different lengths, but you want to apply a function to them that needs a fixed region size.  For
example: you want to predict enhancers and your prediction method is set to learn based on 1kb windows.
You now need to get a set of posible training examples, known enhancers, and make them all 1kb in length.

Installation
============

This code uses the Kent libraries from UCSC so that must be installed to compile bedUniformLength

<ol>
<li> Download and compile the Kent Libraries

<ol>
<li> Check the setting of your machtype variable with:<br />
echo $MACHTYPE<br />
it should be something in this list: i386 i686 sparc alpha x86_64 ppc.  If it is not, set your machtype variable.
<li> Go to a folder on your computer where you want the kent source tree to reside and type:<br />
git clone git://genome-source.cse.ucsc.edu/kent.git<br />
to download the repository onto your own computer.
<li> go to the src/lib directory within the kent source repo that you just cloned:<br />
cd kent/src/lib<br />
<li> Compile the libraries<br />
make
<li> If this was successful, you should have a file here:<br />
kent/src/lib/x86_64/jkweb.a<br />
the x86_64 will be the machtype of your machine.</br />
</ol>

If this was not successful then you should look at the build instructions in the kent repo itself
by looking at this file:<br />
kent/src/product/README.building.source

<li> Compile bedUniformLength
<ol>
<li> Edit the bedUniformLength makefile so that it points to the kent libraries on your system.  These
are the two lines you will have to modify:<br />
HG_INC += -I/home/lowec/kent/src/hg/inc -I/home/lowec/kent/src/inc<br />
L += /home/lowec/kent/src/lib/${MACHTYPE}/jkweb.a<br />

<li> Compile bedUniformLength:<br />
make

<li> Running bedUniformLength with no parameters should give a brief help message.
</ol>
</ol>

References
==========

This code has not yet been described in a publication.

