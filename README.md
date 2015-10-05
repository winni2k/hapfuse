# hapfuse

This is my own fork of the SNPTools::hapfuse tool in [SNPTools v1.0](http://sourceforge.net/projects/snptools).

## Differences to SNPTools::hapfuse

* **Dependency on htslib for parsing input chunks**.
  This allows me to avoid nasty buffer overruns when input VCFs contain too many samples. Also allows for the use of BCF files.

* **Input VCFs must contain the GT, and APP or GP field**.
  APP stands for phred-scaled allelic probability.  SNPTools uses the AP field, which stands for allelic probability.

* **Input VCFs can have fields with arbitrary precision**.  SNPTools::hapfuse uses the AP field with fixed precision floats while hapfuse v0.4 accepts APP and GP fields of arbitrary precision.  

* **Speed and Reliability**.
  Hapfuse is now way faster than SNPTools::hapfuse as it uses htslib for VCF/BCF parsing.  Hapfuse also concurrently loads, fuses and writes chunks, so expect to see CPU usage at about %200.

* **Regression Tests**. Hapfuse comes with regression tests to make sure it works!

* **Concurrency**. Hapfuse now uses concurrency to simultaneously load, merge and write VCF/BCF files.  This gives about a 2 fold speedup.

## WARNING
Hapfuse does not work with data generated by SNPTools::impute.  This
is because I do not use the AP field and am too lazy to implement
support for it.  Patches welcome.

## Installation

This code explains how to build hapfuse v1.1

    # clone the repo
    git clone git@bitbucket.org:wkretzsch/hapfuse.git
    cd hapfuse

    # configure and build hapfuse in hapfuse.1.1
    # the executable will be built at hapfuse.1.1/hapfuse
    ./bootstrap.pl

    # optionally, run regression tests
    cd hapfuse.1.1 && make check

### Static linking

I had to link boost and hts statically.  In order to do this, I ran this command instead of the final linking command:
    /bin/sh ./libtool --tag=CXX   --mode=link g++  -O3 -std=gnu++11 -L/usr/local/lib -Wl,-R,/usr/local/lib /users/flint/winni/lib/libboost_iostreams.a /users/flint/winni/lib/libhts.a /usr/lib64/libz.a /usr/lib64/libbz2.a -o hapfuse hapfuse.o hapSamp.o main.o utils.o -lpthread
libtool: link: g++ -O3 -std=gnu++11 -Wl,-R -Wl,/usr/local/lib -o hapfuse hapfuse.o hapSamp.o main.o utils.o  -L/usr/local/lib /users/flint/winni/lib/libboost_iostreams.a /users/flint/winni/lib/libhts.a /usr/lib64/libz.a /usr/lib64/libbz2.a -lpthread


## Usage

Fusing three bcf files together and save as gzip compressed bcf

    hapfuse file1.bcf file2.bcf file3.bcf -Ob -o fused.bcf.gz

Fusing three WTCCC style hap/sample file pairs together (see
https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html#hapsample)
and saving as gzip compressed BCF

    cat hap.files
    file1.hap.gz
    file2.hap.gz
    file3.hap.gz

    cat sample.files
    file1.sample
    file2.sample
    file3.sample
    
    # output as bcf
    hapfuse -h hap.files -s sample.files -Ob -o fused.bcf.gz
    
    # since version 1.0: output as WTCCC file pair
    hapfuse -h hap.files -s sample.files -Ow -o fused.wtccc.haps.gz,fused.wtccc.sample

## Description

Hapfuse ligates overlapping haplotype chunks with identical samples
by matching haplotypes at overlap sites.  The -o option is
required. Defaults are given in []. If VCF/BCF files are given on
the command line, then hapfuse assumes that the VCF/BCF files
contain haplotype chunks to ligate.  WTCCC style haplotype chunks
can be provided via the -h and -s arguments. The input and output
format tags are restricted to 'GT' when using the -h and -s
arguments.

### Ligation methods

The 'ligation method' is the method used to determine the diplotype
of a sample at overlapping chunk sites.

#### step

The 'step' method discards the outer half of a chunk's overlap
region haplotypes when determining the diplotype of an overlap
region.  This is the most commonly used method for ligating 
haplotypes.

#### average

The 'average' method uses the average of a sample's two overlapping
diplotypes to determine the new diplotypes. This is the method that
SNPTools::hapfuse uses, but it only makes sense when using GP or APP
fields.

#### linear

The 'linear' method is experimental.  It is the same as average, but
it gives more weight to diplotypes that are closer to a chunk's
center. 

## Thanks

I thank the following people for contributing code to this project:

- Joshua Randall for help with the autoconf files
