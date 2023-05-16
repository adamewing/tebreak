## TEBreak 

[![DOI](https://zenodo.org/badge/17894795.svg)](https://zenodo.org/badge/latestdoi/17894795)

Contact: adam.ewing@mater.uq.edu.au

*Tools for analysing insertion mutations*

# Installation

## via conda (or [mamba](https://anaconda.org/conda-forge/mamba)):

```
git clone https://github.com/adamewing/tebreak.git
cd tebreak
conda env create -f tebreak.yml
conda activate tebreak
python setup.py install
cd test && ./test.sh
```
If you use the above method, make sure to activate the Conda environment first with `conda activate tebreak` whenever using tebreak.

## Python libraries (if the conda solution is not an option):
This assumes a working installation of `pip`. Many of these prerequisites can be satisfied through installing [anaconda](https://conda.io/docs/user-guide/install/download.html).

```
pip install pysam
pip install scipy
pip install bx-python
pip install scikit-bio
```

If `pip install bx-python` fails you might need `liblzo2-dev` (via apt: `sudo apt-get install -y liblzo2-dev`).

## LAST (aligner)
```
wget http://last.cbrc.jp/last-716.zip
unzip last-716.zip
make CXXFLAGS=-O3 -C last-716 && sudo make install -C last-716
```

## HTSLIB / SAMtools / BCFtools
```
git clone https://github.com/samtools/htslib.git
git clone https://github.com/samtools/samtools.git
git clone https://github.com/samtools/bcftools.git

make -C htslib && sudo make install -C htslib
make -C samtools && sudo make install -C samtools
make -C bcftools && sudo make install -C bcftools
```

## Minia (sequence assembler)
```
https://github.com/GATB/minia
```

## Exonerate (aligner)
```
git clone https://github.com/adamewing/exonerate.git
cd exonerate
git checkout v2.4.0
autoreconf -i
./configure && make && make check && make install
```
# Install
```
python setup.py install
```

# Test your installation
Assuming `$TB` is the tebreak directory created by `git clone` or unzipping/untarballing an archive:

```
tebreak -b $TB/test/data/example.ins.bam -r $TB/test/data/Homo_sapiens_chr4_50000000-60000000_assembly19.fasta -i $TB/lib/teref.human.fa

```
or
```
cd test && ./test.sh
```

This will generate some output to the terminal and the following files should exist in your working directory:

|filename                         | description |
|---------------------------------|-------------|
|`example.ins.tebreak.detail.out` | Details on all potential insertions detected (probably not a useful final output, used for debugging) |
|`example.ins.tebreak.pickle`     | Raw data on detected insertions. Allows trying multiple parameters via `--use_pickle` without needing to re-run completely. |
|`example.ins.tebreak.resolve.out`| Details on all potential insertions considered (probably not a useful final output, used for debugging) |
|`example.ins.tebreak.table.txt`  | Final output table. Often requires further filtering. |

The file `example.ins.tebreak.table.txt` should contain five insertions.

# Running TEBreak on real data

The parameters for the test run are the bare minimum required to run TEBreak and will be glacially slow on anything larger than the most trivial input. The following is the current recommendation for running TEBreak on WGS data with an average depth > 30x and should also suffice for capture-seq data.

## Generate BAM file(s)
Our recommendation is to use `bwa mem` with the following parameters. Let `$THREADS` be the number of CPU cores available on the system, `$RGID` be a read group id, `$SM` be a sample name, `$RAWBAM` be the output BAM filename, `$REF` be a bwa-indexed reference genome, `$FQ1` and `$FQ2` be .fastq files containing read 1 and read 2, respectively.


```
bwa mem -M -Y -t $THREADS -R "@RG\tID:$RGID\tSM:$BASE\tPL:ILLUMINA" $REF $FQ1 $FQ2 | samtools view -b - > $RAWBAM
```

## Mark Duplicate Reads
There's more than one way to do this, one option is to use [picard](https://broadinstitute.github.io/picard/). Let `$RAWBAM` be from the original alignment and `$BAM` be the BAM file used in subsequent steps.

```
java -jar picard.jar MarkDuplicates I=$RAWBAM O=$BAM M=metrics.out
```

Another good option is to use [samblaster](https://github.com/GregoryFaust/samblaster).

## Build the relevant reference file(s)
This example assumes hg19/GRCh37 without the 'chr' prefix). The following builds a reference containing the locations of relevant human repeatmasker annotations for discordant read pair discovery:
```
cd $TB/lib
./make_discref_hg19.sh
```

Additionally, it may be helpful to build a mappability index using the `./human_mappability.sh` script, but it is not required for this example.

## Run tebreak
Note that the BAM file (`$BAM`) passed to -b can be a comma delimited list of BAM files or a `.txt` file containing a list of BAM files.

```
tebreak -b $BAM -r $REF -p $THREADS -d $TB/lib/hg19.te.disctgt.txt -m $TB/lib/hg19.centromere_telomere.bed --max_ins_reads 500 -i $TB/lib/teref.human.fa 
```

## Filter the output (optional and depreciated - see "Filter" column)
The results table (`$TABLE`) will contain false positives. If desired, it is possible reduce this with an included script at some cost in terms of sensitivity.

```
$TB/scripts/general_filter.py -t $TABLE -i $TB/lib/teref.human.fa -r $REF --numsplit 4 --numdiscord 4 > $FILTEREDTABLE
```

## Values in the Filter column

|Filter         | Description                                                                                                       |
|---------------|-------------------------------------------------------------------------------------------------------------------|
|UnknownInsType | Combination of Superfamily and Subfamily does not appear in insertion library                                     |
|NoConsMapRef   | No insertion consensus sequences match the insertion reference                                                    |
|TotalConsLen   | Total consensus length is less than minimum set by `--min_cons_len`                                               |
|MinEltMatch    | Best match to reference element is less than `--min_ins_match`                                                    |
|MinRefMatch    | Best match to reference genome is less than `--min_ref_match`                                                     |
|MinDiscord     | Number of discordant reads is less than `--min_disc_reads`                                                        |
|MinSplit       | Number of split reads is less than `--min_split_reads`                                                            |
|MismatchTSD    | 5-prime and 3-prime TSD sequences do not match                                                                    |
|LongHomopolTSD | TSD is a long homopolymer (greater than 10 bases)                                                                 |
|MinTELength    | Insertion is shorter than `--minlength`                                                                           |
|MissingVAF     | Genotyping not possible (no TSD), only occurs with `--minvaf` set                                                 |
|MinVAF         | Maximum VAF less than `--minvaf` (only occurs if set)                                                             |
|FracEnd        | 3-prime end of TE alignment not not occur within `--fracend` percent of ref TE end (only occurs with `--fracend`) |
|MaxVars        | More than `--maxvars` variants versus TE reference (only occurs if set)                                           |
|LowMap         | Mappability low (less than 0.5), only occurs if `--map_tabix` is set                                              |
|HomopolSite    | Insertion site contains a homopolymer greater than 20bp                                                           |
|SelfAlign      | Reference genome has a better than 95% match to reference TE via exonerate                                        |
|PASS           | Good to go!                                                                                                       |

## Annotate the output (optional)
Finally, a script is included to annotate the TEBreak table. A useful included annotation source is the list of known non-reference insertions detected in human (hg19/GRCh37 coordinates).

```
$TB/scripts/annotate.py -t $FILTEREDTABLE -x $TB/lib/nonref.collection.hg19.bed.gz -n KnownNonRef --nonref > $FINALTABLE
```

## Getting help

Reporting [issues](https://github.com/adamewing/tebreak/issues) and questions through github is preferred versus e-mail.

For additional documentation, please find the manual in the [doc](https://github.com/adamewing/tebreak/tree/master/doc) subdirectory.
