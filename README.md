## RCSeq

Contact: adam.ewing@mater.uq.edu.au

*Analyse RC-seq data, or in principle any sequence data, for transposable element insertions.*

# Prerequisites:

|what     | where | why |
|---------|-------|-----|
|bwa      | http://bio-bwa.sourceforge.net/  | sequence alignments |
|samtools | http://samtools.sourceforge.net/ | BAM manipulation |
|pysam    | https://github.com/pysam-developers/pysam | parsing SAM/BAM formatted alignments |
|FLASH    | http://ccb.jhu.edu/software/FLASH/ | not needed if not using overlapped libraries, used for assembling read pairs |
|MAFFT    | http://mafft.cbrc.jp/alignment/software/ | needed if using the consensus sequence option |

(latest versions unless otherwise specified)
