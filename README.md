Utility to create a GEX equivalent to atac_fragments.tsv. This makes it easy to perform pileups

== utility ==

java -jar 10xrna2fragments.jar gex_possorted_bam.bam gex_fragments.tsv

Sort it: sort -k 1,1 -k2,2n gex_fragments.tsv > gex_fragments.sorted.tsv
Compress it: bgzip -@ 8 gex_fragments.sorted.tsv
Index it: tabix -p vcf gex_fragments.sorted.tsv.gz
