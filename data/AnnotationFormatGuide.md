

# Annotations Format Guide

### Gene Annotation Structure: `/gene_annotation/gene_structure/refGene.txt`

This file originates from the Gene Transfer Format (GTF) and delineates the structure of gene annotations. The data
columns are separated by a tab (`\t`) and include the following fields:

- `txID`: Transcript ID
- `Chr`: Chromosome
- `Strand`: Strand direction (`+` for positive, `-` for negative)
- `txStart`: Transcription start position
- `txEnd`: Transcription end position
- `cdsStart`: Coding sequence start position
- `cdsEnd`: Coding sequence end position
- `exonCount`: Number of exons
- `exonStarts`: Exon start positions (comma-separated)
- `exonEnds`: Exon end positions (comma-separated)
- `Score`: A score related to the annotation (specific interpretation may vary)
- `geneID`: Gene ID
- `cdsStartStat`: Status of the CDS start (e.g., known, unknown)
- `cdsEndStat`: Status of the CDS end (e.g., known, unknown)
- `exonFrames`: Frame of each exon (comma-separated)

### Regulatory Annotation Structure: `/regulatory_element/{Promoter, Enhancer, DNA_Accessible, TAD}/*.bed`

These files are in the BED format, which is widely used for representing genomic regions. The columns in these files are
delineated by a tab (`\t`) and encompass the following information:

- `Chr`: Chromosome on which the regulatory element is located.
- `Start`: Start position of the regulatory element on the chromosome.
- `End`: End position of the regulatory element on the chromosome.
- `MapGeneName`: Name of the gene associated with this regulatory element. Note that TADs do not specify associated
  genes, therefore, use `TAD_boundary` to populate this field.
- `ID`: Unique identifier for the regulatory element.

Users can generate a new annotation file and replace the original BED file with it.

### Structural Variant Annotation Structure: `/structural_variant/*.bed`

This section details the format of the BED files used for structural variant annotations. The columns in these files are
tab-separated (`\t`) and provide critical information for understanding structural variations in the genome:

- `Chr`: The chromosome where the structural variant is located.
- `Start`: The start position of the structural variant on the chromosome.
- `End`: The end position of the structural variant on the chromosome.
- `SV_ID`: A unique identifier assigned to the structural variant.
- `SVTYPE`: The type of structural variant (e.g., deletion, duplication, inversion).
- `SVLEN`: Length of the structural variant in base pairs.
- `AF`: Allele frequency, indicating how common the variant is in a given population.

To add new public structural variant (SV) sets to the database, it's crucial to ensure that each set name is accurately
paired with its respective BED file path. This can be achieved using two command-line parameters: `-add_sv_set_name`
and `-add_sv_set_dir`. These parameters can be used multiple times within the same command, and they are paired in the
order they appear.

For example, to add two SV sets, you would use the following command structure:

```-add_sv_set_name SET1 -add_sv_set_dir /path/to/dir/SET1.bed -add_sv_set_name SET2 -add_sv_set_dir /path/to/dir/SET2.bed```


In this command:

- `-add_sv_set_name` is used to specify the name of the SV set.
- `-add_sv_set_dir` is used to provide the path to the BED file corresponding to that SV set.

The command ensures that `SET1` is associated with `/path/to/dir/SET1.bed` and `SET2` with `/path/to/dir/SET2.bed`, thereby maintaining a clear and accurate mapping between each SV set name and its file path.
