# SVJudge

We are glad to introduce SVJudge, an innovative method specifically designed to evaluate Structural Variations (SVs) in
the context of their potential association with diseases. Recognizing the inherent challenges posed by the imprecision
of SV breakpoints and their inherent complexity, SVJudge adopts a targeted approach. Instead of analyzing the full
extent of SVs, our method focuses on their intersections with critical genomic functional areas. These include Coding
Sequences (CDS)[1], Untranslated Regions (UTRs), promoters, enhancers, introns, DNA accessible regions, and
Topologically Associating Domain (TAD) boundaries[2].

Our rationale is straightforward yet powerful: if segments of SVs that overlap with these key functional elements are
frequently altered in the general population, it diminishes the likelihood of these SVs being pathogenic. Building on
this concept, we have developed a robust scoring system. This system evaluates the pathogenic potential of SVs by
considering both the prevalence of SVs in the population and the significance of the affected functional genomic
regions.

To ensure the effectiveness and reliability of SVJudge, it is essential to utilize public datasets that satisfy specific
criteria. These include a substantial representation of normal individuals and the provision of allele frequencies.
After careful consideration, we selected comprehensive public SV datasets from sources such as gnomAD[^1],
1000PG_phase3[^2], DGV[^3], Abel et al.[^4], and PGG.SV[^5].

Originally, SVJudge was crafted with a focus on analyzing Structural Variations (SVs) in schizophrenia patients,
utilizing the GRCh37 reference genome. In this context, we prioritized tissue-specific annotations of regulatory
elements pertinent to brain tissue[6,7,8]. For researchers aiming to apply SVJudge in predicting SVs in different
contexts, it
is crucial to customize the annotation file to suit your specific research needs. We provide a comprehensive guide on
creating an appropriate annotation file in [Annotation Format Guide](data/AnnotationFormatGuide.md). By following this guide, researchers can
generate a new, tailored annotation file that aligns with their specific study requirements and seamlessly integrate it
with SVJudge.

## Installation and Usage
### 1. Configure a Conda environment

```shell
conda env create -f environment.yml
```

### 2. Download SVJudge

```
git clone https://github.com/ZhaoXM-Lab/SVJudge
```
### 3. Navigate to the `path/to/SVJudge` and run `SVJudge.py`.
```
cd path/to/SVJudge

python ./SVJudge.py -i ./exmple_input_svs.vcf -o ./test -r GRCh37 -prefix Test`
```
### Required Input/Ouput parameters
```
  -o OUTPUT_PATH        Absolute path to output
  -i INPUT_VCF          Absolute path to vcf file
  -r {GRCh37,GRCh38}    Absolute path to your reference genome (GRCh37 or
                        GRCh38)
  -prefix PREFIX        The prefix of result file (default: SVJudge)
```
For other optional arguments, please use `-h` to view them.
## References
[1]Collins RL, Brand H, Karczewski KJ, Zhao X, Alföldi J, Francioli LC, et al. A structural variation reference for
medical and population genetics. Nature. 2020;581(7809):444-51.

[2]Sudmant PH, Rausch T, Gardner EJ, Handsaker RE, Abyzov A, Huddleston J, et al. An integrated map of structural
variation in 2,504 human genomes. Nature. 2015;526(7571):75-81.

[3]Lappalainen I, Lopez J, Skipper L, Hefferon T, Spalding JD, Garner J, et al. DbVar and DGVa: public archives for
genomic structural variation. Nucleic acids research. 2012;41(D1):D936-D41.

[4]Abel HJ, Larson DE, Regier AA, Chiang C, Das I, Kanchi KL, et al. Mapping and characterization of structural
variation in 17,795 human genomes. Nature. 2020;583(7814):83-9.

[5]Wang Y, Ling Y, Gong J, Zhao X, Zhou H, Xie B, et al. PGG.SV: a whole-genome-sequencing-based structural variant
resource and data analysis platform. Nucleic Acids Research. 2023;51(D1):D1109-D16.

[6]Zhao X, Song L, Yang A, Zhang Z, Zhang J,Yang YT, et al. Prioritizing genes associated with brain disorders
byleveraging enhancer-promoter interactions in diverse neural cells and tissues.Genome Medicine. 2023;15(1):56.

[7]Akbarian S, Liu C, Knowles JA, Vaccarino FM,Farnham PJ, Crawford GE, et al. The PsychENCODE project. Nature
Neuroscience.2015;18(12):1707-12.

[8]Fullard JF, Hauberg ME, BendlJ, Egervari G, Cirnaru MD, Reach SM, et al. An atlas of chromatin accessibility
in the adult human brain. Genome Res. 2018;28(8):1243-52.
