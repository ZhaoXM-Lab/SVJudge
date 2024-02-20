import os

import pandas as pd

import src.annot_func as annot_func
import src.utils as utils


# %%
def preprocess_judge_sv_annot_on_gene(input_sv_bed_path, gene_reference_bed, gene_annotation_dict,
                                      split_gene_public_sv_af_file,
                                      output_dir, overlap_threshold, distance_threshold, out_prefix=''):
    if out_prefix == '':
        out_prefix = 'SVJudge'
    # print('Annote SVs in %s (Gene Structures)...' % out_prefix)
    sv_gene_intersection_bed = annot_func.sv_gene_intersection(sv_bed=input_sv_bed_path,
                                                               gene_annot_bed=gene_reference_bed,
                                                               output_dir=output_dir,
                                                               prefix=out_prefix)
    whole_gene_sv_annot_file = annot_func.annotate_sv_on_whole_gene(sv_gene_intersection_bed,
                                                                    gene_annot_dict=gene_annotation_dict,
                                                                    output_dir=output_dir,
                                                                    file_prefix=out_prefix)
    out_file_whole_gene_annot = os.path.join(output_dir, out_prefix + '.SVJudge.whole_genes.annotated.tsv')
    out_file_gene_split_annot = os.path.join(output_dir, out_prefix + '.SVJudge.gene_elements.annotated.tsv')
    # print('Split SVs in ' + out_prefix)
    split_gene_sv_annot_file = annot_func.split_gene_elements(whole_gene_sv_annot_file=whole_gene_sv_annot_file,
                                                              gene_annot_dict=gene_annotation_dict)
    # print('Add AF for SVs in ' + out_prefix)
    df_split_gene_judge_sv = pd.read_csv(split_gene_sv_annot_file, sep='\t')
    df_split_gene_judge_sv_add_af = annot_func.annotate_judge_sv_af_in_split_gene(
        df_split_gene_judge_sv=df_split_gene_judge_sv,
        split_gene_public_sv_af_file=split_gene_public_sv_af_file,
        overlap_threshold=overlap_threshold,
        distance_threshold=distance_threshold)
    df_whole_gene_judge_sv = pd.read_csv(whole_gene_sv_annot_file, sep='\t')
    df_split_gene_judge_sv_add_af.to_csv(out_file_gene_split_annot, sep='\t', index=False)
    df_whole_gene_judge_sv.to_csv(out_file_whole_gene_annot, sep='\t', index=False)
    os.system('rm %s' % whole_gene_sv_annot_file)
    return out_file_gene_split_annot, out_file_whole_gene_annot


# %%
def preprocess_judge_sv_annot_on_regulation(input_sv_bed_path, regulatory_reference_beds, regulatory_public_sv_af_file,
                                            regulatory_element_types, output_dir,
                                            overlap_threshold, distance_threshold, output_prefix=''):
    """
    Preprocess and annotate structural variants (SVs) related to regulatory elements and add allele frequency (AF) data.

    :param input_sv_bed_path: Path to the input SV BED file.
    :param regulatory_reference_beds: Paths to the regulatory reference BED files.
    :param regulatory_public_sv_af_file: File containing AF data for public SVs.
    :param regulatory_element_types: Types of regulatory elements to consider.
    :param output_dir: Directory for output files.
    :param overlap_threshold: Threshold for determining SV overlap.
    :param distance_threshold: Distance threshold for SV comparison.
    :param output_prefix: Prefix for output file names. Defaults to 'SVJudge'.
    :return: Path to the final annotated file for regulatory elements.
    """
    if output_prefix == 'SVJudge':
        output_prefix = 'SVJudge'
    out_file_regulatory_annot = os.path.join(output_dir, output_prefix + '.SVJudge.regulatory_elements.annotated.tsv')
    # print('Annote SVs in %s (Regulation elements)...' % output_prefix)
    sv_regulation_intersection_bed = annot_func.sv_regulatory_elements_intersection(sv_bed_path=input_sv_bed_path,
                                                                                    regulatory_reference_beds=regulatory_reference_beds,
                                                                                    output_directory=output_dir,
                                                                                    file_prefix=output_prefix)
    regulation_sv_annot_file = annot_func.annotate_sv_on_regulatory_element(
        sv_regulation_intersection_bed=sv_regulation_intersection_bed,
        regulatory_element_types=regulatory_element_types,
        output_directory=output_dir, file_prefix=output_prefix)
    # print('Add AF for SVs in ' + output_prefix)
    df_regulatory_judge_sv = pd.read_csv(regulation_sv_annot_file, sep='\t')
    df_regulatory_judge_sv_add_af = annot_func.annotate_judge_sv_af_in_regulatory(
        df_regulatory_judge_sv=df_regulatory_judge_sv,
        regulatory_public_sv_af_file=regulatory_public_sv_af_file,
        overlap_threshold=overlap_threshold,
        distance_threshold=distance_threshold,
        column_prefix='')
    df_regulatory_judge_sv_add_af.to_csv(out_file_regulatory_annot, sep='\t', index=False)
    return out_file_regulatory_annot


# %%
def run_annotation_only_judge_sv(output_dir, annotation_tuple, input_sv_bed_path, split_gene_public_sv_af_file,
                                 regulatory_public_sv_af_file,
                                 overlap_threshold, distance_threshold, output_prefix=''):
    """
    Annotate judge SV files with gene and regulatory elements information and add allele frequency (AF) data.

    :param output_dir: Directory for output files.
    :param annotation_tuple: tuple of annotation data.
    :param input_sv_bed_path: Path to the input SV BED file.
    :param split_gene_public_sv_af_file: File containing AF data for public SVs related to gene elements.
    :param regulatory_public_sv_af_file: File containing AF data for public SVs related to regulatory elements.
    :param overlap_threshold: Threshold for determining SV overlap.
    :param distance_threshold: Distance threshold for SV comparison.
    :param output_prefix: Prefix for output file names. Defaults to empty.
    :return: Paths to the annotated files for gene split, whole gene, and regulatory elements.
    """
    # os.system("mkdir -p %s" % output_dir)
    # log_file = open(os.path.join(output_dir, 'RUN_SV_Annotation.log'), 'w')
    # start_time = time.time()
    # logging.info('Loading reference and annotation files...')
    # time_point1, log_str = utils.calculate_runtimes(log_str='Loading reference and annotation files...',
    #                                                start_time=start_time)
    gene_annotation_dict, gene_reference_bed, regulatory_reference_beds, regulatory_element_types = annotation_tuple
    # logging.info(f'The types of regulatory elements are {regulatory_element_types}')
    # time_point2, log_str = utils.calculate_runtimes(log_str='Finish Loading', start_time=time_point1)
    # print_str = '使用的公共SV sets：%s' % selected_datasets
    # print(print_str)
    # log_file.write(log_str + '\n')
    out_file_split_gene_annot, out_file_whole_gene_annot = preprocess_judge_sv_annot_on_gene(
        input_sv_bed_path=input_sv_bed_path,
        gene_reference_bed=gene_reference_bed,
        gene_annotation_dict=gene_annotation_dict,
        split_gene_public_sv_af_file=split_gene_public_sv_af_file,
        output_dir=output_dir,
        overlap_threshold=overlap_threshold,
        distance_threshold=overlap_threshold, out_prefix=output_prefix)

    out_file_regulatory_annot = preprocess_judge_sv_annot_on_regulation(
        input_sv_bed_path=input_sv_bed_path, regulatory_reference_beds=regulatory_reference_beds,
        regulatory_public_sv_af_file=regulatory_public_sv_af_file,
        regulatory_element_types=regulatory_element_types,
        output_dir=output_dir,
        overlap_threshold=overlap_threshold, distance_threshold=distance_threshold, output_prefix=output_prefix)
    return out_file_split_gene_annot, out_file_whole_gene_annot, out_file_regulatory_annot
