import logging
import os

import src.annot_func as annot_func


# %%
def preprocess_public_sv_annot_on_gene(input_sv_bed_paths, selected_datasets, gene_reference_bed,
                                       gene_annotation_dict, output_dir,
                                       overlap_threshold, distance_threshold):
    """
    Annotate allele frequencies (AFs) of structural variants (SVs) on gene elements using public SV datasets.

    Args:
        input_sv_bed_paths (dict): Paths to public SV datasets.
        selected_datasets (list): Dataset keys to be used.
        gene_reference_bed (str): Path to gene reference BED file.
        gene_annotation_dict (dict): Dictionary of canonical gene annotations.
        output_dir (str): Directory for output files.
        overlap_threshold (float): Threshold for SV overlap.
        distance_threshold (int): Threshold for breakpoint distance.

    Returns:
        tuple: Path to merged gene-annotated SV file, dictionary of split gene annotations.
    """
    intersection_bed_dict = {}
    whole_gene_sv_annot_dict = {}

    for data_name in selected_datasets:
        if data_name not in selected_datasets:
            continue
        logging.info('Annote SVs in %s (Gene Structures)...' % data_name)
        sv_gene_intersection_bed = annot_func.sv_gene_intersection(sv_bed=input_sv_bed_paths[data_name],
                                                                   gene_annot_bed=gene_reference_bed,
                                                                   output_dir=output_dir,
                                                                   prefix=data_name)
        intersection_bed_dict[data_name] = sv_gene_intersection_bed
        whole_gene_sv_annot_file = annot_func.annotate_sv_on_whole_gene(sv_gene_intersection_bed,
                                                                        gene_annot_dict=gene_annotation_dict,
                                                                        output_dir=output_dir,
                                                                        file_prefix=data_name)
        whole_gene_sv_annot_dict[data_name] = whole_gene_sv_annot_file
    split_gene_sv_annot_dict = {}
    for data_name in whole_gene_sv_annot_dict:
        logging.info('Split SVs in ' + data_name)
        split_gene_sv_annot_dict[data_name] = annot_func.split_gene_elements(
            whole_gene_sv_annot_file=whole_gene_sv_annot_dict[data_name],
            gene_annot_dict=gene_annotation_dict)

    split_gene_sv_annot_merge_sets_file = os.path.join(output_dir,'Public_SV_sets.Merge.gene_annot.split.tsv')
    logging.info('Annote Gene Element by AF...')
    annot_func.merge_filter_sv_sets_in_gene(input_sv_af_dict=split_gene_sv_annot_dict,
                                            output_file=split_gene_sv_annot_merge_sets_file,
                                            overlap_threshold=overlap_threshold, distance_threshold=distance_threshold)
    return split_gene_sv_annot_merge_sets_file, split_gene_sv_annot_dict, whole_gene_sv_annot_dict


# %%
def preprocess_public_sv_annot_on_regulation(input_sv_bed_paths, selected_datasets, regulatory_reference_beds,
                                             regulatory_element_types, output_dir,
                                             overlap_threshold, distance_threshold):
    """
    Annotate allele frequencies (AFs) of structural variants (SVs) on regulatory elements using public SV datasets.

    Args:
        input_sv_bed_paths (dict): Paths to public SV datasets.
        selected_datasets (list): Dataset keys to be used.
        regulatory_reference_beds (list): Paths to regulatory element annotation files.
        regulatory_element_types (list): List of regulatory element type keys.
        output_dir (str): Directory for output files.
        overlap_threshold (float): Threshold for SV overlap.
        distance_threshold (int): Threshold for breakpoint distance.
    Returns:
        tuple: Path to merged regulatory element annotated SV file, dictionary of regulatory element annotations.
    """
    intersection_bed_dict = {}
    regulation_sv_annot_dict = {}
    for data_name in input_sv_bed_paths:
        if data_name not in selected_datasets:
            continue
        logging.info('Annote SVs in %s (Regulation elements)...' % data_name)
        sv_regulation_intersection_bed = annot_func.sv_regulatory_elements_intersection(
            sv_bed_path=input_sv_bed_paths[data_name]
            , regulatory_reference_beds=regulatory_reference_beds,
            output_directory=output_dir,
            file_prefix=data_name)
        intersection_bed_dict[data_name] = sv_regulation_intersection_bed
        regulation_sv_annot_file = annot_func.annotate_sv_on_regulatory_element(
            sv_regulation_intersection_bed=sv_regulation_intersection_bed,
            regulatory_element_types=regulatory_element_types,
            output_directory=output_dir, file_prefix=data_name)
        regulation_sv_annot_dict[data_name] = regulation_sv_annot_file
    regulation_sv_annot_merge_sets_file = os.path.join(output_dir, 'Public_SV_sets.Merge.RE_annot.tsv')
    logging.info('Annote Regulation elements by AF...')
    annot_func.merge_filter_sv_sets_in_regulatory(input_files_dict=regulation_sv_annot_dict,
                                                  output_file=regulation_sv_annot_merge_sets_file,
                                                  overlap_threshold=overlap_threshold,
                                                  distance_threshold=distance_threshold)
    return regulation_sv_annot_merge_sets_file, regulation_sv_annot_dict


# %%
def run_annotation_only_public_sv(annotation_tuple, output_dir, input_sv_bed_paths, selected_datasets,
                                  overlap_threshold, distance_threshold):
    """
    Run annotation using only public structural variant datasets.

    Args:
        annotation_tuple (tuple): tuple of annotation data.
        output_dir (str): Directory to store output files.
        input_sv_bed_paths (dict): Paths to input SV bed files.
        selected_datasets (list): List of selected public SV datasets.
        overlap_threshold (float): Threshold for overlap in annotation.
        distance_threshold (int): Distance threshold for annotation.

    Returns:
        tuple: Paths to merged gene and regulatory element SV annotation files.
    """
    gene_annotation_dict, gene_reference_bed, regulatory_reference_beds, regulatory_element_types = annotation_tuple
    logging.info(f'The types of regulatory elements are {regulatory_element_types}')
    split_gene_sv_annot_merge_sets_file, split_gene_sv_annot_dict, whole_gene_sv_annot_dict = preprocess_public_sv_annot_on_gene(
        input_sv_bed_paths=input_sv_bed_paths,
        selected_datasets=selected_datasets,
        gene_reference_bed=gene_reference_bed,
        gene_annotation_dict=gene_annotation_dict,
        overlap_threshold=overlap_threshold,
        distance_threshold=distance_threshold,
        output_dir=output_dir)

    regulation_sv_annot_merge_sets_file, regulation_sv_annot_dict = preprocess_public_sv_annot_on_regulation(
        input_sv_bed_paths=input_sv_bed_paths,
        selected_datasets=selected_datasets,
        regulatory_element_types=regulatory_element_types,
        regulatory_reference_beds=regulatory_reference_beds,
        overlap_threshold=overlap_threshold,
        distance_threshold=distance_threshold,
        output_dir=output_dir)
    return split_gene_sv_annot_merge_sets_file, split_gene_sv_annot_dict, regulation_sv_annot_merge_sets_file, regulation_sv_annot_dict
