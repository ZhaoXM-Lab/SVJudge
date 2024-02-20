import logging
import os
import subprocess
import sys
from typing import Optional

import numpy as np
import pandas as pd
from communities.algorithms import louvain_method
from tqdm import tqdm

import src.utils as utils


# %%
def sv_gene_intersection(sv_bed: str, gene_annot_bed: str, output_dir: str, prefix: str) -> Optional[str]:
    """
    Annotate structural variants (SV) with gene using bedtools.

    :param sv_bed: Path to the SV BED file that needs annotation.
    :param gene_annot_bed: Path to the processed gene annotation BED file.
    :param output_dir: Directory for the output file.
    :param prefix: Prefix for the output file name.
    :return: Path to the annotated output file.
    """
    # Setup logging
    logging.basicConfig(filename='annotation_log.log', level=logging.INFO,
                        format='%(asctime)s:%(levelname)s:%(message)s')

    # Generate the output file name
    file_name = f'{prefix}.gene_annot.bedtools.temp.bed'
    output_file = os.path.join(output_dir, file_name)

    # Check if the input files exist
    if not os.path.exists(sv_bed):
        logging.error("Input BED files: %s not found." % sv_bed)
        sys.exit()
    if not os.path.exists(gene_annot_bed):
        logging.error("Annotation BED files: %s not found." % gene_annot_bed)
        sys.exit()
    # Run bedtools command to intersect the SV BED file with the gene annotation BED file
    bedtools_command = ['bedtools', 'intersect', '-a', sv_bed, '-b', gene_annot_bed, '-r', '-wo']
    try:
        with open(output_file, 'w') as out_file:
            subprocess.run(bedtools_command, stdout=out_file, check=True)
            return output_file
    except subprocess.CalledProcessError as e:
        logging.error(f"Error in running bedtools command: {e}")
        sys.exit()


# %%
def sv_regulatory_elements_intersection(sv_bed_path, regulatory_reference_beds, output_directory, file_prefix):
    """
    Annotate structural variants (SV) with regulatory elements (RE) information using bedtools.

    :param sv_bed_path: str, Path to the SV BED file that needs annotation.
    :param regulatory_reference_beds: List[str], Paths to the processed regulatory elements annotation BED files.
    :param output_directory: str, Directory for the output file.
    :param file_prefix: str, Prefix for the output file name.
    :return: str, Path to the annotated output file.
    """
    # Initialize logging

    logging.basicConfig(filename='annotation_log.log', level=logging.INFO,
                        format='%(asctime)s:%(levelname)s:%(message)s')

    # Check if input files exist
    if not os.path.exists(sv_bed_path):
        logging.error(f"SV BED file not found at {sv_bed_path}")
        return None

    for bed_path in regulatory_reference_beds:
        if not os.path.exists(bed_path):
            logging.error(f"Regulatory BED file not found at {bed_path}")
            return None

    # Create the output file name
    output_file_name = f'{file_prefix}.RE_annot.bedtools.temp.bed'
    output_file_path = os.path.join(output_directory, output_file_name)

    # Construct the bedtools command
    bedtools_command = [
                           'bedtools', 'intersect', '-a', sv_bed_path, '-b'
                       ] + regulatory_reference_beds + ['-r', '-wo']

    # Execute the bedtools command and handle potential errors
    try:
        with open(output_file_path, 'w') as output_file:
            subprocess.run(bedtools_command, stdout=output_file, check=True)
    except subprocess.CalledProcessError as e:
        logging.error(f"bedtools command failed: {e}")
        return None

    return output_file_path


# %%
# 前后断点所在位置
def annotate_sv_breakpoints_exon_intron(breakend_left, breakend_right, element_dict, gene_strand):
    """
    Annotate SV breakpoints location from the perspective of exons and introns.

    Args:
        breakend_left (int): Left breakpoint of the SV.
        breakend_right (int): Right breakpoint of the SV.
        element_dict (dict): Dictionary of elements (exons and introns) with element name as key and tuple of start and end as value.
        gene_strand (str): Strand of the gene, can be '+' or '-'.

    Returns:
        Tuple: containing locations of the left and right breakpoints.
    """
    # Initialize location variables
    location_left = ''
    location_right = ''

    # Adjust right breakpoint if it is the same as the left
    if breakend_left == breakend_right:
        breakend_right = breakend_left + 1

    # Annotate the breakpoint locations
    for element_name, (start, end) in element_dict.items():
        if start <= breakend_left < end:
            location_left = element_name
        if start < breakend_right <= end:
            location_right = element_name
        if location_left and location_right:
            break

    # Adjust the order of locations based on gene strand
    return (location_left, location_right) if gene_strand == '+' else (location_right, location_left)


# %%
def annotate_sv_breakpoints_cds_utr(breakend_left, breakend_right, cds, tx, gene_strand):
    """
    Annotate SV breakpoints location from the perspective of CDS and UTR.

    Args:
    - breakend_left: Left breakpoint of the SV.
    - breakend_right: Right breakpoint of the SV.
    - cds: CDS range.
    - tx: Transcript range.
    - gene_strand: Strand of the gene.

    Returns:
    Tuple containing locations of the left and right breakpoints.
    """
    # Initialize variables
    location_left = ''
    location_right = ''
    utr5 = '5\'UTR'
    utr3 = '3\'UTR'

    # Handle empty CDS case
    if abs(cds[1] - cds[0]) == 0:
        return 'UTR', 'UTR'

    # Determine UTRs based on gene strand
    utr1, utr2 = (utr5, utr3) if gene_strand == '+' else (utr3, utr5)

    # Annotate the breakpoint locations
    if tx[0] <= breakend_left < cds[0]:
        location_left = utr1
    elif cds[0] <= breakend_left < cds[1]:
        location_left = 'CDS'
    elif cds[1] <= breakend_left < tx[1]:
        location_left = utr2

    if tx[0] < breakend_right <= cds[0]:
        location_right = utr1
    elif cds[0] < breakend_right <= cds[1]:
        location_right = 'CDS'
    elif cds[1] < breakend_right <= tx[1]:
        location_right = utr2

    return location_left, location_right


# %%
def determine_impacted_elements(start_loc, end_loc):
    """
    Determine all exons and introns impacted by a structural variant (SV).

    Args:
    start_loc (str): The location impacted by the left breakpoint, e.g., 'exon1'.
    end_loc (str): The location impacted by the right breakpoint, e.g., 'exon2'.

    Returns:
    A string representing all impacted locations.
    """
    # Extract the element number from the start location
    element_number = int(''.join(filter(str.isdigit, start_loc)))

    impacted_elements = [start_loc]

    # Generate the list of impacted elements between start_loc and end_loc
    while impacted_elements[-1] != end_loc:
        if 'exon' in impacted_elements[-1]:
            impacted_elements.append(f'intron{element_number}')
        else:
            element_number += 1
            impacted_elements.append(f'exon{element_number}')

    return ';'.join(impacted_elements)


# %%
def determine_impacted_regions(start_region, end_region):
    """
    Determine the impacted coding and non-coding regions by a structural variant (SV).

    Args:
    - start_region: The location impacted by the left breakpoint, e.g., 'CDS'.
    - end_region: The location impacted by the right breakpoint, e.g., 'UTR'.

    Returns:
    A string representing the impacted coding and non-coding regions.
    """
    # Return the start_region if it is equal to the end_region
    if start_region == end_region:
        return start_region
    else:
        impacted_regions = set()

        # Add regions based on the start and end locations
        if '5\'UTR' in {start_region, end_region}:
            impacted_regions.update(['5\'UTR', 'CDS'])
        if '3\'UTR' in {start_region, end_region}:
            impacted_regions.add('3\'UTR')
        if 'CDS' in {start_region, end_region}:
            impacted_regions.add('CDS')

        return ';'.join(sorted(impacted_regions))


# %%
def is_svtype_match(input_type, target_type):
    """
    Determine if the SVTYPE matches the target type.

    Args:
    input_type (str): The SVTYPE to match.
    target_type (str): The target SVTYPE, which can be one of ['DEL', 'DUP', 'INS', 'INV', 'TRA'].

    Returns:
    bool: True if the types match.
    """
    # Check if the target type is 'INS'
    if target_type == 'INS':
        return input_type in ['ALU', 'LINE1', 'SVA', 'INS']  # Return True if input_type is in the specified list
    # Check if the target type is 'DEL'
    elif target_type == 'DEL':
        return 'DEL' in input_type or 'CNV' in input_type  # Return True if 'DEL' or 'CNV' is in input_type
    # Check if the target type is 'DUP'
    elif target_type == 'DUP':
        return input_type == 'DUP' or 'CNV' in input_type  # Return True if input_type is 'DUP' or 'CNV' is in input_type
    # Check if the target type is 'TRA'
    elif target_type == 'TRA':
        return input_type in ['TRA', 'BND', 'CTX']  # Return True if input_type is in the specified list
    # If target type is not any of the specified types, check if it matches input type
    else:
        return target_type == input_type  # Return True if target_type matches input_type


# %%
def standardize_sv_type(input_type):
    """
    Standardize the structural variant type to a uniform format.

    Args:
    input_type (str): The input structural variant type.

    Returns:
    str: One of the six standardized structural variant types:
    ['DEL', 'DUP', 'INS', 'INV', 'TRA', 'CNV'].
    """
    # Standardize the input_type to a uniform format
    if input_type in ['ALU', 'LINE1', 'SVA', 'INS']:
        return 'INS'  # Insertion
    elif 'DEL' in input_type:
        return 'DEL'  # Deletion
    elif input_type == 'DUP':
        return 'DUP'  # Duplication
    elif 'CNV' in input_type:
        return 'CNV'  # Copy number variation
    elif input_type in ['TRA', 'BND', 'CTX']:
        return 'TRA'  # Translocation
    else:
        return input_type


# %%

def calculate_sv_impact_on_cds_length(sv_breakend_start, sv_breakend_end, exon_introns_dict,
                                      impacted_exons_str, cds_range, sv_data):
    """
    Calculate the length of the impact of a structural variant (SV) on the coding sequence (CDS).

    Args:
    - sv_breakend_start: Start position of the SV breakend.
    - sv_breakend_end: End position of the SV breakend.
    - exon_introns_dict: Dictionary of exons and introns.
    - impacted_exons_str: String of impacted exons separated by semicolons.
    - cds_range: Range of the coding sequence (CDS).
    - sv_data: Data array containing information about the SV.

    Returns:
    Length of SV impact on the CDS.
    """
    # Initialize impact length
    impact_length = 0

    # Get the SV type
    sv_type = sv_data[4]

    # If SV type is insertion and the SV impacts the CDS
    if is_svtype_match(sv_type, 'INS'):
        if 'exon' in impacted_exons_str and cds_range[0] <= sv_breakend_start < cds_range[1]:
            # If the impact length is unknown, return 'unknown'
            if sv_data[5] == 'unknown':
                return 'unknown'
            # Otherwise, calculate impact length
            impact_length = abs(int(sv_data[5]))
        return impact_length

    # If the CDS range has zero length, return 0
    if abs(cds_range[0] - cds_range[1]) == 0:
        return 0

    # Iterate through impacted exons
    for exon in impacted_exons_str.split(';'):
        # If the exon is impacted
        if 'exon' in exon:
            # Get start and end positions of the exon
            exon_start, exon_end = exon_introns_dict[exon]
            # Find the overlap between the SV breakend and CDS
            breakend_max = max(exon_start, sv_breakend_start)
            breakend_min = min(exon_end, sv_breakend_end)

            # If the breakend is outside of the CDS range, continue to the next exon
            if breakend_max >= cds_range[1] or breakend_min <= cds_range[0]:
                continue

            # Calculate the impact length
            impact_length += min(breakend_min, cds_range[1]) - max(breakend_max, cds_range[0])

    return impact_length


# %%
def determine_sv_gene_relationship(sv_start, sv_end, breakpoints_exon_intron, breakpoints_cds_utr, cds_range,
                                   gene_info):
    """
    Determine the relationship of a structural variant (SV) to a gene.

    Args:
        sv_start (int): Start position of the SV.
        sv_end (int): End position of the SV.
        breakpoints_exon_intron (tuple): Tuple containing the left and right breakpoint annotations from the exon/intron perspective.
        breakpoints_cds_utr (tuple): Tuple containing the left and right breakpoint annotations from the CDS/UTR perspective.
        cds_range (tuple): Start and end positions of the gene's coding sequence (CDS).
        gene_info (dict): A dictionary containing gene data.

    Returns:
        str: Relationship of the SV to the gene or None in case of an error.
    """
    sv_exon_intron_left, sv_exon_intron_right = breakpoints_exon_intron
    sv_cds_utr_left, sv_cds_utr_right = breakpoints_cds_utr

    if sv_start <= gene_info['txStart'] <= gene_info['txEnd'] <= sv_end:
        if cds_range[0] == cds_range[1]:
            return 'No CDS in the gene'
        else:
            return 'SV contains all Tx'
    elif sv_exon_intron_left == sv_exon_intron_right and 'intron' in sv_exon_intron_left:
        return 'SV in Intronic'
    elif cds_range[0] == cds_range[1]:
        return 'No CDS in the gene'
    elif sv_cds_utr_left == sv_cds_utr_right and 'UTR' in sv_cds_utr_left:
        return 'SV in UTR'
    elif sv_cds_utr_left == sv_cds_utr_right and 'CDS' in sv_cds_utr_left:
        return 'SV in Coding or Splice'
    elif sv_cds_utr_left != sv_cds_utr_right:
        return 'SV Partial overlap Tx'
    else:
        return None


def annotate_sv_on_whole_gene(sv_gene_intersection_bed, gene_annot_dict, output_dir, file_prefix):
    """
    Generate a comprehensive annotation file for SV impact on whole genes.
    :param sv_gene_intersection_bed: Path to the BED file with overlaps of SV and gene elements annotations from bedtools.
    :param gene_annot_dict: Dictionary of canonical gene annotations.
    :param output_dir: Directory for output files.
    :param file_prefix: Prefix for the output file name.
    :return: Path to the output file containing gene annotation information.
    """
    # Open the input BED file
    with open(sv_gene_intersection_bed, 'r') as input_bed_file:
        output_file_path = os.path.join(output_dir, f'{file_prefix}.gene_annot.full.tsv')
        with open(output_file_path, 'w') as output_file:
            # Write header line to the output file
            header_columns = ['#CHR', 'START', 'END', 'ID', 'SVTYPE', 'SVLEN', 'AF',
                              'txID', 'geneID', 'breakend1', 'breakend2', 'relationship', 'OverlapElements',
                              'OverlapRegions', 'Overlap_CDS_Length', 'Frameshit']
            output_file.write('\t'.join(header_columns) + '\n')
            for line_index, line in enumerate(input_bed_file):
                line_data = line.strip().split('\t')
                chr_name, sv_start, sv_end, tx_id = line_data[0], int(line_data[1]), int(line_data[2]), line_data[10]

                gene_info = gene_annot_dict[str(chr_name)].loc[tx_id]
                sv_breakend_start, sv_breakend_end = max(sv_start, int(gene_info['txStart'])), min(sv_end, int(
                    gene_info['txEnd']))

                # Process elements and regions impacted by SV
                elements_dict = {}  # Processed elements dictionary
                gene_strand = gene_info['Strand']
                exon_starts = gene_info['exonStarts'].split(',')[:-1]
                exon_ends = gene_info['exonEnds'].split(',')[:-1]
                for index_k in range(int(gene_info['exonCount'])):
                    if gene_strand == '+':
                        exon_k = index_k + 1
                        intron_k = index_k + 1
                    else:
                        exon_k = int(gene_info['exonCount']) - index_k
                        intron_k = int(gene_info['exonCount']) - index_k - 1
                    elements_dict['exon' + str(exon_k)] = [int(exon_starts[index_k]), int(exon_ends[index_k])]
                    if 0 < intron_k < int(gene_info['exonCount']):
                        elements_dict['intron' + str(intron_k)] = [int(exon_ends[index_k]),
                                                                   int(exon_starts[index_k + 1])]
                breakpoint_exon_intron_left, breakpoint_exon_intron_right = annotate_sv_breakpoints_exon_intron(
                    sv_breakend_start, sv_breakend_end, elements_dict, gene_strand)
                if breakpoint_exon_intron_left == breakpoint_exon_intron_right == '':
                    continue
                impact_element_str = determine_impacted_elements(breakpoint_exon_intron_left,
                                                                 breakpoint_exon_intron_right)
                cds_range = (int(gene_info['cdsStart']), int(gene_info['cdsEnd']))
                tx_region = (int(gene_info['txStart']), int(gene_info['txEnd']))
                breakpoint_cds_utr_left, breakpoint_cds_utr_right = annotate_sv_breakpoints_cds_utr(
                    sv_breakend_start, sv_breakend_end, cds_range, tx_region, gene_strand)
                impact_regions_str = determine_impacted_regions(breakpoint_cds_utr_left,
                                                                breakpoint_cds_utr_right)
                overlap_cds_length = calculate_sv_impact_on_cds_length(sv_breakend_start, sv_breakend_end,
                                                                       elements_dict, impact_element_str, cds_range,
                                                                       line_data)
                relationship = determine_sv_gene_relationship(sv_start, sv_end,
                                                              breakpoints_exon_intron=(breakpoint_exon_intron_left,
                                                                                       breakpoint_exon_intron_right),
                                                              breakpoints_cds_utr=(
                                                                  breakpoint_cds_utr_left, breakpoint_cds_utr_right)
                                                              , cds_range=cds_range, gene_info=gene_info)
                if relationship is None:
                    logging.error("Location con not identified")

                w_lines = line_data[0:7]
                w_lines.extend([line_data[10], gene_info['geneID']
                                   , str(sv_breakend_start)
                                   , str(sv_breakend_end)
                                   , str(relationship)
                                   , str(impact_element_str)
                                   , str(impact_regions_str)
                                   , str(overlap_cds_length)])
                if overlap_cds_length == 'unknown':
                    w_lines.append('')
                elif overlap_cds_length < 0:
                    logging.error('error! length')
                    break
                elif overlap_cds_length % 3 != 0:
                    w_lines.append('Yes')
                else:
                    w_lines.append('No')

                w_line = '\t'.join(w_lines)
                output_file.write('%s\n' % w_line)
            return output_file_path


# %%


def annotate_sv_on_regulatory_element(sv_regulation_intersection_bed, regulatory_element_types, output_directory,
                                      file_prefix):
    """
    Generate a comprehensive annotation file for the impact of structural variants on regulatory elements.

    Args:
    - sv_regulation_intersection_bed: Path to the BED file with overlaps of SV and regulatory elements annotations.
    - regulatory_element_types: List of regulatory elements keywords like promoter, enhancer, TAD boundary.
    - output_directory: Directory for output files.
    - file_prefix: Prefix for the output file name.

    Returns:
    - Path to the annotated output file.
    """

    # Open the input BED file
    with open(sv_regulation_intersection_bed, 'r') as input_bed_file:
        # Define the output file name and path
        output_file_name = f'{file_prefix}.RE_annot.tsv'
        output_file_path = os.path.join(output_directory, output_file_name)

        # Write the header columns to the output file
        with open(output_file_path, 'w') as output_file:
            header_columns = ['#CHR', 'START', 'END', 'ID', 'SVTYPE', 'SVLEN', 'AF',
                              'RE_ID', 'RE_type', 'geneID', 'breakend1', 'breakend2', 'Overlap_rate', 'relationship']
            output_file.write('\t'.join(header_columns) + '\n')

            # Iterate through each line in the input BED file
            for line in tqdm(input_bed_file):
                columns = line.strip().split('\t')
                gene_id = columns[11]
                sv_start, sv_end = int(columns[1]), int(columns[2])
                re_start, re_end = int(columns[9]), int(columns[10])
                re_id = columns[12]

                breakend_start, breakend_end = max(sv_start, re_start), min(sv_end, re_end)
                overlap_rate = abs(breakend_start - breakend_end) / (re_end - re_start)
                re_type = regulatory_element_types[int(columns[7]) - 1]
                relationship = f'SV overlap {re_type}'

                new_line = columns[:7] + [re_id, re_type, gene_id, str(breakend_start), str(breakend_end),
                                          str(overlap_rate), relationship]
                output_file.write('\t'.join(new_line) + '\n')
    return output_file_path


# %%
def split_gene_elements(whole_gene_sv_annot_file, gene_annot_dict):
    """
    Split SV annotations on a gene into specific regions (exon, intron, CDS, UTR).
    :param whole_gene_sv_annot_file: File with SV annotations on the whole gene.
    :param gene_annot_dict: Dictionary of gene annotations.
    :return: Path to the output file with split annotations.
    """
    df_full = pd.read_csv(whole_gene_sv_annot_file, sep='\t')
    out_split_file = whole_gene_sv_annot_file.replace('.gene_annot.full.tsv', '.gene_annot.split.tsv')
    out_h = open(out_split_file, 'w')
    w_s_col = ['#CHR', 'START', 'END', 'ID', 'SVTYPE', 'SVLEN', 'txID', 'geneID', 'Element', 'Region',
               'breakend1', 'breakend2', 'Overlap_rate', 'AF']
    out_h.write('%s\n' % '\t'.join(w_s_col))
    for index, row in tqdm(df_full.iterrows()):
        tx_id = row['txID']
        sv_start = int(row['START'])
        sv_end = int(row['END'])
        overlap_elements = row['OverlapElements'].split(';')
        overlap_regions = row['OverlapRegions'].split(';')
        chr_name = row['#CHR']
        gene_info = gene_annot_dict[str(chr_name)].loc[tx_id]
        exons_start_list = gene_info['exonStarts'].split(',')
        exons_end_list = gene_info['exonEnds'].split(',')
        cds_start = int(gene_info['cdsStart'])
        cds_end = int(gene_info['cdsEnd'])
        strand = gene_info['Strand']
        exon_count = int(gene_info['exonCount'])
        sv_af = row['AF']
        for exon_name in overlap_elements:
            if strand == '+':
                exon_index = int(''.join(filter(str.isdigit, exon_name))) - 1
                if 'exon' in exon_name:
                    exon_breakend_start = int(exons_start_list[exon_index])
                    exon_breakend_end = int(exons_end_list[exon_index])
                else:
                    exon_breakend_start = int(exons_end_list[exon_index])
                    exon_breakend_end = int(exons_start_list[exon_index + 1])
            else:
                exon_index = exon_count - int(''.join(filter(str.isdigit, exon_name)))
                if 'exon' in exon_name:
                    exon_breakend_start = int(exons_start_list[exon_index])
                    exon_breakend_end = int(exons_end_list[exon_index])
                else:
                    exon_breakend_start = int(exons_end_list[exon_index - 1])
                    exon_breakend_end = int(exons_start_list[exon_index])
            breakpoint_start = max(exon_breakend_start, sv_start)
            breakpoint_end = min(exon_breakend_end, sv_end)
            if len(overlap_regions) == 1:
                new_breakpoint_start = breakpoint_start
                new_breakpoint_end = breakpoint_end
                region = overlap_regions[0]
                exon_overlap_rate = abs(new_breakpoint_start - new_breakpoint_end) / abs(
                    exon_breakend_start - exon_breakend_end)
                w_lines = [chr_name, row['START'], row['END'], row['ID'], row['SVTYPE'], row['SVLEN'], row['txID'],
                           row['geneID'], exon_name, region, new_breakpoint_start, new_breakpoint_end,
                           exon_overlap_rate, sv_af]
                w_lines = [str(x) for x in w_lines]
                out_h.write('%s\n' % '\t'.join(w_lines))
            else:
                if strand == '+':
                    utr1 = '5\'UTR'
                    utr2 = '3\'UTR'
                else:
                    utr1 = '3\'UTR'
                    utr2 = '5\'UTR'
                if cds_start < breakpoint_end and breakpoint_start < cds_end:
                    new_breakpoint_start = max(breakpoint_start, cds_start)
                    new_breakpoint_end = min(breakpoint_end, cds_end)
                    region = 'CDS'
                    exon_overlap_rate = abs(new_breakpoint_start - new_breakpoint_end) / abs(
                        max(exon_breakend_start, cds_start) - min(exon_breakend_end, cds_end))
                    w_lines = [chr_name, row['START'], row['END'], row['ID'], row['SVTYPE'], row['SVLEN'], row['txID'],
                               row['geneID'], exon_name, region, new_breakpoint_start, new_breakpoint_end,
                               exon_overlap_rate, sv_af]
                    w_lines = [str(x) for x in w_lines]
                    out_h.write('%s\n' % '\t'.join(w_lines))
                    if breakpoint_start < cds_start:
                        new_breakpoint_start = breakpoint_start
                        new_breakpoint_end = cds_start
                        region = utr1
                        exon_overlap_rate = abs(new_breakpoint_start - new_breakpoint_end) / abs(
                            exon_breakend_start - cds_start)
                        w_lines = [chr_name, row['START'], row['END'], row['ID'], row['SVTYPE'], row['SVLEN'],
                                   row['txID'],
                                   row['geneID'], exon_name, region, new_breakpoint_start, new_breakpoint_end,
                                   exon_overlap_rate, sv_af]
                        w_lines = [str(x) for x in w_lines]
                        out_h.write('%s\n' % '\t'.join(w_lines))
                    if cds_end < breakpoint_end:
                        new_breakpoint_start = cds_end
                        new_breakpoint_end = breakpoint_end
                        region = utr2
                        exon_overlap_rate = abs(new_breakpoint_start - new_breakpoint_end) / abs(
                            exon_breakend_end - cds_end)
                        w_lines = [chr_name, row['START'], row['END'], row['ID'], row['SVTYPE'], row['SVLEN'],
                                   row['txID'],
                                   row['geneID'], exon_name, region, new_breakpoint_start, new_breakpoint_end,
                                   exon_overlap_rate, sv_af]
                        w_lines = [str(x) for x in w_lines]
                        out_h.write('%s\n' % '\t'.join(w_lines))
                else:
                    new_breakpoint_start = breakpoint_start
                    new_breakpoint_end = breakpoint_end
                    exon_overlap_rate = abs(new_breakpoint_start - new_breakpoint_end) / abs(
                        exon_breakend_start - exon_breakend_end)
                    if cds_end <= breakpoint_start:
                        region = utr2
                    else:
                        region = utr1
                    w_lines = [chr_name, row['START'], row['END'], row['ID'], row['SVTYPE'], row['SVLEN'], row['txID'],
                               row['geneID'], exon_name, region, new_breakpoint_start, new_breakpoint_end,
                               exon_overlap_rate, sv_af]
                    w_lines = [str(x) for x in w_lines]
                    out_h.write('%s\n' % '\t'.join(w_lines))
    out_h.close()
    return out_split_file


# %%
def _calculate_overlap_similarity(start1, end1, start2, end2, overlap_threshold):
    """
     Calculate the overlap similarity between two ranges based on their start and end points and an overlap threshold.
     Parameters:
         start1 (int): The start point of the first range.
         end1 (int): The end point of the first range.
         start2 (int): The start point of the second range.
         end2 (int): The end point of the second range.
         overlap_threshold (float): The threshold for considering overlap similarity.
     Returns:
         float: The calculated overlap similarity between the two ranges.
     """
    overlap = max(0, min(end1, end2) - max(start1, start2))
    s1 = overlap / (end1 - start1)
    s2 = overlap / (end2 - start2)
    if s1 >= overlap_threshold and s2 >= overlap_threshold or max(s1, s2) >= overlap_threshold and abs(s1 - s2) <= 0.1:
        return s1 + s2
    return 0


def calculate_sv_similarity(start1, end1, start2, end2, sv_type=None, sv_len1=None, sv_len2=None,
                            distance_threshold=100, overlap_threshold=0.5):
    """
    Calculate the similarity score between two structural variants (SVs).

    Prioritizes mutual coverage. For translocations (TRA) and SVs without length, considers the distance between breakpoints.

    Args:
        start1 (int): Left breakpoint of the first SV.
        end1 (int): Right breakpoint of the first SV.
        start2 (int): Left breakpoint of the second SV.
        end2 (int): Right breakpoint of the second SV.
        sv_type (str, optional): Type of the structural variant. Defaults to None.
        sv_len1 (int, optional): Length of the first SV. Defaults to None.
        sv_len2 (int, optional): Length of the second SV. Defaults to None.
        distance_threshold (int, optional): Maximum threshold for breakpoint distance, used for TRA or unannotated INS. Defaults to 100.
        overlap_threshold (float, optional): Minimum threshold for SV overlap. Defaults to 0.5.

    Returns:
        float: The calculated similarity score between the two SVs.
    """
    # Handle potential non-integer inputs for sv_len1 and sv_len2
    try:
        sv_len1 = abs(int(sv_len1)) if sv_len1 is not None else 0
        sv_len2 = abs(int(sv_len2)) if sv_len2 is not None else 0
    except ValueError:
        sv_len1, sv_len2 = 0, 0

    # Convert inputs to integers
    start1, end1, start2, end2 = map(int, [start1, end1, start2, end2])

    # Check if sv_type is not INS, BND, or TRA, then calculate overlap similarity
    if sv_type not in ['INS', 'BND', 'TRA']:
        return _calculate_overlap_similarity(start1, end1, start2, end2, overlap_threshold)

    # Check if sv_len1 and sv_len2 are greater than 1, then calculate overlap similarity
    if sv_len1 > 1 and sv_len2 > 1:
        return _calculate_overlap_similarity(start1, start1 + sv_len1, start2, start2 + sv_len2, overlap_threshold)

    # Calculate similarity score based on distance between breakpoints
    return max(distance_threshold - abs(start2 - start1), 0) / distance_threshold * 2


# %%
def filter_sv_by_overlap(df_data, sv_type, overlap_threshold=None, distance_threshold=None):
    """
    Filter structural variants (SVs) based on mutual overlap.

    Args:
        df_data (DataFrame): Input data containing SV information.
        sv_type (str): Type of the structural variant.
        overlap_threshold (float, optional): Minimum threshold for considering overlap. Defaults to None.
        distance_threshold (float, optional): Maximum threshold for breakpoint distance. Defaults to None.

    Returns:
        DataFrame: Filtered DataFrame after applying the overlap criteria.
    """
    # Make a copy of the input DataFrame
    df_filtered = df_data.copy()
    df_filtered.reset_index(drop=True, inplace=True)

    # Get the number of records in the DataFrame
    num_records = len(df_filtered)

    # Create an empty matrix to store overlap information
    overlap_matrix = np.zeros((num_records, num_records))

    # Calculate similarity matrix
    for i in range(num_records):
        for j in range(num_records):
            if i != j:
                # Calculate overlap between SVs
                overlap_matrix[i, j] = calculate_sv_similarity(
                    start1=df_filtered.at[i, 'breakend1'],
                    end1=df_filtered.at[i, 'breakend2'],
                    start2=df_filtered.at[j, 'breakend1'],
                    end2=df_filtered.at[j, 'breakend2'],
                    sv_type=sv_type,
                    sv_len1=df_filtered.at[i, 'SVLEN'],
                    sv_len2=df_filtered.at[j, 'SVLEN'],
                    distance_threshold=distance_threshold,
                    overlap_threshold=overlap_threshold
                )

    # Apply Louvain method to identify clusters
    clusters, _ = louvain_method(overlap_matrix)

    # Iterate over clusters and remove non-maximum allele frequency records
    for cluster in clusters:
        if len(cluster) > 1:
            df_cluster = df_filtered.loc[list(cluster)]
            index_of_max_af = df_cluster['AF'].idxmax()
            indices_to_remove = list(set(cluster) - {index_of_max_af})
            df_filtered.drop(indices_to_remove, inplace=True)

    # Reset the index of the filtered DataFrame
    df_filtered.reset_index(drop=True, inplace=True)

    return df_filtered


# %%
def filter_sv_by_containment(df_data, sv_type):
    """
    Filter structural variants (SVs) based on containment relationships.

    Args:
        df_data (DataFrame): Input DataFrame with SV data.
        sv_type (str): Type of the structural variant.

    Returns:
        DataFrame: Filtered DataFrame after applying containment criteria.
    """
    # Make a copy of the input DataFrame
    df_filtered = df_data.copy()

    # Reset the index of the DataFrame
    df_filtered.reset_index(drop=True, inplace=True)

    # Sort the DataFrame based on 'Overlap_rate' and 'AF' columns
    df_filtered.sort_values(by=['Overlap_rate', 'AF'], ascending=(False, False), ignore_index=True, inplace=True)

    # Get the number of records in the DataFrame
    num_records = len(df_filtered)

    # Get the start positions and SV lengths as lists
    starts = df_filtered['breakend1'].tolist()
    sv_lengths = df_filtered['SVLEN'].tolist()

    # Calculate the end positions based on the SV type
    if sv_type == 'INS' and 'unknown' not in sv_lengths:
        ends = [start + int(length) for start, length in zip(starts, sv_lengths)]
    else:
        ends = df_filtered['breakend2'].tolist()

    # Get the AF values as a list
    af_values = df_filtered['AF'].tolist()

    # Initialize a list to store indices of records to remove
    indices_to_remove = []

    # Iterate through the records to find contained SVs
    for i in range(num_records):
        for j in range(i + 1, num_records):
            if starts[i] <= starts[j] and ends[j] <= ends[i]:
                if af_values[i] >= af_values[j]:
                    indices_to_remove.append(j)

    # Drop the records based on the indices to remove
    df_filtered.drop(indices_to_remove, inplace=True)

    # Reset the index of the filtered DataFrame
    df_filtered.reset_index(drop=True, inplace=True)

    return df_filtered


# %%
def dataframe_gene_split_filter(df_input, out_file, overlap_threshold, distance_threshold):
    """
    Merge and filter public data structural variation (SV) sets to obtain Allele Frequency (AF) for gene elements.

    Args:
        df_input (DataFrame): Input DataFrame containing SV data.
        out_file (str): Path to the output file where the filtered SV data will be saved.
        overlap_threshold (float): Threshold for mutual overlap of SVs.
        distance_threshold (float): Threshold for distance between SV breakpoints.

    Returns:
        str: Path to the output file containing the filtered SV data.
    """
    # Standardize SVTYPE column
    df_input['U_SVTYPE'] = df_input['SVTYPE'].apply(standardize_sv_type)

    # Group by gene elements and region
    df_grouped = df_input.groupby(by=['txID', 'U_SVTYPE', 'Element', 'Region'])

    merged_values = []

    # Iterate through each group
    for _, group_df in tqdm(df_grouped):
        if len(group_df) == 1:
            # If only one SV in the group, add it to the merged values
            merged_values.extend(group_df.values)
        else:
            # If multiple SVs in the group, filter and merge them
            filtered_df = filter_sv_by_overlap(group_df, group_df['U_SVTYPE'].iloc[0], overlap_threshold,
                                               distance_threshold)
            if len(filtered_df) > 1:
                # Further filter if still multiple SVs
                filtered_df = filter_sv_by_containment(filtered_df, filtered_df['U_SVTYPE'].iloc[0])
            # Add filtered SVs to the merged values
            merged_values.extend(filtered_df.values)

    # Create a new DataFrame with the merged values and save it to the output file
    df_final = pd.DataFrame(merged_values, columns=df_input.columns)
    df_final.to_csv(out_file, sep='\t', index=False)

    return out_file


# %%
def dataframe_re_split_filter(df_input, out_file, overlap_threshold, distance_threshold):
    """
    Merge and filter public data SV sets to obtain Allele Frequency (AF) for regulatory elements.

    Args:
        df_input (DataFrame): Input DataFrame with SV data.
        out_file (str): Path to the output file.
        overlap_threshold (float): Threshold for mutual overlap.
        distance_threshold (float): Threshold for breakpoint distance.

    Returns:
        str: Path to the output file.
    """
    # Standardize SV type
    df_input['U_SVTYPE'] = df_input['SVTYPE'].apply(standardize_sv_type)

    # Group by regulatory element ID, type, and standardized SV type
    df_grouped = df_input.groupby(by=['RE_ID', 'RE_type', 'U_SVTYPE'])

    # Initialize list for merged values
    merged_values = []

    # Iterate through grouped DataFrame
    for _, group_df in tqdm(df_grouped):
        # If only one row in group, add it to merged values
        if len(group_df) == 1:
            merged_values.extend(group_df.values)
        # If multiple rows in group, filter and merge
        else:
            filtered_df = filter_sv_by_overlap(group_df, group_df['U_SVTYPE'].iloc[0], overlap_threshold,
                                               distance_threshold)
            # If multiple filtered rows, further filter by containment
            if len(filtered_df) > 1:
                filtered_df = filter_sv_by_containment(filtered_df, filtered_df['U_SVTYPE'].iloc[0])
            merged_values.extend(filtered_df.values)

    # Create final DataFrame from merged values
    df_final = pd.DataFrame(merged_values, columns=df_input.columns)

    # Write final DataFrame to output file
    df_final.to_csv(out_file, sep='\t', index=False)

    # Return path to the output file
    return out_file


# %%
def merge_filter_sv_sets_in_gene(input_sv_af_dict, output_file, overlap_threshold, distance_threshold):
    """
    Merge and filter public database SV sets to obtain Allele Frequencies (AF) for gene body.

    Args:
        input_sv_af_dict (dict): Dictionary of paths to public database SV sets containing AF (annotated for genes).
        output_file (str): Path for the output file.
        overlap_threshold (float): Threshold for mutual overlap.
        distance_threshold (float): Threshold for breakpoint distance.

    Returns:
        str: Path to the output file.
    """
    # Initialize an empty dataframe
    merged_df = pd.DataFrame()

    # Define the columns for the merged dataframe
    columns = ['#CHR', 'ID', 'SVTYPE', 'SVLEN', 'txID', 'geneID', 'Element', 'Region', 'breakend1',
               'breakend2', 'Overlap_rate', 'AF']

    # Merge the datasets and assign dataset name to each entry
    for dataset in input_sv_af_dict.values():
        temp_df = pd.read_csv(dataset, sep='\t')[columns]
        temp_df['Datasets'] = dataset  # Assigning dataset name to each entry
        merged_df = merged_df.append(temp_df, ignore_index=True)

    # Call the function to filter the merged dataframe and write to the output file
    dataframe_gene_split_filter(df_input=merged_df, out_file=output_file,
                                overlap_threshold=overlap_threshold, distance_threshold=distance_threshold)

    return output_file


# %%

def merge_filter_sv_sets_in_regulatory(input_files_dict, output_file, overlap_threshold, distance_threshold):
    """
    Merge and filter structural variation (SV) sets to obtain Allele Frequencies (AF) for regulatory elements.

    Args:
        input_files_dict (dict): Dictionary of paths to SV sets annotated for regulatory elements.
        output_file (str): Path for the output file.
        overlap_threshold (float): Threshold for mutual overlap.
        distance_threshold (float): Threshold for breakpoint distance.

    Returns:
        str: Path to the output file.
    """
    # Create an empty DataFrame
    merged_df = pd.DataFrame()

    # Define the columns for the merged DataFrame
    columns = ['#CHR', 'ID', 'SVTYPE', 'SVLEN', 'RE_ID', 'RE_type', 'geneID', 'breakend1', 'breakend2',
               'Overlap_rate', 'AF']

    # Iterate through each dataset in the input files dictionary
    for dataset in input_files_dict.values():
        # Read the dataset into a temporary DataFrame and select the specified columns
        temp_df = pd.read_csv(dataset, sep='\t')[columns]
        # Assign the dataset name to each entry
        temp_df['Datasets'] = dataset
        # Append the temporary DataFrame to the merged DataFrame
        merged_df = merged_df.append(temp_df, ignore_index=True)

    # Call the function to filter and split the merged DataFrame, and write the result to the output file
    dataframe_re_split_filter(df_input=merged_df, out_file=output_file,
                              overlap_threshold=overlap_threshold, distance_threshold=distance_threshold)

    # Return the path to the output file
    return output_file


# %%
def _get_union(range_a, range_b):
    """
    Calculate the union of two ranges.

    Args:
        range_a (list/tuple): First range as [lower_bound, upper_bound].
        range_b (list/tuple): Second range as [lower_bound, upper_bound].

    Returns:
        list: A list representing the union of the two ranges.
              Format: [status_code, lower_bound, upper_bound].
              Status code 1 means a single continuous range,
              and 2 means two distinct ranges.
    """
    lower_a, upper_a = range_a
    lower_b, upper_b = range_b

    if upper_a + 1 == lower_b:
        # Ranges are consecutive, merge them into a single continuous range.
        return [1, lower_a, upper_b]
    if upper_a < lower_b:
        # Ranges are distinct, return them as separate.
        return [2, range_a, range_b]
    if lower_a <= lower_b:
        # Ranges overlap, find the combined range.
        combined_upper = max(upper_a, upper_b)
        return [1, lower_a, combined_upper]


def get_nset_union(section_lists):
    """
    Calculate the union of multiple ranges.

    Args:
        section_lists (list of lists/tuples): List of ranges, each represented as [lower_bound, upper_bound].

    Returns:
        list: List of merged ranges after calculating the union.
    """
    # Sort the ranges based on their start points
    section_lists.sort(key=lambda range_item: range_item[0])

    i = 0
    while i < len(section_lists) - 1:
        current_range, next_range = section_lists[i], section_lists[i + 1]
        union_result = _get_union(current_range, next_range)

        if union_result[0] == 1:
            # Merge the ranges and replace the current range with the merged range
            section_lists[i] = [union_result[1], union_result[2]]
            del section_lists[i + 1]  # Remove the next range as it's merged
        elif union_result[0] == 2:
            # Move to the next range if no merge is possible
            i += 1

    return section_lists


# %%
def is_contained(bench_start, bench_end, compare_start, compare_end, bench_sv_len, compare_sv_len, sv_type):
    """
    Determine if one range is completely contained within another range.

    Args:
        bench_start (int): Start position of the benchmark range.
        bench_end (int): End position of the benchmark range.
        compare_start (int): Start position of the compared range.
        compare_end (int): End position of the compared range.
        bench_sv_len (int): Length of the benchmark structural variation.
        compare_sv_len (int): Length of the compared structural variation.
        sv_type (str): Type of structural variation (e.g., INS, BND, TRA).

    Returns:
        int: 1 if the compared range is completely contained within the benchmark range, else 0.
    """
    # Ensure sv_len values are positive integers
    if sv_type in ['INS', 'BND', 'TRA']:
        try:
            bench_sv_len = abs(int(bench_sv_len))
            compare_sv_len = abs(int(compare_sv_len))
        except ValueError:
            bench_sv_len = 0
            compare_sv_len = 0
        # Adjust end positions based on sv_len
        if bench_sv_len != 0 and compare_sv_len != 0:
            bench_end = int(bench_start) + bench_sv_len
            compare_end = int(compare_start) + compare_sv_len

    # Check if compare range is completely contained within benchmark range
    if int(bench_start) <= int(compare_start) <= int(compare_end) <= int(bench_end):
        return 1
    else:
        return 0


# %%
def compare_and_obtain_af(row, df_annot_slice, overlap_threshold, distance_threshold):
    """
    Compare a given structural variant (SV) against annotated data to determine the most relevant allele frequency (AF).

    Args:
        row (pd.Series): The row representing the SV to be compared.
        df_annot_slice (pd.DataFrame): The slice of annotated data for comparison.
        overlap_threshold (float): Threshold for overlap rate.
        distance_threshold (int): Threshold for distance in SV comparison.

    Returns:
        tuple: The most relevant AnnotSV_ID, Dataset name, and AF.
    """
    # If the slice of annotated data is empty, return unknown values
    if df_annot_slice.empty:
        return 'unknown', 'unknown', 0

    # Rule 1: Calculate similarity rate and find the highest similarity
    bk1, bk2, sv_len, sv_type = row['breakend1'], row['breakend2'], row['SVLEN'], row['SVTYPE']
    df_annot_slice['Similarity_rate'] = df_annot_slice.apply(
        lambda x: calculate_sv_similarity(
            start1=bk1, end1=bk2, start2=x['breakend1'], end2=x['breakend2'],
            sv_type=sv_type, sv_len1=sv_len, sv_len2=x['SVLEN'],
            distance_threshold=distance_threshold, overlap_threshold=overlap_threshold
        ), axis=1)

    # Filter rows with non-zero similarity rate
    df_rule1 = df_annot_slice[df_annot_slice['Similarity_rate'] != 0]
    if not df_rule1.empty:
        # Sort by similarity rate and allele frequency (AF), and return the top row
        df_rule1.sort_values(by=['Similarity_rate', 'AF'], ascending=[False, False], inplace=True)
        top_row = df_rule1.iloc[0]
        return top_row['ID'], top_row['Datasets'], top_row['AF']

    # Rule 2: Check for containment
    df_annot_slice['Contain_target_SV'] = df_annot_slice.apply(
        lambda x: is_contained(
            bench_start=x['breakend1'], bench_end=x['breakend2'],
            compare_start=bk1, compare_end=bk2, bench_sv_len=x['SVLEN'],
            compare_sv_len=sv_len, sv_type=sv_type
        ), axis=1)

    # Filter rows with non-zero containment
    df_rule2 = df_annot_slice[df_annot_slice['Contain_target_SV'] != 0]
    if not df_rule2.empty:
        # Sort by AF and return the top row
        df_rule2.sort_values(by='AF', ascending=False, inplace=True)
        top_row = df_rule2.iloc[0]
        return top_row['ID'], top_row['Datasets'], top_row['AF']

    # Rule 3: Check if the SV contains more than one annotated SV
    df_annot_slice['Be_Contained'] = df_annot_slice.apply(
        lambda x: is_contained(
            bench_start=bk1, bench_end=bk2, compare_start=x['breakend1'],
            compare_end=x['breakend2'], bench_sv_len=sv_len, compare_sv_len=x['SVLEN'],
            sv_type=sv_type
        ), axis=1)

    # Filter rows with non-zero be-contained value
    df_rule3 = df_annot_slice[df_annot_slice['Be_Contained'] != 0]
    if len(df_rule3) > 1:
        # Check for overlap and return combined AnnotSV_ID, Datasets, and minimum AF
        bk1s, bk2s = df_rule3['breakend1'].tolist(), df_rule3['breakend2'].tolist()
        num_list = [[int(bk1s[i]), int(bk2s[i])] for i in range(len(bk1s))]
        num_list = get_nset_union(num_list)
        overlap_rate = sum([abs(x[1] - x[0]) for x in num_list]) / abs(bk2 - bk1)
        if overlap_rate >= overlap_threshold:
            return ';'.join(df_rule3['ID'].tolist()), ';'.join(df_rule3['Datasets'].tolist()), min(
                df_rule3['AF'].tolist())

    return 'unknown', 'unknown', 0


# %%
def annotate_judge_sv_af_in_split_gene(df_split_gene_judge_sv, split_gene_public_sv_af_file,
                                       overlap_threshold=None, distance_threshold=None,
                                       column_prefix=''):
    '''
    Annotate new structural variations (SVs) in the need  data frame with gene body information.

    :param df_split_gene_judge_sv: DataFrame containing annotated diagnostic SVs.
    :param split_gene_public_sv_af_file: File path to the gene elements annotation file with allele frequencies (AFs).
    :param overlap_threshold: Threshold for overlap in SV annotations.
    :param distance_threshold: Threshold for breakpoint distance in SV annotations.
    :param column_prefix: Prefix for the new columns to be added in the DataFrame.
    :return: DataFrame with annotated SVs.
    '''
    # Make a copy of the input dataframe
    df_input = df_split_gene_judge_sv.copy()
    split_gene_public_sv_af = pd.read_csv(split_gene_public_sv_af_file, sep='\t')
    # Create a dictionary by column from the gene annotation file
    merge_dict = utils.create_dict_by_column(split_gene_public_sv_af, key_column='#CHR', index_column='txID')

    # Initialize lists to store annotated SV ids, dataset names, and allele frequencies
    annotated_sv_ids = []
    dataset_names = []
    afs = []

    # Iterate through each row in the input dataframe
    for index, row in tqdm(df_input.iterrows()):
        chr_id = row['#CHR']
        tx_id = row['txID']

        # Try to get the annotation slice for the current chromosome and transcript id
        try:
            annot_slice = merge_dict[chr_id].loc[[tx_id]]
        except KeyError:
            # If the annotation slice is not found, store 'unknown' values and continue to the next row
            annotated_sv_ids.append('unknown')
            dataset_names.append('unknown')
            afs.append(0)
            continue
        # Get the SV type from the current row
        sv_type = row['SVTYPE']

        # Add a new column to the annotation slice to check for SV type equality
        annot_slice['SV_TYPE_equal'] = annot_slice['SVTYPE'].apply(
            lambda x: str(is_svtype_match(input_type=x, target_type=sv_type)))

        # Filter the annotation slice based on SV type, region, and element
        annot_slice = annot_slice.loc[(annot_slice['SV_TYPE_equal'] == 'True')
                                      & (annot_slice['Region'] == row['Region'])
                                      & (annot_slice['Element'] == row['Element']), :]

        # Compare and obtain annotated SV id, dataset name, and allele frequency
        asv_id, d_name, af = compare_and_obtain_af(row=row, df_annot_slice=annot_slice,
                                                   overlap_threshold=overlap_threshold,
                                                   distance_threshold=distance_threshold)

        # Append the obtained values to the respective lists
        annotated_sv_ids.append(asv_id)
        dataset_names.append(d_name)
        afs.append(af)

    # Add new columns to the input dataframe with the annotated SV ids, dataset names, and allele frequencies
    df_input[column_prefix + 'AnnotSV_ID'] = annotated_sv_ids
    df_input[column_prefix + 'DataSet_names'] = dataset_names
    df_input[column_prefix + 'AF'] = afs

    return df_input


# %%


def annotate_judge_sv_af_in_regulatory(df_regulatory_judge_sv, regulatory_public_sv_af_file, overlap_threshold=None,
                                       distance_threshold=None,
                                       column_prefix=''):
    """
    Annotate new structural variation (SV) sets with Allele Frequencies (AF) for regulatory elements (RE).

    :param df_regulatory_judge_sv: DataFrame containing the SVs to be annotated.
    :param regulatory_public_sv_af_file: File path to the public database SV sets annotation with RE component AFs.
    :param overlap_threshold: Threshold for mutual overlap.
    :param distance_threshold: Threshold for breakpoint distance.
    :param column_prefix: Prefix for the new columns to be added.
    :return: DataFrame with annotated SVs.
    """
    # Make a copy of the input DataFrame
    df_input = df_regulatory_judge_sv.copy()
    regulatory_public_sv_af = pd.read_csv(regulatory_public_sv_af_file, sep='\t')
    # Create a dictionary by column from the provided file path
    merge_dict = utils.create_dict_by_column(regulatory_public_sv_af, key_column='#CHR', index_column='RE_ID')

    # Initialize lists to store the annotated SVs, dataset names, and allele frequencies
    annot_sv_ids = []
    dataset_names = []
    afs = []

    # Iterate through each row of the input DataFrame
    for index, row in tqdm(df_input.iterrows()):
        chr_id = row['#CHR']
        re_id = row['RE_ID']

        # Attempt to retrieve annotation data for the current chromosome and RE ID
        try:
            df_annot_slice = merge_dict[chr_id].loc[[re_id]]
        except:
            # If data retrieval fails, append default values and continue to the next iteration
            annot_sv_ids.append('unknown')
            dataset_names.append('unknown')
            afs.append(0)
            continue

        # Extract the SV type from the current row
        sv_type = row['SVTYPE']
        # Filter the annotation slice based on SV type and RE type
        df_annot_slice['SV_TYPE_equal'] = df_annot_slice['SVTYPE'].apply(
            lambda x: str(is_svtype_match(input_type=x, target_type=sv_type)))
        df_annot_slice = df_annot_slice.loc[(df_annot_slice['SV_TYPE_equal'] == 'True')
                                            & (df_annot_slice['RE_type'] == row['RE_type']), :]

        # Compare and obtain the annotated SV ID, dataset name, and allele frequency
        asv_id, d_name, af = compare_and_obtain_af(row=row, df_annot_slice=df_annot_slice,
                                                   overlap_threshold=overlap_threshold,
                                                   distance_threshold=distance_threshold)

        # Append the obtained values to the respective lists
        annot_sv_ids.append(asv_id)
        dataset_names.append(d_name)
        afs.append(af)

    # Add new columns to the input DataFrame containing the annotated SVs
    df_input[column_prefix + 'PublicSets_SVID'] = annot_sv_ids
    df_input[column_prefix + 'DataSet_names'] = dataset_names
    df_input[column_prefix + 'AF'] = afs

    # Return the updated DataFrame with annotated SVs
    return df_input


def extract_gene_element_details(ref_gene_canonical_annot_dict):
    """
    Extract and annotate gene element details (like exons, introns, UTRs) from the canonical gene annotation dictionary.

    :param ref_gene_canonical_annot_dict: Dictionary containing reference gene canonical annotations.
    :return: DataFrame with detailed annotations of gene elements.
    """
    out_data = []
    out_col = ['#CHR', 'txID', 'Element', 'Region', 'start', 'end', 'length']
    for chr_i in ref_gene_canonical_annot_dict:
        for r_i, row in ref_gene_canonical_annot_dict[chr_i].iterrows():
            exon_starts = row['exonStarts'].split(',')
            exon_ends = row['exonEnds'].split(',')
            cds_start = int(row['cdsStart'])
            cds_end = int(row['cdsEnd'])
            strand = row['Strand']
            exon_count = int(row['exonCount'])
            gene_start = int(row['txStart'])
            gene_end = int(row['txEnd'])
            all_positions = [gene_start, gene_end]
            all_positions += (exon_starts + exon_ends)
            all_positions += [cds_start, cds_end]
            all_positions = list(filter(None, all_positions))
            all_positions = [int(i) for i in all_positions]  # Remove empty values if any
            all_positions.sort()
            tx_length = 0
            element = ''
            ele_index = 0
            tx_id = row.name
            for l_i in range(len(all_positions) - 1):
                if all_positions[l_i] == all_positions[l_i + 1]:
                    continue
                else:
                    if strand == "+":
                        if str(all_positions[l_i]) in exon_starts:
                            ele_index = str(exon_starts.index(str(all_positions[l_i])) + 1)
                            element = 'exon' + ele_index
                        elif str(all_positions[l_i]) in exon_ends:
                            ele_index = str(exon_ends.index(str(all_positions[l_i])) + 1)
                            element = 'intron' + ele_index
                        elif str(all_positions[l_i + 1]) in exon_ends:
                            ele_index = str(exon_ends.index(str(all_positions[l_i + 1])) + 1)
                            element = 'exon' + ele_index
                        elif str(all_positions[l_i + 1]) in exon_starts:
                            ele_index = str(exon_starts.index(str(all_positions[l_i + 1])) + 1)
                            element = 'intron' + ele_index
                        else:
                            if 'exon' in element and all_positions[l_i] > int(exon_starts[int(ele_index) - 1]) and \
                                    all_positions[
                                        l_i + 1] < int(exon_ends[int(ele_index) - 1]):
                                element = 'exon' + ele_index
                            else:
                                logging.error('Error tx_id: %s' % tx_id)
                                element = 'error'
                                ele_index = -1
                        if int(cds_start) == int(cds_end):
                            region = 'UTR'
                        elif all_positions[l_i + 1] <= int(cds_start):
                            region = '5\'UTR'
                        elif all_positions[l_i] >= int(cds_end):
                            region = '3\'UTR'
                        else:
                            region = 'CDS'
                    else:
                        if str(all_positions[l_i]) in exon_starts:
                            ele_index = str(exon_count - exon_starts.index(str(all_positions[l_i])))
                            element = 'exon' + ele_index
                        elif str(all_positions[l_i]) in exon_ends:
                            ele_index = str(exon_count - 1 - exon_ends.index(str(all_positions[l_i])))
                            element = 'intron' + ele_index
                        elif str(all_positions[l_i + 1]) in exon_ends:
                            ele_index = str(exon_count - exon_ends.index(str(all_positions[l_i + 1])))
                            element = 'exon' + ele_index
                        elif str(all_positions[l_i + 1]) in exon_starts:
                            ele_index = str(exon_count - 1 - exon_starts.index(str(all_positions[l_i + 1])))
                            element = 'intron' + ele_index
                        else:
                            if 'exon' in element and all_positions[l_i] > int(
                                    exon_starts[exon_count - int(ele_index)]) and \
                                    all_positions[
                                        l_i + 1] < int(exon_ends[exon_count - int(ele_index)]):
                                element = 'exon' + ele_index
                            else:
                                logging.error('Error txID: %s' % tx_id)
                                element = 'error'
                                ele_index = -1
                        if int(cds_start) == int(cds_end):
                            region = 'UTR'
                        elif all_positions[l_i + 1] <= int(cds_start):
                            region = '3\'UTR'
                        elif all_positions[l_i] >= int(cds_end):
                            region = '5\'UTR'
                        else:
                            region = 'CDS'
                    length = abs(int(all_positions[l_i + 1]) - int(all_positions[l_i]))
                    out_value = [row['Chr'], tx_id, element, region, all_positions[l_i], all_positions[l_i + 1],
                                 length]
                    out_data.append(dict(zip(out_col, out_value)))
                    tx_length += length
            if tx_length != abs(int(row['txStart']) - int(row['txEnd'])):
                logging.error('Length check failed:  %s' % tx_id)
    df_results = pd.DataFrame(out_data)
    return df_results


# %%
def calculate_end_position_for_sv(row):
    """
    Calculate the end position for a structural variant (SV).

    :param row: A row of SV data.
    :return: The calculated end position for insertion type SV, or the original end position for other types.
    """
    if is_svtype_match(input_type=row['SVTYPE'], target_type='INS'):
        try:
            svlen = abs(int(row['SVLEN']))
            return int(row['START']) + svlen
        except ValueError:
            return row['END']
    else:
        return row['END']
