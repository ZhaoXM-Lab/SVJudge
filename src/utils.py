# from multiprocessing.dummy import Pool as ThreadPool

import logging
import os
import time

import pandas as pd
import vcf


# %%
def vcf2bed(vcf_dir, output_bed, public_sv=False, include_TRA=False):
    vcf_reader = vcf.Reader(open(vcf_dir, 'r'), )
    TRA_record = []
    with open(output_bed, 'w') as bed:
        for record in vcf_reader:
            sv_id = record.ID
            chrom = str(record.CHROM)
            if not chrom.startswith('chr'):
                chrom = 'chr' + chrom
            pos = int(record.POS)

            sv_type = record.INFO.get('SVTYPE')
            if type(sv_type) == list:
                sv_type = sv_type[0]
            sv_length = record.INFO.get('SVLEN', 1)
            if type(sv_length) == list:
                sv_length = sv_length[0]
            else:
                sv_length = abs(sv_length)
            end = record.INFO.get('END', pos + abs(sv_length))
            if type(end) == list:
                end = int(end[0])
            else:
                end = int(end)
            if public_sv:
                af = record.INFO.get('AF', 'unknown')
                if type(af) == list:
                    af = af[0]
            else:
                af = 'unknown'
            if sv_type == 'TRA' or sv_type == 'BND':
                if include_TRA:
                    if 'CHR2' in record.INFO:
                        chrom2 = str(record.INFO.get('CHR2'))
                        if not chrom2.startswith('chr'):
                            chrom2 = 'chr' + chrom2
                        pos1 = pos
                        pos2 = end
                        end1 = pos1 + 1
                        end2 = pos2 + 1

                    else:
                        alt = str(record.ALT[0])
                        bnd_2 = alt.replace('[', '').replace(']', '').replace('N', '').split(':')
                        chrom2 = bnd_2[0]
                        pos2 = int(bnd_2[1])
                        pos1 = int(pos)
                        end1 = pos1 + 1
                        end2 = pos2 + 1
                    record1 = '%s:%s-%s:%s' % (chrom, pos1, chrom2, pos2)
                    record2 = '%s:%s-%s:%s' % (chrom2, pos2, chrom, pos1)
                    if record1 not in TRA_record and record2 not in TRA_record:
                        TRA_record.append(record1)
                        bed.write('%s\t%d\t%d\t%s\t%s\t%s\t%s\n' % (chrom, pos1, end1, sv_id, sv_type, sv_length, af))
                        bed.write('%s\t%d\t%d\t%s\t%s\t%s\t%s\n' % (chrom2, pos2, end2, sv_id, sv_type, sv_length, af))
            else:
                bed.write('%s\t%d\t%d\t%s\t%s\t%s\t%s\n' % (chrom, pos, end, sv_id, sv_type, sv_length, af))
    return output_bed


# %%

def load_annotation_files(annotation_files_dict, output_dir):
    """
    Load gene and regulatory regions annotations.

    Args:
    annotation_files_dict (dict): A dictionary of file paths.
    output_dir (str): Directory for output files.

    Returns:
    tuple: A tuple containing gene information dictionary, path to gene reference temp BED file,
           regulatory elements reference BED files, and regulatory elements keys.
    """

    # Read gene annotation file and create a temporary BED file
    gene_annotation = pd.read_csv(annotation_files_dict['Gene'], header=None, sep='\t')
    gene_annotation.columns = ['txID', 'Chr', 'Strand', 'txStart', 'txEnd', 'cdsStart',
                               'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'Score',
                               'geneID', 'cdsStartStat', 'cdsEndStat', 'exonFrames']

    df_gene_reference_bed = gene_annotation[['Chr', 'txStart', 'txEnd', 'txID']]
    gene_reference_bed = os.path.join(output_dir, 'GeneRef.temp.bed')
    df_gene_reference_bed.to_csv(gene_reference_bed, sep='\t', index=False, header=False)

    # Create a dictionary from the DataFrame for faster lookup
    gene_annotation_dict = create_dict_by_column(df=gene_annotation, key_column='Chr', index_column='txID')

    # Load regulatory annotations files
    regulatory_reference_beds, regulatory_element_types = [], []
    for key in annotation_files_dict:
        if key != "Gene":
            regulatory_reference_beds.append(annotation_files_dict[key])
            regulatory_element_types.append(key)

    return gene_annotation_dict, gene_reference_bed, regulatory_reference_beds, regulatory_element_types


# %%

def calculate_runtimes(log_str, start_time):
    """
    Calculate and print the runtime from the start to the current time.

    :param log_str: str - The keyword for the operation, used for identifying log information.
    :param start_time: float - The starting point of time, usually the timestamp when the process started.
    :return: tuple - Current timestamp and the formatted string of elapsed time.
    """
    # Calculate elapsed time
    elapsed_time = time.time() - start_time
    # Convert elapsed time to hours, minutes, and seconds
    hours, remainder = divmod(elapsed_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    # Format the elapsed time as a string
    formatted_str = f"{log_str} (Runtime Duration:{int(hours):02d} h {int(minutes):02d} min {int(seconds):02d} s)"
    # Center the formatted string
    formatted_str = formatted_str.center(50, '*')
    # Print the formatted string
    logging.info(formatted_str)
    # Return current timestamp and formatted string
    return time.time(), formatted_str


# %%
def create_dict_by_column(df, key_column=None, index_column=None):
    """
    Create a dictionary from a DataFrame for faster lookup, split by a specified column.

    Args:
    df (DataFrame): The input DataFrame.
    key_column (str): The column to be used as a key for the dictionary.
    index_column (str): The column to set as the index of the DataFrame.

    Returns:
    dict: A dictionary where each key corresponds to a subset of the DataFrame.
    """
    # Convert key_column to string type for consistent comparison
    if key_column:
        df[key_column] = df[key_column].astype(str)

        # Set index_column as the index of the DataFrame if provided
        if index_column:
            df.set_index(index_column, inplace=True)

        # Create a dictionary where each key corresponds to a subset of the DataFrame
        unique_values = df[key_column].unique()
        dict_by_column = {key: df[df[key_column] == key] for key in unique_values}
    else:
        # If no key_column provided, use 0 as the key for the entire DataFrame
        dict_by_column = {0: df}
    return dict_by_column
#%%