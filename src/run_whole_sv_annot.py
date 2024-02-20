import os

import pandas as pd

import src.annot_func as annot_func


def whole_sv_annotate(judge_sv_file, public_sv_sets, selected_datasets, overlap_threshold, distance_threshold, output_dir,
                      col_prefix=''):
    """
    Annotates structural variations (SVs) in the target file with allele frequencies (AF) from public SV datasets.

    Args:
    - judge_sv_file (str): Path to the target file containing SVs
    - public_sv_sets (dict): Dictionary mapping dataset keys to file paths of public SV datasets
    - used_datasets (list): List of dataset keys to be used for annotation
    - overlap_threshold (float): Minimum overlap threshold for filtering
    - distance_threshold (int): Maximum distance threshold for filtering
    - output_dir (str): Output directory for temporary and final files
    - col_prefix (str, optional): Prefix for the overall AF column name

    Returns:
    - df_out_file (DataFrame): Annotated SVs dataframe with overall AF column
    """

    # Define the column names for the input BED file
    bed_cols = ['CHR', 'START', 'END', 'ID', 'SVTYPE', 'SVLEN', 'AF']

    # Read the target file into a dataframe
    df_target = pd.read_csv(judge_sv_file, sep='\t', header=None, names=bed_cols)

    # Calculate the END position for each SV and update the dataframe
    df_target['END'] = df_target.apply(lambda x: annot_func.calculate_end_position_for_sv(x), axis=1)

    # Generate a temporary target file
    temp_target_file = os.path.join(output_dir, 'Target.temp.bed')
    df_target.to_csv(temp_target_file, index=False, header=False, sep='\t')

    # Define the column names for the output BED file
    bedtools_cols = ['CHR_1', 'START_1', 'END_1', 'ID_1', 'SVTYPE_1', 'SVLEN_1', 'AF_1',
                     'CHR_2', 'START_2', 'END_2', 'ID_2', 'SVTYPE_2', 'SVLEN_2', 'AF_2',
                     'overlap_bp']

    # Create the output dataframe with selected columns and rename them
    df_out_file = df_target[['ID', 'CHR', 'SVTYPE', 'START', 'END', 'SVLEN']]
    df_out_file.columns = ['ID', '#CHR', 'SVTYPE', 'START', 'END', 'SVLEN']
    df_out_file.drop_duplicates(subset='ID')

    # Iterate through the used datasets for annotation
    for selected_key in selected_datasets:
        # Read the dataset file into a dataframe
        df_key_file = pd.read_csv(public_sv_sets[selected_key], sep='\t', header=None, names=bed_cols)
        df_key_file['END'] = df_key_file.apply(lambda x: annot_func.calculate_end_position_for_sv(x), axis=1)

        # Generate a temporary file for the dataset
        temp_key_file = os.path.join(output_dir, '%s.temp.bed' % selected_key)
        df_key_file.to_csv(temp_key_file, index=False, header=False, sep='\t')

        # Run bedtools intersect to calculate overlaps and create a temporary output file
        bedtools_key_file = os.path.join(output_dir, 'Overall.bedtools.overlap.bed')
        os.system('bedtools intersect -a  %s  -b %s -f %s -r -wo > %s'
                  % (temp_target_file, temp_key_file, overlap_threshold, bedtools_key_file))

        # Read the bedtools output into a dataframe and perform additional filtering
        df_bedtools_temp = pd.read_csv(bedtools_key_file, sep='\t', header=None, names=bedtools_cols)
        #print(df_bedtools_temp)
        if len(df_bedtools_temp)!= 0:
            df_bedtools_temp['SVTYPE_equal'] = df_bedtools_temp.apply(
                lambda x: annot_func.is_svtype_match(input_type=x['SVTYPE_2'], target_type=x['SVTYPE_1']), axis=1)
            df_bedtools_temp['Max_distance'] = df_bedtools_temp.apply(
                lambda x: max(abs(x['END_1'] - x['END_2']), abs(x['START_1'] - x['START_2'])), axis=1)
            df_bedtools_temp = df_bedtools_temp.loc[(df_bedtools_temp['Max_distance'] <= distance_threshold)
                                                    | (df_bedtools_temp['overlap_bp'] == df_bedtools_temp['SVLEN_1']), :]
            df_bedtools_temp = df_bedtools_temp.loc[df_bedtools_temp['SVTYPE_equal'] == True, :]
            df_bedtools_temp = df_bedtools_temp[['ID_1', 'AF_2']]
            df_bedtools_temp.sort_values(by=['AF_2'], ascending=False, inplace=True, ignore_index=True)
            df_bedtools_temp.drop_duplicates(subset=['ID_1'], keep='first', inplace=True, ignore_index=True)
            df_bedtools_temp.columns = ['ID', selected_key + '_AF']
            df_bedtools_temp['ID'] = df_bedtools_temp['ID'].astype(str)
        else:
            df_bedtools_temp = pd.DataFrame(columns=['ID', selected_key + '_AF'])
        # Merge the annotated AF results into the output dataframe
        df_out_file['ID'] = df_out_file['ID'].astype(str)
        df_out_file = pd.merge(df_out_file, df_bedtools_temp, on='ID', how='left')
        # Clean up temporary files
        os.system('rm %s' % temp_key_file)
        os.system('rm %s' % bedtools_key_file)

    # Fill in missing values with 0 and calculate the overall AF
    df_out_file.fillna(0, inplace=True)
    af_cols = [i + '_AF' for i in selected_datasets]
    if col_prefix == '':
        df_out_file['Whole_SV_AF'] = df_out_file.apply(lambda x: max(x[af_cols]), axis=1)
    else:
        df_out_file[col_prefix + '_Whole_SV_AF'] = df_out_file.apply(lambda x: max(x[af_cols]), axis=1)

    # Return the annotated SVs dataframe
    return df_out_file


def run_whole_sv_annotate(judge_sv_file, public_sv_sets, selected_datasets,
                          overlap_threshold, distance_threshold, output_dir, out_prefix=''):
    df_out = whole_sv_annotate(judge_sv_file=judge_sv_file,
                               public_sv_sets=public_sv_sets,
                               selected_datasets=selected_datasets,
                               overlap_threshold=overlap_threshold,
                               distance_threshold=distance_threshold,
                               output_dir=output_dir,
                               col_prefix='')
    out_file = os.path.join(output_dir, out_prefix + '.Whole_SV_AF.tsv')
    df_out.to_csv(out_file, sep='\t', index=False)
    return out_file
