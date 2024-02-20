import logging

import math
import pandas as pd


# %%
def weight_by_af(af, common_af=0.01, max_weight_by_af=1, max_decline_rate=2000):
    """
    Calculate the weight based on allele frequency (AF) with a logistic function.
    Ensures that math range errors are avoided for large AF values.

    :param af: Allele Frequency (AF).
    :param common_af: AF threshold where the decline becomes significant.
    :param max_weight_by_af: Maximum weight when AF is 0.
    :param max_decline_rate: Maximum decline rate of the curve.
    :return: Calculated weight based on AF.
    """
    # Adjust the formula to make sure the value is 1 at af=0
    adjusted_af = af - common_af * 0.9
    decline_factor = adjusted_af * max_decline_rate
    try:
        return max_weight_by_af / (1 + math.exp(decline_factor))
    except OverflowError:
        return 0


# %%
def calculate_utr_lengths(gene_annotation_dict):
    """
    Calculate the lengths of the 5' UTR and 3' UTR for each gene.

    :param gene_annotation_dict: A dictionary containing gene reference annotations.
    :return: A dictionary mapping transcript IDs to their respective UTR lengths,
             start and end positions.
    """
    # Initialize a dictionary to store the UTR lengths
    utr_lengths_dict = {}

    # Iterate through each chromosome in the gene annotations
    for chromosome in gene_annotation_dict:
        # Iterate through each transcript ID and its corresponding row in the gene annotations
        for txID, row in gene_annotation_dict[chromosome].iterrows():
            # Check if the transcript has a coding sequence
            if row['cdsStart'] != row['cdsEnd']:
                # Calculate variable UTR length if the transcript has a coding sequence
                utr_5_length, utr_3_length = calculate_variable_utr_length(row)
            else:
                # Handle the case where there is no coding sequence
                utr_5_length = int(row['cdsStart']) - int(row['txStart'])
                utr_3_length = utr_5_length

            # Determine the UTR lengths based on the strand of the transcript
            if row['Strand'] == '+':
                utr_lengths = [utr_5_length, utr_3_length]
            else:
                utr_lengths = [utr_3_length, utr_5_length]

            # Store the UTR lengths and transcript start and end positions in the dictionary
            utr_lengths_dict[txID] = utr_lengths + [row['txStart'], row['txEnd']]

    return utr_lengths_dict


def calculate_variable_utr_length(row):
    """
    Calculate UTR lengths for genes with variable UTRs.

    :param row: A row from the gene annotations dataframe.
    :return: Lengths of the 5' and 3' UTR.
    """
    utr_1_length = 0
    utr_2_length = 0
    found_cds_start = False
    found_cds_end = False
    exon_starts = row['exonStarts'].split(',')[:-1]
    exon_ends = row['exonEnds'].split(',')[:-1]

    for exon_start, exon_end in zip(exon_starts, exon_ends):
        if not found_cds_start:
            if int(exon_start) <= int(row['cdsStart']) <= int(exon_end):
                utr_1_length += int(row['cdsStart']) - int(exon_start)
                found_cds_start = True
            else:
                utr_1_length += int(exon_end) - int(exon_start)

        if int(exon_start) <= int(row['cdsEnd']) <= int(exon_end):
            utr_2_length += int(exon_end) - int(row['cdsEnd'])
            found_cds_end = True
        elif found_cds_end:
            utr_2_length += int(exon_end) - int(exon_start)

    return utr_1_length, utr_2_length


# %%
def calculate_sv_utr_length(sv_annot_row, gene_annotation_dict):
    """
    Calculate the length of the UTR covered by a structural variant (SV).

    :param sv_annot_row: A row from a dataframe representing SV annotation.
    :param gene_annotation_dict: A dictionary of gene annotations.
    :return: Length of the UTR covered by the SV.
    """
    sv_type = sv_annot_row['SVTYPE']
    if sv_type == 'INS':
        return int(sv_annot_row['SVLEN'])

    breakend1, breakend2 = int(sv_annot_row['breakend1']), int(sv_annot_row['breakend2'])
    overlap_regions = sv_annot_row['OverlapRegions']
    gene_row = gene_annotation_dict.loc[sv_annot_row['txID']]

    utr_bk1, utr_bk2 = determine_utr_boundaries(gene_row, overlap_regions)

    if not (utr_bk1 <= breakend1 <= breakend2 <= utr_bk2):
        logging.error(f'Error! SV ID: {sv_annot_row["ID"]}, Transcript ID: {sv_annot_row["txID"]}')
        return 0

    sv_utr_length = calculate_overlapped_utr_length(gene_row, sv_annot_row, utr_bk1, utr_bk2)
    return sv_utr_length


def determine_utr_boundaries(gene_row, overlap_regions):
    """
    Determine the UTR boundaries for a gene.

    :param gene_row: A row from the gene annotations dataframe.
    :param overlap_regions: The region of overlap with the SV.
    :return: The start and end positions of the UTR.
    """
    if "5'UTR" == overlap_regions:
        utr_start = gene_row['txStart'] if gene_row['Strand'] == '+' else gene_row['cdsEnd']
        utr_end = gene_row['cdsStart'] if gene_row['Strand'] == '+' else gene_row['txEnd']
    else:
        utr_start = gene_row['txStart'] if gene_row['Strand'] == '-' else gene_row['cdsEnd']
        utr_end = gene_row['cdsStart'] if gene_row['Strand'] == '-' else gene_row['txEnd']
    return int(utr_start), int(utr_end)


def calculate_overlapped_utr_length(gene_row, sv_annot_row, utr_start, utr_end):
    """
    Calculate the length of the UTR overlapped by the SV.

    :param gene_row: A row from the gene annotations dataframe.
    :param sv_annot_row: A row from the SV annotations dataframe.
    :param utr_start: Start position of the UTR.
    :param utr_end: End position of the UTR.
    :return: Length of the overlapped UTR.
    """
    sv_utr_length = 0
    breakend1, breakend2 = int(sv_annot_row['breakend1']), int(sv_annot_row['breakend2'])
    for name in sv_annot_row['OverlapElements'].split(';'):
        if 'exon' in name:
            exon_start, exon_end = gene_row['elements_dict'][name]
            count_start = max(exon_start, breakend1)
            count_end = min(exon_end, breakend2)
            sv_utr_length += min(count_end, utr_end) - max(count_start, utr_start)
    return sv_utr_length


# %%
def apply_utr_score(x, utr_length, weight_for_utr, common_af):
    """
    Calculate the UTR score based on various factors.

    :param x: A dictionary or similar object with SV data.
    :param utr_length: Length of the UTR region.
    :param weight_for_utr: Weighting factor for the score calculation.
    :param common_af: Common Allele Frequency
    :return: Calculated UTR score.
    """
    # Default value for Function_Element_Final_weight
    weight_fge = x.get('Function_Element_Final_weight', 1)

    # Check for high allele frequency
    if x.get('High_AF(LRS_sets)', 'No') == 'Yes':
        return 0

    # Calculate the score based on SV length in the region
    sv_len_in_region = int(x.get('SV_LEN_in Region', 0))
    af = x.get('AF', 0)  # Default allele frequency

    # Adjust score calculation based on SV length in the region
    if sv_len_in_region != 0:
        score = abs(weight_for_utr * sv_len_in_region) / utr_length * weight_by_af(af, common_af=common_af) * weight_fge
    else:
        score = abs(weight_for_utr * utr_length / 2) / utr_length * weight_by_af(af, common_af=common_af) * weight_fge

    return score


# %%
def give_utr_sv_scores(df_sv_split_gene_utr, sv_type, utr_length_dict):
    """
    Calculate scores based on UTR lengths.

    :param df_sv_split_gene_utr: Input DataFrame.
    :param sv_type: Type of structural variant.
    :param utr_length_dict: Dictionary of gene UTR lengths.
    :return: DataFrame with calculated scores.
    """
    # Initialize output data dictionary
    out_data = {
        'txID': [], 'geneID': [], 'Region': [],
        'GeneBody_AF': [], 'GeneBody_weight': [], 'GeneBody_Score': []
    }

    utr_index = {'5\'UTR': 0, '3\'UTR': 1, 'UTR': 0}

    for tx_id, group_df in df_sv_split_gene_utr.groupby('txID'):
        for utr_class in group_df['Region'].unique():
            utr_class_df = group_df[(group_df['Region'] == utr_class) & (group_df['Element'].str.contains('exon'))]
            if utr_class_df.empty:
                append_empty_class_data(out_data, group_df, tx_id, utr_class)
                continue

            utr_len = utr_length_dict[tx_id][utr_index[utr_class]]
            calculate_class_scores(utr_class_df, utr_len, sv_type, out_data, tx_id, utr_class)

    return pd.DataFrame(out_data).set_index('txID')


def append_empty_class_data(out_data, group_df, tx_id, utr_class):
    """ Append data for classes with no exon elements. """
    out_data['txID'].append(tx_id)
    out_data['geneID'].append(group_df['geneID'].values[0])
    out_data['Region'].append(utr_class)
    out_data['GeneBody_AF'].append(group_df['AF'].min())
    out_data['GeneBody_weight'].append(group_df['Function_Element_Final_weight'].values[0])
    out_data['GeneBody_Score'].append(0)


def calculate_class_scores(utr_class_df, utr_len, sv_type, out_data, tx_id, utr_class):
    """ Calculate scores for each class. """
    if sv_type in ['DEL', 'DUP', 'INV']:
        utr_class_df['SV_LEN_in Region'] = utr_class_df.apply(lambda x: abs(x['breakend1'] - x['breakend2']), axis=1)
    else:
        utr_class_df['SV_LEN_in Region'] = utr_class_df['SVLEN'].apply(lambda x: int(x) if x != 'unknown' else 0)

    utr_class_df['Scores for UTR'] = utr_class_df.apply(
        lambda x: apply_utr_score_func(x, utr_len, 1 if sv_type in ['DEL', 'DUP', 'INV'] else 1 / 2), axis=1)

    out_data['txID'].append(tx_id)
    out_data['geneID'].append(utr_class_df['geneID'].values[0])
    out_data['Region'].append(utr_class)
    out_data['GeneBody_AF'].append(utr_class_df['AF'].min())
    out_data['GeneBody_weight'].append(
        utr_class_df.loc[utr_class_df['AF'] == utr_class_df['AF'].min(), 'Function_Element_Final_weight'].values[0])
    out_data['GeneBody_Score'].append(min(utr_class_df['Scores for UTR'].sum(), 0.5))


def apply_utr_score_func(x, utr_len, weight):
    # Implement the logic of Apply_UTR_Score_Func function here
    pass


# %%
def give_split_gene_sv_scores(df_sv_gene_split, sv_type, common_af):
    """
    Assign scores to genes based on the split data.
    :param df_sv_gene_split: DataFrame containing gene split data.
    :param sv_type: Type of structural variant.
    :param common_af: Common AF threshold.
    :return: DataFrame with calculated scores.
    """
    df_sv_gene_split['GeneBody_Score'] = df_sv_gene_split.apply(process_split_gene_sv_row, sv_type=sv_type,
                                                                common_af=common_af, axis=1)
    df_sv_gene_split['GeneBody_AF'] = df_sv_gene_split['AF'].apply(float)
    df_sv_gene_split['GeneBody_weight'] = df_sv_gene_split['Function_Element_Final_weight'].fillna(1)

    return df_sv_gene_split


def process_split_gene_sv_row(row, sv_type, common_af):
    """
    Calculate the score for a row in the DataFrame.

    :param row: A row from DataFrame.
    :param sv_type: Type of structural variant.
    :param common_af: Common AF threshold.
    :return: Score calculated based on the row data.
    """
    element = row['Element']
    region = row['Region']
    af = float(row.get('AF', 0))
    w_af = weight_by_af(af, common_af=common_af)
    w_fge = row.get('Function_Element_Final_weight', 1)
    if sv_type in ['TRA', 'BND']:
        score = 1 if region == 'UTR' else 1
    else:
        if 'intron' in element:
            score = 0.1
        elif 'UTR' in region:
            score = 0.5
        else:
            score = 1
    return score * w_af * w_fge


# %%
def integrate_regulatory_sv_score(df_regulatory_sv, regulatory_element_types, regulatory_element_weight, common_af):
    """
    Integrate the impact of SV on regulatory regions.

    :param df_regulatory_sv: DataFrame with regulatory element data.
    :param regulatory_element_types: List of regulatory element types to consider.
    :param regulatory_element_weight: Dictionary of weights for each regulatory element type.
    :param common_af: Common AF threshold.
    :return: DataFrame with integrated data.
    """
    # Initialize output DataFrame
    df_integrate_regulatory_sv_score = pd.DataFrame()
    gene_ids = df_regulatory_sv['geneID'].unique().tolist()
    df_integrate_regulatory_sv_score['geneID'] = gene_ids
    for re_key in regulatory_element_types:
        df_integrate_regulatory_sv_score[f'Overlap_{re_key}'] = 'No'
        df_integrate_regulatory_sv_score[f'{re_key}_AF'] = 0
        df_integrate_regulatory_sv_score[f'{re_key}_weight'] = 1
        df_integrate_regulatory_sv_score[f'{re_key}_Score'] = 0

    # Update df_integrate_regulatory_sv_score based on df_regulatory_sv
    for index, row in df_regulatory_sv.iterrows():
        gene_id = row['geneID']
        re_type = row['RE_type']
        af = row.get('AF', 0)
        w_fre = row.get('Function_Element_Final_weight', 1)

        for re_key in regulatory_element_types:
            if re_key in re_type:
                df_integrate_regulatory_sv_score.loc[
                    df_integrate_regulatory_sv_score['geneID'] == gene_id, f'Overlap_{re_key}'] = 'Yes'
                df_integrate_regulatory_sv_score.loc[
                    df_integrate_regulatory_sv_score['geneID'] == gene_id, f'{re_key}_AF'] = af
                df_integrate_regulatory_sv_score.loc[
                    df_integrate_regulatory_sv_score['geneID'] == gene_id, f'{re_key}_weight'] = w_fre
                score = regulatory_element_weight[re_key] * weight_by_af(af, common_af=common_af) * w_fre
                df_integrate_regulatory_sv_score.loc[
                    df_integrate_regulatory_sv_score['geneID'] == gene_id, f'{re_key}_Score'] = score
    return df_integrate_regulatory_sv_score


# %%
def give_whole_gene_conservation_scores(df_data, sv_type, gene_hi_lof, regulatory_element_types):
    """
    Calculate composite scores and additional metrics for SV conservation,
    rounding scores in the DataFrame.

    Args:
    - df_data: DataFrame containing the data
    - sv_type: type of structural variation
    - gene_hi_lof: indicator for gene high or low impact
    - consider_re_keys: indicator for considering regulatory elements keys

    Returns:
    - DataFrame with composite scores and additional metrics
    """

    # Round scores in the DataFrame
    score_columns = [col for col in df_data.columns if col.endswith('_Score')]
    for col in score_columns:
        df_data[col] = df_data[col].astype(float).apply(lambda x: round(x, 2))

    # Calculate composite scores and additional metrics
    scores = []
    loeuf_bin = []
    hi_bin = []
    ti_bin = []
    for r_i, row in df_data.iterrows():
        # Calculate the scores for each regulatory element type
        re_scores = [row.get(f'{re_key}_Score', 0) for re_key in regulatory_element_types]
        # Calculate the transcript score based on frameshift and SV type
        tx_s = row['GeneBody_Score'] - 0.1 if row['Frameshit'] == 'No' and sv_type != 'INS' else row['GeneBody_Score']
        tx_s = max(tx_s - 0.4, 0) if sv_type == 'INS' else tx_s

        # Get LoF, hi, and tolerated scores for the gene
        lof_hi_ti = gene_hi_lof.get(row['geneID'], [9, 1, 1])

        # Calculate the sum of all scores
        all_s = tx_s + sum(re_scores)

        # If sum of all scores is very low, return 0
        if all_s <= 0.01:
            scores.append(0)
            loeuf_bin.append(lof_hi_ti[0])
            hi_bin.append(lof_hi_ti[1])
            ti_bin.append(lof_hi_ti[2])
            continue
        # Additional score adjustments based on transcript score
        if tx_s > 0:
            all_s += 0.45 - 0.05 * lof_hi_ti[0]
        elif row['Frameshit'] == 'Yes' and row['relationship'] == 'SV in Intronic':
            all_s += 0.27 - 0.03 * lof_hi_ti[0]

        # Adjustments specific to SVTYPE
        if sv_type == 'DUP':
            all_s += 0.5 - 0.5 * lof_hi_ti[2]
        else:
            all_s += 0.5 - 0.5 * lof_hi_ti[1]

        # Return the maximum of 0 and the composite score, and the LoF, hi, and tolerated scores
        scores.append(max(0, all_s))
        loeuf_bin.append(lof_hi_ti[0])
        hi_bin.append(lof_hi_ti[1])
        ti_bin.append(lof_hi_ti[2])

    df_data['Score'] = scores
    df_data['LOEUF_bin'] = loeuf_bin
    df_data['HI_bin'] = hi_bin
    df_data['TI_bin'] = ti_bin
    return df_data


def adjust_scores_for_sv(df_data, chrom, sv_type, sv_breakend1, sv_breakend2, overlap_tad, gene_annotation_dict):
    """
    Adjust scores for special structural variation (SV) conditions.

    Args:
    df_data (pd.DataFrame): The input dataframe containing the data to be adjusted.
    chrom (str): The chromosome of the gene.
    sv_type (str): The type of structural variation (e.g., DUP, INV).
    sv_breakend1 (int): The start position of the structural variation.
    sv_breakend2 (int): The end position of the structural variation.
    overlap_tad (str): Whether the structural variation overlaps with a TAD (Topologically Associating Domain).
    gene_annotation_dict (dict): A dictionary containing gene annotations for each chromosome.

    Returns:
    pd.DataFrame: The input dataframe with adjusted scores.
    """

    # Initialize lists to store adjusted scores
    tx_scores = []
    op_scores = []

    # Iterate through each row in the dataframe
    for index, row in df_data.iterrows():
        # Extract relevant data from the row
        tx_id = row['txID']
        relationship = row['relationship']
        try:
            gene_data = gene_annotation_dict[chrom].loc[tx_id]
        except KeyError:
            tx_scores.append(row['GeneBody_Score'])
            op_scores.append(row['Promoter_Score'])
            continue
        tx_bk1, tx_bk2 = int(gene_data['txStart']), int(gene_data['txEnd'])

        # Extend the boundary for strand-specific cases
        if gene_data['Strand'] == '+':
            tx_bk1 -= 2000
        else:
            tx_bk2 += 2000

        # Adjust scores based on SV type and relationship
        if sv_type == 'DUP' and relationship == 'SV Partial overlap Tx':
            # Check if the SV partially overlaps with the gene
            is_within_bounds = sv_breakend1 <= tx_bk1 < sv_breakend2 < tx_bk2 or tx_bk1 < sv_breakend1 < tx_bk2 <= sv_breakend2
            # Set score adjustment based on overlap condition
            score_adjustment = 0 if is_within_bounds else row['GeneBody_Score']

        elif sv_type == 'INV' and relationship == 'SV contains all Tx':
            # Check if the SV contains the entire gene and overlaps with TAD
            is_within_bounds = sv_breakend1 <= tx_bk1 <= tx_bk2 <= sv_breakend2
            # Set score adjustment based on overlap and TAD condition
            score_adjustment = 0 if is_within_bounds and overlap_tad == 'No' else 0.3

        else:
            # Default score adjustment
            score_adjustment = row['GeneBody_Score']

        # Append adjusted scores to the lists
        tx_scores.append(score_adjustment)
        op_scores.append(row['Promoter_Score'])

    # Update the dataframe with the adjusted scores
    df_data['GeneBody_Score'] = tx_scores
    df_data['Promoter_Score'] = op_scores

    return df_data


# %%
def sum_score(df_data):
    """
    Calculate weighted sum of scores based on certain conditions.

    Args:
    df_data: DataFrame containing 'Score', 'relationship', 'Overlap_Promoter', 'Overlap_Enhancer' columns

    Returns:
    Weighted sum of scores
    """
    # Exclude rows with Score 0
    df_data_no_0 = df_data.loc[df_data['Score'] != 0, :]

    # If no rows with non-zero Score, return 0
    if len(df_data_no_0) == 0:
        return 0

    # Group by relationship and Overlap columns, then calculate count and sum of scores
    df_count = df_data_no_0.groupby(by=['relationship', 'Overlap_Promoter', 'Overlap_Enhancer'])['Score'].count()
    df_sum = df_data_no_0.groupby(by=['relationship', 'Overlap_Promoter', 'Overlap_Enhancer'])['Score'].sum()

    # Calculate weighted sum based on count and sum of scores
    sum_s = 0
    index_list = df_sum.index.tolist()
    for in_i in index_list:
        N = df_count.loc[in_i]
        if N > 5:
            sum_s = (df_sum.loc[in_i] / N) * 5 + sum_s
        else:
            sum_s = sum_s + df_sum.loc[in_i]

            # Round the sum to 2 decimal places
            sum_s = round(sum_s, 2)
    return sum_s


def max_score(df_data):
    """
    Calculate the maximum score from the given dataframe.

    Args:
    df_data (pandas.DataFrame): The input dataframe.

    Returns:
    int: The maximum score from the dataframe, or 0 if the dataframe is empty.
    """
    #print(df_data)
    return df_data['Score'].max() if not df_data.empty else 0


# %%
def calculate_score_func(row, intergenic_score, common_af):
    """
    Calculate a score based on the allele frequency and high allele frequency status.

    Args:
    row (pandas.Series): A row from a DataFrame.
    intergenic_score (float): Base score for intergenic regions.
    common_af (float): Common (overall) Allele Frequency.
    Returns:
    float: Calculated score.
    """
    # Get the overall allele frequency from the row, default to 0 if not present
    overall_af = row.get('Whole_SV_AF', 0)

    # Calculate the score by weighting the intergenic score by the overall allele frequency
    return intergenic_score * weight_by_af(overall_af, common_af=common_af)


def give_intergenic_sv_score(df_sv_scores, df_out_sv_score_info, target_overall_af, intergenic_score, common_af):
    # Identify SVs that are only present in intergenic regions
    existing_ids = set(df_sv_scores['ID'])
    target_overall_af['Only_in_Intergenic(No Function)'] = target_overall_af['ID'].apply(
        lambda x: 'Yes' if x not in existing_ids else 'No'
    )

    # Filter for intergenic SVs and calculate their scores
    intergenic_svs = target_overall_af[target_overall_af['Only_in_Intergenic(No Function)'] == 'Yes']
    intergenic_svs['Score'] = intergenic_svs.apply(
        lambda x: calculate_score_func(x, intergenic_score, common_af), axis=1
    )

    # Merge the new scores with the existing score information
    df_out_sv_score_info_updated = pd.merge(
        df_out_sv_score_info,
        target_overall_af[['ID', 'Only_in_Intergenic(No Function)', 'Whole_SV_AF']],
        on='ID', how='left'
    )

    # Concatenate the new intergenic SV scores
    df_out_sv_score_info_updated = pd.concat([
        df_out_sv_score_info_updated,
        intergenic_svs[['ID', 'Only_in_Intergenic(No Function)', 'Whole_SV_AF', 'Score']]
    ])
    # Fill NaN values appropriately
    fill_values = {'Overlap': 'No', 'Score': 0, 'AF': 0, 'relationship': 'No Overlap'}
    for col in df_out_sv_score_info_updated.columns:
        fill_na = False
        for col_sub, fill_val in fill_values.items():
            if col_sub in col:
                df_out_sv_score_info_updated[col] = df_out_sv_score_info_updated[col].fillna(fill_val)
                fill_na=True
        if not fill_na:
            df_out_sv_score_info_updated[col] = df_out_sv_score_info_updated[col].fillna('NaN')


    # Update the final score
    intergenic_svs['Final_Score'] = intergenic_svs.pop('Score')

    # Append the new SV scores to the existing scores
    df_sv_scores_updated = df_sv_scores.append(intergenic_svs[df_sv_scores.columns])

    return df_sv_scores_updated, df_out_sv_score_info_updated
