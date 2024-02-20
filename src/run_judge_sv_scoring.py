# %%
import logging
import os

import pandas as pd
from tqdm import tqdm

import src.scoring_func as scoring_func
import src.utils as utils


# %%
def give_score(df_judge_sv_whole_gene, df_judge_sv_split_gene, df_judge_sv_regulatory, df_judge_whole_sv_af,
               gene_hi_lof, gene_element_sv_weight_dict, regulatory_element_sv_weight_dict, gene_annotation_dict,
               regulatory_element_types, regulatory_element_weight,
               element_sv_weight_source, common_af, mode='Max'):
    """
    Calculates the overall scores for structural variants (SVs), taking into account whole-gene annotations, split-gene annotations, regulatory element annotations, and allele frequencies for the whole SV.

    Parameters:
    df_judge_sv_whole_gene (DataFrame): DataFrame containing whole-gene annotation information.
    df_judge_sv_split_gene (DataFrame): DataFrame containing split-gene annotation information.
    df_judge_sv_regulatory (DataFrame): DataFrame containing regulatory element annotation information.
    df_judge_whole_sv_af (DataFrame): DataFrame containing allele frequency information for the whole SV.
    gene_hi_lof (Dict/File): Dictionary or file path containing gene HI, TI, LOF information.
    gene_element_sv_weight_dict (Dict): Dictionary of gene element weight data.
    regulatory_element_sv_weight_dict (Dict): Dictionary of regulatory element weight data.
    gene_annotation_dict (Dict): Dictionary of gene annotation information.
    regulatory_element_types (List/Set): Types of regulatory elements considered.
    regulatory_element_weight (Dict): Weight dictionary for regulatory elements.
    element_sv_weight_source (List, optional): Sources of element weight data, default is ['PGG_SV_LRS'].
    common_af (float): Common allele frequency
    mode (str, optional): Integration mode, either 'Max' or 'Sum', default is 'Max'.

    Returns:
    Tuple: Contains two DataFrames, the first with SV score information, and the second with SV scores.

    This function evaluates the potential impact of structural variants by merging and processing datasets from various sources. It utilizes multiple weights and annotation information to calculate the final score for each SV, aiding further analysis and interpretation.
    """
    df_gene_element_sv_weight = pd.DataFrame()
    df_regulatory_element_sv_weight = pd.DataFrame()
    mean_col = []
    for data in element_sv_weight_source:
        mean_col.append('%s_weight' % data)
        gene_element_sv_weight_dict[data]['#CHR'] = gene_element_sv_weight_dict[data]['#CHR'].astype(str)
        regulatory_element_sv_weight_dict[data]['#CHR'] = regulatory_element_sv_weight_dict[data]['#CHR'].astype(str)
        gene_element_sv_weight_dict[data].rename(columns={'weight': '%s_weight' % data}, inplace=True)
        regulatory_element_sv_weight_dict[data].rename(columns={'weight': '%s_weight' % data}, inplace=True)
        if df_gene_element_sv_weight.empty:
            df_gene_element_sv_weight = gene_element_sv_weight_dict[data][
                ['#CHR', 'geneID', 'txID', 'Element', 'Region', '%s_weight' % data]]
            df_regulatory_element_sv_weight = regulatory_element_sv_weight_dict[data][
                ['#CHR', 'geneID', 'RE_ID', '%s_weight' % data]]
        else:
            df_gene_element_sv_weight = pd.merge(df_gene_element_sv_weight, gene_element_sv_weight_dict[data][
                ['#CHR', 'geneID', 'txID', 'Element', 'Region', '%s_weight' % data]],
                                                 on=['#CHR', 'geneID', 'txID', 'Element', 'Region'],
                                                 how='outer').reset_index(drop=True)
            df_regulatory_element_sv_weight = pd.merge(df_regulatory_element_sv_weight,
                                                       regulatory_element_sv_weight_dict[data][
                                                           ['#CHR', 'geneID', 'RE_ID', '%s_weight' % data]],
                                                       on=['#CHR', 'geneID', 'RE_ID'],
                                                       how='outer').reset_index(drop=True)
    df_regulatory_element_sv_weight.fillna(1, inplace=True)
    df_gene_element_sv_weight.fillna(1, inplace=True)
    df_gene_element_sv_weight['Function_Element_Final_weight'] = df_gene_element_sv_weight[mean_col].mean(axis=1).apply(
        lambda x: round(x, 2))
    df_regulatory_element_sv_weight['Function_Element_Final_weight'] = df_regulatory_element_sv_weight[mean_col].mean(
        axis=1).apply(
        lambda x: round(x, 2))
    df_whole_gene_sv = df_judge_sv_whole_gene.copy()
    # Gene_elements
    df_split_gene_sv = df_judge_sv_split_gene.copy()
    df_split_gene_sv['#CHR'] = df_split_gene_sv['#CHR'].astype(str)
    df_gene_element_sv_weight['#CHR'] = df_gene_element_sv_weight['#CHR'].astype(str)
    df_split_gene_sv = pd.merge(df_split_gene_sv,
                                df_gene_element_sv_weight[
                                    ['#CHR', 'txID', 'Element', 'Region', 'Function_Element_Final_weight']],
                                on=['#CHR', 'txID', 'Element', 'Region'], how='left')
    df_split_gene_sv['Function_Element_Final_weight'] = df_split_gene_sv['Function_Element_Final_weight'].fillna(1)
    # RE_elements
    df_regulatory_sv = df_judge_sv_regulatory.copy()
    df_regulatory_element_sv_weight['#CHR'] = df_regulatory_element_sv_weight['#CHR'].astype(str)
    df_regulatory_sv = pd.merge(df_regulatory_sv, df_regulatory_element_sv_weight[
        ['#CHR', 'geneID', 'RE_ID', 'Function_Element_Final_weight']],
                                on=['#CHR', 'geneID', 'RE_ID'], how='left')
    df_regulatory_sv['Function_Element_Final_weight'] = df_regulatory_sv['Function_Element_Final_weight'].fillna(1)
    # Get_Score
    df_sv = pd.concat([df_whole_gene_sv[['ID', '#CHR', 'SVTYPE', 'START', 'END', 'SVLEN']],
                       df_regulatory_sv[['ID', '#CHR', 'SVTYPE', 'START', 'END', 'SVLEN']]]
                      , axis=0)
    df_sv.drop_duplicates(subset=['ID', '#CHR', 'SVTYPE'], inplace=True)
    df_sv['#CHR'] = df_sv['#CHR'].astype(str)
    df_sv.reset_index(drop=True, inplace=True)
    utr_lengths_dict = scoring_func.calculate_utr_lengths(gene_annotation_dict)
    df_whole_gene_sv = utils.create_dict_by_column(df_whole_gene_sv, key_column='#CHR', index_column='ID')
    df_split_gene_sv = utils.create_dict_by_column(df_split_gene_sv, key_column='#CHR', index_column='ID')
    df_regulatory_sv = utils.create_dict_by_column(df_regulatory_sv, key_column='#CHR', index_column='ID')
    # 10.28修改
    tad_score = regulatory_element_weight.get('TAD', 0.3)
    if 'TAD' in regulatory_element_types:
        regulatory_element_types.remove('TAD')
        del regulatory_element_weight['TAD']
    # 12.18修改
    info_col = ['#CHR', 'ID', 'geneID', 'txID', 'relationship', 'Frameshit', 'GeneBody_AF', 'GeneBody_weight',
                'GeneBody_Score']
    for r_key in regulatory_element_types:
        info_col.append('Overlap_%s' % r_key)
        info_col.append('%s_AF' % r_key)
        info_col.append('%s_weight' % r_key)
        info_col.append('%s_Score' % r_key)
    info_col += ['LOEUF_bin', 'HI_bin', 'TI_bin', 'Score']
    sv_score_info_dict = {}
    for c_i in info_col:
        sv_score_info_dict[c_i] = []
    final_scores = []
    tra_id_index = {}
    sv_impact_tad = []
    sv_impact_tad_af = []
    fill_na_dict = {'relationship': 'No Overlap', 'GeneBody_Score': 0, 'GeneBody_AF': 0, 'GeneBody_weight': 1,
                    'Frameshit': 'NaN', 'txID': 'NaN'}
    regulatory_info_cols = ['geneID']
    for r_key in regulatory_element_types:
        regulatory_info_cols.append('Overlap_%s' % r_key)
        regulatory_info_cols.append('%s_AF' % r_key)
        regulatory_info_cols.append('%s_weight' % r_key)
        regulatory_info_cols.append('%s_Score' % r_key)
    for r_key in regulatory_element_types:
        fill_na_dict['Overlap_%s' % r_key] = 'No'
        fill_na_dict['%s_AF' % r_key] = '0'
        fill_na_dict['%s_weight' % r_key] = '1'
        fill_na_dict['%s_Score' % r_key] = '0'
    for index, row in tqdm(df_sv.iterrows()):
        sv_id = row['ID']
        chrom = row['#CHR']
        sv_type = row['SVTYPE']
        overlap_tad = 'No'
        try:
            df_whole_gene_sv_temp = df_whole_gene_sv[chrom].loc[[sv_id]]
            df_split_gene_sv_temp = df_split_gene_sv[chrom].loc[[sv_id]]
        except KeyError:
            df_whole_gene_sv_temp = pd.DataFrame()
            df_split_gene_sv_temp = pd.DataFrame()
        try:
            df_regulatory_sv_temp = df_regulatory_sv[chrom].loc[[sv_id]]
        except KeyError:
            df_regulatory_sv_temp = pd.DataFrame()
        if len(df_whole_gene_sv_temp) > 0:
            df_split_gene_sv_temp = scoring_func.give_split_gene_sv_scores(df_sv_gene_split=df_split_gene_sv_temp,
                                                                           sv_type=sv_type, common_af=common_af)
            df_split_gene_sv_utr = df_split_gene_sv_temp.loc[df_split_gene_sv_temp['Region'].str.contains('UTR'), :]
            if len(df_split_gene_sv_utr) != 0:
                df_split_gene_sv_temp = df_split_gene_sv_temp.loc[
                                        ~(df_split_gene_sv_temp['Region'].str.contains('UTR')), :]
                # give_utr_sv_scores(df_utr_sv, sv_type, utr_length_dict)
                df_sv_utr_score = scoring_func.give_utr_sv_scores(df_sv_split_gene_utr=df_split_gene_sv_utr,
                                                                  sv_type=sv_type,
                                                                  utr_length_dict=utr_lengths_dict)
                df_split_gene_sv_temp = df_split_gene_sv_temp.append(df_sv_utr_score, ignore_index=False)
            # df_gene_group = df_split_gene_sv_temp.groupby(by = ['geneID']).max()[['tx_AF','tx_Score']]
            # df_gene_group = df_gene_group.reset_index(drop=False)
            df_split_gene_sv_temp.sort_values(by=['GeneBody_Score'], inplace=True, ascending=False)
            df_gene_group = df_split_gene_sv_temp.drop_duplicates(subset=['geneID'], keep='first')[
                ['geneID', 'GeneBody_AF', 'GeneBody_weight', 'GeneBody_Score']]
            df_gene_group = pd.merge(df_gene_group,
                                     df_whole_gene_sv_temp[['geneID', 'txID', 'relationship', 'Frameshit']],
                                     on='geneID')
        else:
            df_gene_group = pd.DataFrame(
                columns=['geneID', 'txID', 'GeneBody_AF', 'GeneBody_weight', 'GeneBody_Score', 'relationship',
                         'Frameshit'])
        if len(df_regulatory_sv_temp) > 0:
            if 'TAD' in df_regulatory_sv_temp['RE_type'].values.tolist():
                overlap_tad = 'Yes'
                sv_impact_tad.append(sv_id)
                df_regulatory_sv_only_tad = df_regulatory_sv_temp.loc[df_regulatory_sv_temp['RE_type'] == 'TAD']
                overlap_tad_af = df_regulatory_sv_only_tad['AF'].max()
                try:
                    overlap_tad_weight_for_element_sv = df_regulatory_sv_only_tad['Function_Element_Final_weight'].max()
                except KeyError:
                    overlap_tad_weight_for_element_sv = 1
                    logging.error("Overlap_TAD error: Could not find Function_Element_Final_weight")
                sv_impact_tad_af.append(overlap_tad_af)
                df_regulatory_sv_no_tad = df_regulatory_sv_temp.loc[df_regulatory_sv_temp['RE_type'] != 'TAD']
                if len(df_regulatory_sv_no_tad) > 0:
                    df_regulatory_sv_scores = scoring_func.integrate_regulatory_sv_score(
                        df_regulatory_sv=df_regulatory_sv_no_tad,
                        regulatory_element_types=regulatory_element_types,
                        regulatory_element_weight=regulatory_element_weight, common_af=common_af)
                else:
                    df_regulatory_sv_scores = pd.DataFrame(columns=regulatory_info_cols)
            else:
                df_regulatory_sv_scores = scoring_func.integrate_regulatory_sv_score(
                    df_regulatory_sv=df_regulatory_sv_temp,
                    regulatory_element_types=regulatory_element_types,
                    regulatory_element_weight=regulatory_element_weight, common_af=common_af)
                overlap_tad = 'No'
        else:
            df_regulatory_sv_scores = pd.DataFrame(columns=regulatory_info_cols)
        df_gene_group = pd.merge(df_gene_group, df_regulatory_sv_scores, on=['geneID'], how='outer')

        df_gene_group.fillna(fill_na_dict, inplace=True)

        #print('#debug0')
        #print(df_gene_group)
        if sv_type == 'DUP' and 'SV Partial overlap Tx' in df_gene_group['relationship'].values.tolist():
            df_gene_group = scoring_func.adjust_scores_for_sv(df_data=df_gene_group, sv_type=sv_type,
                                                              gene_annotation_dict=gene_annotation_dict,
                                                              chrom=chrom, sv_breakend1=row['START'],
                                                              sv_breakend2=row['END'], overlap_tad=overlap_tad)
        elif sv_type == 'INV' and 'SV contains all Tx' in df_gene_group['relationship'].values.tolist():
            df_gene_group = scoring_func.adjust_scores_for_sv(df_data=df_gene_group, sv_type=sv_type,
                                                              gene_annotation_dict=gene_annotation_dict,
                                                              chrom=chrom, sv_breakend1=row['START'],
                                                              sv_breakend2=row['END'], overlap_tad=overlap_tad)
        #print('#debug1')
        #print(df_gene_group)
        df_gene_group = scoring_func.give_whole_gene_conservation_scores(df_data=df_gene_group, sv_type=sv_type,
                                                                         regulatory_element_types=regulatory_element_types,
                                                                         gene_hi_lof=gene_hi_lof)
        #print('#debug2')
        #print(df_gene_group)
        # TAD权重
        if mode == 'Max':
            if overlap_tad == 'Yes':
                sum_s = scoring_func.max_score(df_data=df_gene_group) + tad_score * scoring_func.weight_by_af(
                    overlap_tad_af, common_af=common_af) * overlap_tad_weight_for_element_sv  # 10.28修改
            else:
                sum_s = scoring_func.max_score(df_data=df_gene_group)
            if sv_type in ['TRA', 'BND'] and sv_id not in tra_id_index:
                tra_id_index[sv_id] = index
                final_scores.append(sum_s)
            elif sv_type in ['TRA', 'BND'] and id in tra_id_index:
                sum_s = max(sum_s, final_scores[tra_id_index[sv_id]])
                final_scores.append(sum_s)
                final_scores[tra_id_index[sv_id]] = sum_s
            else:
                final_scores.append(sum_s)
        elif mode == 'Sum':
            if overlap_tad == 'Yes':
                sum_s = scoring_func.sum_score(df_data=df_gene_group) + tad_score * scoring_func.weight_by_af(
                    overlap_tad_af, common_af=common_af) * overlap_tad_weight_for_element_sv  # 10.28修改
            else:
                sum_s = scoring_func.sum_score(df_data=df_gene_group)
            if sv_type in ['TRA', 'BND'] and sv_id not in tra_id_index:
                tra_id_index[sv_id] = index
                final_scores.append(sum_s)
            elif sv_type in ['TRA', 'BND'] and sv_id in tra_id_index:
                sum_s = sum_s + final_scores[tra_id_index[sv_id]]
                final_scores.append(sum_s)
                final_scores[tra_id_index[sv_id]] = sum_s
            else:
                final_scores.append(sum_s)
        df_gene_group['#CHR'] = [chrom] * len(df_gene_group)
        df_gene_group['ID'] = [sv_id] * len(df_gene_group)
        for c_i in info_col:
            sv_score_info_dict[c_i] = sv_score_info_dict[c_i] + df_gene_group[c_i].values.tolist()
    df_sv_score_info = pd.DataFrame(sv_score_info_dict)
    sv_impact_tad = set(sv_impact_tad)
    df_sv_score_info['Overlap_TAD'] = df_sv_score_info['ID'].apply(
        lambda x: 'Yes' if x in sv_impact_tad else 'No')
    df_sv['Final_Score'] = final_scores
    df_sv_scores = df_sv.drop_duplicates(subset='ID', keep='first')
    df_sv_scores.reset_index(drop=True, inplace=True)

    # 10.28修改：给Intergenic(No Function)分值
    intergenic_score = 0.1
    df_sv_scores, df_sv_score_info = scoring_func.give_intergenic_sv_score(df_sv_scores, df_sv_score_info,
                                                                           df_judge_whole_sv_af, intergenic_score,
                                                                           common_af=common_af)
    # 12.18修改
    df_sv_score_info.drop_duplicates(subset=['#CHR', 'ID', 'geneID'], keep='first', inplace=True, ignore_index=True)
    return df_sv_score_info, df_sv_scores


# %%
def run_scoring(judge_sv_whole_gene_file, judge_sv_split_gene_file, annotation_tuple, judge_sv_regulatory_file,
                judge_whole_sv_af, mode, gene_hi_lof_file, output_dir, gene_element_sv_weight_files,
                regulatory_element_sv_weight_files, weight_source_dataset, common_af,
                prefix=''):
    """
    Run the scoring for structural variants.
    Args:
    judge_sv_whole_gene_file (str): Path to the whole gene structural variant file.
    judge_sv_split_gene_file (str): Path to the split gene structural variant file.
    annotation_tuple (tuple): tuple of annotation data.
    judge_sv_regulatory_file (str): Path to the regulatory element file.
    judge_whole_sv_af (str): Path to the whole structural variant allele frequency file.
    mode (str): Mode of scoring.
    gene_hi_lof_file (str): Path to the gene HI, TI, LOF file.
    output_dir (str): Path to the output directory.
    gene_element_sv_weight_files (dict): Dictionary of gene element structural variant weight path.
    regulatory_element_sv_weight_files (dict): Dictionary of regulatory element structural variant weight path.
    weight_source_dataset (str): Source dataset for the weight.
    common_af (float): Common AF Frequency.
    prefix (str, optional): Prefix for the scoring files.
    Returns:
    tuple: Paths to the scoring final file and annotation info file.
    """
    gene_annotation_dict, _, _, _ = annotation_tuple
    regulatory_element_types = ['Promoter', 'Enhancer', 'DNA_Accessible', 'TAD']
    regulatory_element_weight = {'Promoter': 0.5, 'Enhancer': 0.4, 'DNA_Accessible': 0.2, 'TAD': 0.3}
    df_judge_sv_whole_gene = pd.read_csv(judge_sv_whole_gene_file, sep='\t')
    df_judge_sv_split_gene = pd.read_csv(judge_sv_split_gene_file, sep='\t')
    df_judge_sv_regulatory = pd.read_csv(judge_sv_regulatory_file, sep='\t')
    df_judge_whole_sv_af = pd.read_csv(judge_whole_sv_af, sep='\t')
    gene_hi_lof = pd.read_csv(gene_hi_lof_file, sep='\t')
    gene_hi_lof.set_index('geneID', inplace=True)
    gene_element_sv_weight_dict, regulatory_element_sv_weight_dict = {}, {}
    for data in gene_element_sv_weight_files:
        gene_element_sv_weight_dict[data] = pd.read_csv(gene_element_sv_weight_files[data])
    for data in regulatory_element_sv_weight_files:
        regulatory_element_sv_weight_dict[data] = pd.read_csv(regulatory_element_sv_weight_files[data])
    df_out_sv_score_info, df_out_sv_scores = give_score(df_judge_sv_whole_gene=df_judge_sv_whole_gene,
                                                        df_judge_sv_split_gene=df_judge_sv_split_gene,
                                                        df_judge_sv_regulatory=df_judge_sv_regulatory,
                                                        df_judge_whole_sv_af=df_judge_whole_sv_af,
                                                        gene_hi_lof=gene_hi_lof,
                                                        gene_element_sv_weight_dict=gene_element_sv_weight_dict,
                                                        regulatory_element_sv_weight_dict=regulatory_element_sv_weight_dict,
                                                        gene_annotation_dict=gene_annotation_dict,
                                                        regulatory_element_types=regulatory_element_types,
                                                        regulatory_element_weight=regulatory_element_weight,
                                                        element_sv_weight_source=weight_source_dataset,
                                                        common_af=common_af,
                                                        mode=mode)
    if prefix == '':
        prefix = 'SVJudge'
    annotation_info_file = os.path.join(output_dir, prefix + '.SVJudge.Scores.detail.tsv')
    scoring_final_file = os.path.join(output_dir, prefix + '.SVJudge.Scores.' + mode.lower() + '.summary.tsv')
    df_out_sv_score_info.to_csv(annotation_info_file, sep='\t', index=False)
    df_out_sv_scores.to_csv(scoring_final_file, sep='\t', index=False)
    return scoring_final_file, annotation_info_file
