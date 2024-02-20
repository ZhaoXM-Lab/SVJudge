import logging
import os
import re

import pandas as pd
import statsmodels.api as sm

import src.annot_func as annot_func


def sv_counts_by_gene_element(gene_split_annot_dict, ref_gene_canonical_annot_dict, common_af, weight_source_dataset):
    """
    Calculate Structural Variations (SV) counts based on gene elements.

    :param gene_split_annot_dict: Dictionary containing gene split annotation data.
    :param ref_gene_canonical_annot_dict: Dictionary containing reference gene canonical annotations.
    :param common_af: Threshold value to identify common structural variations.
    :param weight_source_dataset: List of public SV sets to calculate the weights
    :return: Dictionary with counts of structural variations for each gene element.
    """
    # Extract gene element details from the reference gene canonical annotation dictionary
    ref_gene_canonical_split = annot_func.extract_gene_element_details(ref_gene_canonical_annot_dict)

    # Initialize an empty dictionary to store the SV counts for each gene element
    public_gene_split_sv_counts = {}

    # Iterate through the gene split annotation data
    for dataset_name, file_path in gene_split_annot_dict.items():
        # Read the data from the file into a pandas DataFrame
        if dataset_name not in weight_source_dataset:
            continue
        df_temp = pd.read_csv(file_path, sep='\t')

        # Calculate the total count and common count of SVs based on specified threshold
        total_count = df_temp.groupby(['#CHR', 'geneID', 'txID', 'Element', 'Region'])[
            'AF'].count().reset_index().rename(columns={'AF': 'total_count'})
        common_count = \
            df_temp[df_temp['AF'] >= common_af].groupby(['#CHR', 'geneID', 'txID', 'Element', 'Region'])[
                'AF'].count().reset_index().rename(columns={'AF': 'common_count'})

        # Merge the total count and common count DataFrames
        count_merge = pd.merge(total_count, common_count, how='left',
                               on=['#CHR', 'geneID', 'txID', 'Element', 'Region']).fillna(0)

        # Merge with the reference gene canonical split details
        count_merge = pd.merge(count_merge, ref_gene_canonical_split[['txID', 'Element', 'Region', 'length']],
                               on=['txID', 'Element', 'Region'], how='left')

        # Calculate densities
        count_merge['total_density'] = count_merge['total_count'] / count_merge['length']
        count_merge['common_density'] = count_merge['common_count'] / count_merge['length']

        # Store the SV counts for the current dataset in the dictionary
        public_gene_split_sv_counts[dataset_name] = count_merge
    return public_gene_split_sv_counts


# %%
def count_sv_by_regulatory_element(re_split_annot_dict, re_ref_bed_files, re_keys, common_af,weight_source_dataset):
    """
    Calculate Structural Variations (SV) counts based on regulatory elements.

    :param re_split_annot_dict: Dictionary containing RE split annotation data.
    :param re_ref_bed_files: List of file paths for RE reference bed files.
    :param re_keys: List of keys corresponding to the RE reference bed files.
    :param common_af: Threshold value to identify common structural variations.
    :param weight_source_dataset: List of public SV sets to calculate the weights
    :return: Dictionary with counts of structural variations for each regulatory element.
    """
    public_re_split_sv_counts = {}
    df_re_annotation = pd.DataFrame()
    columns = ['#CHR', 'RE_Start', 'RE_end', 'Map_gene', 'RE_ID']

    # Aggregate RE reference data from multiple files
    for key, file in zip(re_keys, re_ref_bed_files):
        df_re_annotation = df_re_annotation.append(pd.read_csv(file, sep='\t', header=None, names=columns))

    # Calculate length of each RE
    df_re_annotation['length'] = df_re_annotation.apply(lambda x: abs(int(x['RE_Start']) - int(x['RE_end'])), axis=1)

    # Process each annotation dataset
    for dataset_name, file_path in re_split_annot_dict.items():
        if dataset_name not in weight_source_dataset:
            continue
        df_temp = pd.read_csv(file_path, sep='\t')

        # Count total and common SVs
        total_count = df_temp.groupby(['#CHR', 'geneID', 'RE_ID'])['AF'].count().reset_index().rename(
            columns={'AF': 'total_count'})
        common_count = df_temp[df_temp['AF'] >= common_af].groupby(['#CHR', 'geneID', 'RE_ID'])[
            'AF'].count().reset_index().rename(columns={'AF': 'common_count'})

        # Merge counts and calculate densities
        count_merge = pd.merge(total_count, common_count, how='left', on=['#CHR', 'geneID', 'RE_ID']).fillna(0)
        count_merge = pd.merge(count_merge, df_re_annotation[['RE_ID', 'length']], on=['RE_ID'], how='left')
        count_merge['total_density'] = count_merge['total_count'] / count_merge['length']
        count_merge['common_density'] = count_merge['common_count'] / count_merge['length']
        public_re_split_sv_counts[dataset_name] = count_merge

    return public_re_split_sv_counts


# %%
def split_class_by_resource(row, resource=''):
    """
    Generate a label based on the specified resource type.

    :param row: A dictionary or row containing the data to process.
    :param resource: The type of resource ('gene' or 'RE').
    :return: A string label derived from the specified fields in the row.
    """
    # Initialize the label_class
    label_class = ''

    # Check the resource type
    if resource == 'gene':
        # Process the element and region fields
        element = re.sub(r'\d+', '', row['Element'])
        region = re.sub(r'\d+', '', row['Region']).replace('\'', '')
        # Generate label for gene resource
        label_class = f"{element}_{region}"
    elif resource == 'RE':
        # Generate label for RE resource
        label_class = re.sub(r'\d+', '', row['RE_ID'])

    # Return the generated label
    return label_class


# %%
def calculate_weight_from_studentized_residual(studentized_residual):
    """
    Calculate a weight based on the studentized residual.

    :param studentized_residual: The studentized residual value.
    :return: The calculated weight, rounded to 2 decimal places.
    """
    weight = min(2, 2 ** abs(studentized_residual / 3))
    return round(weight if studentized_residual <= 0 else 1 / weight, 2)


# %%
def calculate_weights(df):
    """
    Calculate weights using Ordinary Least Squares (OLS) model.

    :param df: DataFrame with density data.
    :return: DataFrame with calculated weights.
    """
    # Fit the OLS model
    model = sm.OLS(df['common_density'].values,
                   df['total_density'].values.reshape(-1, 1)).fit()

    # Add predicted values to the DataFrame
    df['predict'] = model.predict(df['total_density'])

    # Add residuals to the DataFrame
    df['residuals'] = model.resid

    # Calculate influence and add studentized residuals to the DataFrame
    influence = model.get_influence()
    df['resid_studentized'] = influence.resid_studentized_external

    # Calculate weights from studentized residuals and add to the DataFrame
    df['weight'] = df['resid_studentized'].apply(calculate_weight_from_studentized_residual)

    return df


# %%
def process_data_split(df, resource_type):
    # Add a new column 'Split_Class' based on the resource type
    df['Split_Class'] = df.apply(lambda x: split_class_by_resource(x, resource_type), axis=1)

    # Create an empty DataFrame to store weights
    weight_df = pd.DataFrame()

    # Get the unique split classes from the DataFrame
    split_classes = df['Split_Class'].value_counts().index.tolist()

    # Iterate through the split classes and calculate weights for each
    for sc in tqdm(split_classes):
        # Filter the DataFrame based on the split class
        df_temp = df[df['Split_Class'] == sc].reset_index(drop=True)
        # Calculate weights for the filtered DataFrame and append to weight_df
        weight_df = weight_df.append(calculate_weights(df_temp), ignore_index=True)

    return weight_df


def process_data_all(df_gene, df_re):
    df_gene['general_ID'] = ['G_%s' % i for i in df_gene.index]
    df_re['general_ID'] = ['R_%s' % i for i in df_re.index]
    combined_df = df_gene[['general_ID', 'total_density', 'common_density']].append(
        df_re[['general_ID', 'total_density', 'common_density']]).reset_index(drop=True)
    gene_weight_df = pd.merge(df_gene, combined_df, on=['general_ID', 'total_density', 'common_density'], how='left')
    re_weight_df = pd.merge(df_re, combined_df, on=['general_ID', 'total_density', 'common_density'], how='left')
    del gene_weight_df['general_ID']
    del re_weight_df['general_ID']
    return gene_weight_df, re_weight_df


from tqdm import tqdm


def calculate_weight_for_elements(gene_sv_counts, re_sv_counts, output_dir, used_data='', method='split'):
    """
     Process the data for a specific resource type (gene or RE).

     :param df: DataFrame of either gene or RE data.
     :param resource_type: 'gene' or 'RE'
     :param method: 'split' or 'All'
     :return: DataFrame after processing.
     """
    gene_weights, re_weights = {}, {}

    used_data = used_data or list(gene_sv_counts.keys())
    for data in tqdm(used_data):
        if method == 'split':
            gene_weights[data] = os.path.join(output_dir, '%s.sv_gene_elements.conservation_weights.temp.csv' % data)
            df_gene_weights = process_data_split(gene_sv_counts[data], 'gene')
            df_gene_weights.to_csv(gene_weights[data], index=False)
            re_weights[data] = os.path.join(output_dir, '%s.sv_re_elements.conservation_weights.temp.csv' % data)
            df_re_weights = process_data_split(re_sv_counts[data], 'RE')
            df_re_weights.to_csv(re_weights[data], index=False)
        elif method == 'All':
            gene_weights[data], re_weights[data] = process_data_all(gene_sv_counts[data], re_sv_counts[data])
        else:
            logging.error(f"Invalid method in calculating conserveration weights for gene and RE elements : {method}")
    return gene_weights, re_weights


def run_calculate_conservation_weight(sv_gene_split_annotation, gene_annotation_dict, sv_regulation_annot_dict,
                                      regulatory_reference_beds, weight_source_dataset,
                                      regulatory_element_types, common_af, output_dir, ):
    """
    Calculate the weights for structural variations in gene and regulatory elements.

    :param sv_gene_split_annotation: Dictionary with gene split annotation data.
    :param gene_annotation_dict: Dictionary with canonical gene annotations.
    :param sv_regulation_annot_dict: Dictionary with regulatory elements split annotation data.
    :param regulatory_reference_beds: List of bed files for regulatory elements.
    :param regulatory_element_types: List of keys for regulatory elements reference.
    :param common_af: Common threshold for structural variation analysis.
    :param output_dir: Output file path
    :param weight_source_dataset: List of public SV sets to calculate the weights
    :return: Tuple of two dictionaries containing weights for gene and RE elements respectively.
    """
    logging.info("calculate sv counts in every gene element")
    gene_sv_counts = sv_counts_by_gene_element(sv_gene_split_annotation, gene_annotation_dict, common_af,
                                               weight_source_dataset)
    logging.info("calculate sv counts in every regulatory element")
    re_sv_counts = count_sv_by_regulatory_element(sv_regulation_annot_dict, regulatory_reference_beds,
                                                  regulatory_element_types,
                                                  common_af, weight_source_dataset)

    gene_sv_weights, re_sv_weights = calculate_weight_for_elements(
        gene_sv_counts=gene_sv_counts,
        re_sv_counts=re_sv_counts,
        used_data='',
        method='split',
        output_dir=output_dir
    )
    return gene_sv_weights, re_sv_weights
