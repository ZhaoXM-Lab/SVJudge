#!/usr/bin/env python3
import argparse
import logging
import os
import sys
from time import strftime, localtime

import src.utils as utils
from src.element_sv_conservation_func import run_calculate_conservation_weight
from src.run_judge_sv_annot import run_annotation_only_judge_sv
from src.run_judge_sv_scoring import run_scoring
from src.run_public_sv_annot import run_annotation_only_public_sv
from src.run_whole_sv_annot import run_whole_sv_annotate
from src.version import __version__
import shutil
import warnings
# 屏蔽所有警告信息
warnings.filterwarnings('ignore')

def parse_arguments(arguments=sys.argv[1:]):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="""SVJudge {0} \n \nShort Usage: SVJudge [parameters] -o <output path> -i <input vcf path> -r <reference> -t <task>""".format(
                                         __version__))

    required_params = parser.add_argument_group("Input/Output parameters")
    required_params.add_argument('-o', dest="output_path", type=os.path.abspath, required=True,
                                 help='Absolute path to output ')
    required_params.add_argument('-i', dest='input_vcf', type=os.path.abspath, required=True,
                                 help='Absolute path to vcf file')
    required_params.add_argument('-r', dest="reference", choices=['GRCh37', 'GRCh38'], type=str,
                                 required=True
                                 , help='Absolute path to your reference genome (GRCh37 or  GRCh38)')

    optional_params = parser.add_argument_group("Optional parameters")
    optional_params.add_argument('-op', dest="overlap_threshold", type=int, default=0.5,
                                 help='SV overlap threshold(default: %(default)s)')
    optional_params.add_argument('-dis', dest="distance_threshold", type=int, default=1000,
                                 help='SV distance threshold (default: %(default)s)')
    optional_params.add_argument('-common_af', dest="common_af", type=float, default=0.01,
                                 help='AF threshold of common SVs(default: %(default)s)')
    optional_params.add_argument('-prefix', dest="prefix", type=str, default='SVJudge',
                                 help='The prefix of result file (default: %(default)s)')
    optional_params.add_argument('-annot_sv_dir', dest="prefix", type=os.path.abspath, default='',
                                 help='Only used in task score (default: %(default)s)')
    optional_params.add_argument('-use_sv_set', dest='use_sv_set',
                                 default='gnomAD, 1000PG_phase3, DGV, IMH, PGG_LRS_SV, PGG_SRS_SV',
                                 help='Public SV sets will use in task, default: %(default)s'
                                      '\nMultiple inputs separated by commas. For example: -use_sv_set 1000PG_phase3,gnomAD')
    optional_params.add_argument('-add_sv_set_name', dest='add_sv_set_name', action="append",
                                 help='The name of new public SV set added in task (default: %(default)s)'
                                      '\nThis parameter can be specified multiple times to add multiple SV sets')
    optional_params.add_argument('-add_sv_set_dir', dest='add_sv_set_dir', action="append",
                                 help='The name of new public SV set added in task (default: %(default)s)'
                                      '\nThis parameter can be specified multiple times to add multiple SV sets')
    optional_params.add_argument('-integrate_mode', dest='integrate_mode', default='Max', choices=['Max', 'Sum'],
                                 help='How to integrate the score of SV impacted multiple function regions (Max or Sum), default: %(default)s)')
    optional_params.add_argument('-sv_set_for_weight', dest='weight_source_dataset', default='PGG_LRS_SV',
                                 help='SV sets of element weight data, default: %(default)s)'
                                      '\nMultiple inputs separated by commas. For example: -sv_set_for_weight PGG_LRS_SV, gnomAD')

    options = parser.parse_args(arguments)
    return options


if __name__ == '__main__':

    options = parse_arguments()
    work_dir = options.output_path
    temp_dir = os.path.join(work_dir, 'temp')
    if not os.path.exists(work_dir):
        print(f'Create the output directory {work_dir}')
        os.mkdir(work_dir)

    log_format = logging.Formatter("%(asctime)s [%(levelname)-7.7s]  %(message)s")
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)

    fileHandler = logging.FileHandler("{0}/SVJudge_{1}.log".format(work_dir, strftime("%y%m%d_%H%M%S", localtime())),
                                      mode="w")
    fileHandler.setFormatter(log_format)

    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(log_format)
    root_logger.addHandler(fileHandler)
    logging.info('******************** Start SVJudge, version {0} ********************'.format(__version__))
    logging.info("CMD: {0}".format(" ".join(sys.argv)))
    logging.info("WORKDIR DIR: {0}".format(os.path.abspath(work_dir)))
    # logging.info("TASK: {0}".format(options.task))
    logging.info("INPUT VCF: {0}".format(os.path.abspath(options.input_vcf)))
    logging.info("REFERENCE GENOME: {0}".format(options.reference))
    script_path = os.path.abspath(__file__)
    data_dir = os.path.join(os.path.dirname(script_path), 'data')
    if options.reference.upper() == 'GRCh37'.upper():
        data_dir = os.path.join(data_dir, 'GRCh37')
    elif options.reference.upper() == 'GRCh38'.upper():
        data_dir = os.path.join(data_dir, 'GRCh38')
    # annotation files
    annotation_files_dict = {'Gene': os.path.join(data_dir, 'gene_annotation/gene_structure/refGene.txt')}
    regulation_dir = os.path.join(data_dir, 'regulatory_element')
    for element in os.listdir(regulation_dir):
        if element not in ['Promoter','Enhancer','DNA_Accessible','TAD']:
            logging.error("Cannot identified regulatory elements %s" %element)
        subdir = os.path.join(regulation_dir,element)
        annot_file = os.listdir(subdir)[0]
        annotation_files_dict[element] = os.path.join(subdir,annot_file)

    # public sv sets
    public_sv_set_paths = {'gnomAD': os.path.join(data_dir, 'public_sv_set/gnomAD-sv.preprocessed.bed'),
                           '1000PG_phase3': os.path.join(data_dir, 'public_sv_set/1000GP_phase3-sv.preprocessed.bed'),
                           'PGG_LRS_SV': os.path.join(data_dir, 'public_sv_set/PGG_SV.LRS.preprocessed.bed'),
                           'PGG_SRS_SV': os.path.join(data_dir, 'public_sv_set/PGG_SV.SRS.preprocessed.bed'),
                           'IMH': os.path.join(data_dir, 'public_sv_set/IMH-sv.preprocessed.bed'),
                           'DGV': os.path.join(data_dir, 'public_sv_set/DGV-sv.preprocessed.bed')}
    gene_hi_lof_file = os.path.join(data_dir, 'gene_annotation/gene_intolerance_annotations/GENE_HI_LOF.merge.tsv')
    use_public_sv_sets = options.use_sv_set.split(',')
    use_public_sv_sets = [i.strip() for i in use_public_sv_sets]
    if options.add_sv_set_name is not None and options.add_sv_set_dir is not None:
        if len(options.add_sv_set_name) != len(options.add_sv_set_dir):
            logging.error('The length of add_sv_set_name and add_sv_set_dir should be the same')
            sys.exit()
        else:
            for add_name, add_dir in zip(options.add_sv_set_name, options.add_sv_set_dir):
                public_sv_set_paths[add_name] = add_dir
                if not os.path.exists(add_dir):
                    logging.error('The directory for {0} not exist'.format(add_dir))
                    sys.exit()
        use_public_sv_sets = list(set(use_public_sv_sets + options.add_sv_set_name))
    elif options.add_sv_set_name is not None or options.add_sv_set_dir is not None:
        logging.error('The length of add_sv_set_name and add_sv_set_dir should be the same')
        sys.exit()
    weight_source_dataset = [i.strip() for i in options.weight_source_dataset.split(',')]
    if not set(weight_source_dataset).issubset(set(use_public_sv_sets)):
        logging.info('weight_source_dataset:%s' % options.weight_source_dataset)
        logging.info('use_public_sv_sets:%s' % options.use_sv_set)
        logging.error(
            'The public SV set specified in \"-sv_set_for_weight\" should also be included in \"-use_sv_set.\"')
        sys.exit()
    for use_set in use_public_sv_sets:
        if use_set not in public_sv_set_paths:
            logging.error('The public SV set {0} was not found'.format(use_set))
            sys.exit()
    logging.info('"Importing input VCF...')
    temp_dir = os.path.join(work_dir, 'temp')
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)
    input_bed = '.'.join(options.input_vcf.split('/')[-1].split('.')[:-1]) + 'SVJudge.input.bed'
    input_bed_path = os.path.join(temp_dir, input_bed)
    judge_sv_file = utils.vcf2bed(options.input_vcf, input_bed_path, include_TRA=False)

    logging.info('Loading and processing annotation files...')
    gene_annotation_dict, gene_reference_bed, regulatory_reference_beds, regulatory_element_types = utils.load_annotation_files(
        annotation_files_dict=annotation_files_dict, output_dir=temp_dir)
    annotation_tuple = (gene_annotation_dict, gene_reference_bed, regulatory_reference_beds, regulatory_element_types)
    split_gene_public_sv_annot_merge_sets_file, split_gene_sv_annot_dict, regulation_sv_annot_merge_sets_file, regulation_sv_annot_dict = run_annotation_only_public_sv(
        annotation_tuple=annotation_tuple,
        input_sv_bed_paths=public_sv_set_paths,
        selected_datasets=use_public_sv_sets,
        overlap_threshold=options.overlap_threshold,
        distance_threshold=options.distance_threshold,
        output_dir=temp_dir)
    #debug
    """
    split_gene_sv_annot_dict = {}
    regulation_sv_annot_dict = {}
    gene_element_sv_conservation = {}
    regulation_element_sv_conservation = {}
    for file in os.listdir(temp_dir):
        if file.endswith('.gene_annot.split.tsv'):
            data = file.replace('.gene_annot.split.tsv','')
            if data not in use_public_sv_sets:
                continue
            split_gene_sv_annot_dict[data] = os.path.join(temp_dir,file)
        elif file.endswith('.RE_annot.tsv'):
            data = file.replace('.RE_annot.tsv', '')
            if data not in use_public_sv_sets:
                continue
            regulation_sv_annot_dict[data] = os.path.join(temp_dir, file)
        elif file.endswith('.sv_gene_elements.conservation_weights.temp.csv'):
            data = file.replace('.sv_gene_elements.conservation_weights.temp.csv', '')
            if data not in weight_source_dataset:
                continue
            gene_element_sv_conservation[data] = os.path.join(temp_dir,file)
        elif file.endswith('.sv_re_elements.conservation_weights.temp.csv'):
            data = file.replace('.sv_re_elements.conservation_weights.temp.csv', '')
            if data not in weight_source_dataset:
                continue
            regulation_element_sv_conservation[data] = os.path.join(temp_dir, file)
    split_gene_public_sv_annot_merge_sets_file = os.path.join(temp_dir, 'Public_SV_sets.Merge.gene_annot.split.tsv')
    regulation_sv_annot_merge_sets_file = os.path.join(temp_dir, 'Public_SV_sets.Merge.RE_annot.tsv')
    """

    gene_element_sv_conservation, regulation_element_sv_conservation = run_calculate_conservation_weight(
        sv_gene_split_annotation=split_gene_sv_annot_dict,
        gene_annotation_dict=gene_annotation_dict,
        sv_regulation_annot_dict=regulation_sv_annot_dict,
        regulatory_reference_beds=regulatory_reference_beds,
        regulatory_element_types=regulatory_element_types,
        common_af=options.common_af,
        output_dir=temp_dir,
        weight_source_dataset=weight_source_dataset)

    logging.info('Annotating input SVs...')
    judge_sv_split_gene_annot_file, judge_sv_whole_gene_annot_file, judge_sv_regulatory_annot_file = run_annotation_only_judge_sv(
        annotation_tuple=annotation_tuple,
        input_sv_bed_path=judge_sv_file,
        split_gene_public_sv_af_file=split_gene_public_sv_annot_merge_sets_file,
        regulatory_public_sv_af_file=regulation_sv_annot_merge_sets_file,
        overlap_threshold=options.overlap_threshold,
        distance_threshold=options.distance_threshold,
        output_dir=work_dir,
        output_prefix=options.prefix)
    whole_judge_sv_annot_af = run_whole_sv_annotate(judge_sv_file=judge_sv_file,
                                                    public_sv_sets=public_sv_set_paths,
                                                    selected_datasets=use_public_sv_sets,
                                                    overlap_threshold=options.overlap_threshold,
                                                    distance_threshold=options.distance_threshold,
                                                    output_dir=temp_dir,
                                                    out_prefix=options.prefix)
    print("whole_judge_sv_annot_af: %s" % whole_judge_sv_annot_af)
    logging.info('Annotating input SVs...')

    run_scoring(judge_sv_whole_gene_file=judge_sv_whole_gene_annot_file,
                judge_sv_split_gene_file=judge_sv_split_gene_annot_file, annotation_tuple=annotation_tuple,
                judge_sv_regulatory_file=judge_sv_regulatory_annot_file,
                judge_whole_sv_af=whole_judge_sv_annot_af,
                mode=options.integrate_mode,
                gene_hi_lof_file=gene_hi_lof_file,
                output_dir=work_dir,
                gene_element_sv_weight_files=gene_element_sv_conservation,
                regulatory_element_sv_weight_files=regulation_element_sv_conservation,
                weight_source_dataset=weight_source_dataset,
                common_af=options.common_af,
                prefix=options.prefix)
    shutil.rmtree(temp_dir)
    suffixes_to_keep = ['max.summary.tsv','.annotated.tsv','.log','.detail.tsv']
    for filename in os.listdir(work_dir):
        file_path = os.path.join(work_dir, filename)
        if os.path.isfile(file_path) and not any(filename.endswith(suffix) for suffix in suffixes_to_keep):
            os.remove(file_path)
