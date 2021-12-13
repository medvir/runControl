# -*- coding: utf-8 -*-

"""Estimate a quality score for the run control file with respect to
the training run control files.

    filein -- csv file in the following format:
        freq| mut| pos

    mut: the nc in the given position
    freq: the frequency of the nc at the given position
    pos: the absolute position in HXB2 reference

Choose the variable sites with high frequency [0.2, 0.8], calculate score as the ratio of
number of sites present in the training set and in the sample. For example:
variable sites in the sample 10
training sites 20
score_sample = 10/20
"""
import logging
import sys
import os
from functools import reduce
import pandas as pd
import glob

MEAN_MIN_THRESHOLD = 0.2
MEAN_MAX_THRESHOLD = 0.8
MAX_STD = 0.25
MIN_STD = 0.00001


def calc_mean_std_population(df_merged, discard_file=None):
    """Calculate mean and STD of the traning datasets
    """
    A_freq_cols = [col for col in df_merged.columns if 'A' in col]
    C_freq_cols = [col for col in df_merged.columns if 'C' in col]
    G_freq_cols = [col for col in df_merged.columns if 'G' in col]
    T_freq_cols = [col for col in df_merged.columns if 'T' in col]

    # Make sure at each position, there are at least THRESHOLD_NT_PRESENT identified nc
    THRESHOLD_NT_PRESENT = 4
    df_merged.loc[:, 'A_mean'] = df_merged.loc[df_merged[A_freq_cols].count(axis=1) >
                                               THRESHOLD_NT_PRESENT, A_freq_cols].mean(axis=1)
    df_merged.loc[:, 'A_std'] = df_merged.loc[df_merged[A_freq_cols].count(axis=1) >
                                              THRESHOLD_NT_PRESENT, A_freq_cols].std(axis=1)
    df_merged.loc[:, 'C_mean'] = df_merged.loc[df_merged[C_freq_cols].count(axis=1) >
                                               THRESHOLD_NT_PRESENT, C_freq_cols].mean(axis=1)
    df_merged.loc[:, 'C_std'] = df_merged.loc[df_merged[C_freq_cols].count(axis=1) >
                                              THRESHOLD_NT_PRESENT, C_freq_cols].std(axis=1)
    df_merged.loc[:, 'G_mean'] = df_merged.loc[df_merged[G_freq_cols].count(axis=1) >
                                               THRESHOLD_NT_PRESENT, G_freq_cols].mean(axis=1)
    df_merged.loc[:, 'G_std'] = df_merged.loc[df_merged[G_freq_cols].count(axis=1) >
                                              THRESHOLD_NT_PRESENT, G_freq_cols].std(axis=1)
    df_merged.loc[:, 'T_mean'] = df_merged.loc[df_merged[T_freq_cols].count(axis=1) >
                                               THRESHOLD_NT_PRESENT, T_freq_cols].mean(axis=1)
    df_merged.loc[:, 'T_std'] = df_merged.loc[df_merged[T_freq_cols].count(axis=1) >
                                              THRESHOLD_NT_PRESENT, T_freq_cols].std(axis=1)

    df_mean_std = df_merged.loc[:, 'A_mean':]
    nt_m_std = df_mean_std.copy()

    col_mean_name = [col for col in nt_m_std if 'mean' in col]
    col_std_name = [col for col in nt_m_std if 'std' in col]

    # Select positions with average variation between 0.2 and 0.8
    nt_m_std_no_ones = nt_m_std[~((nt_m_std[col_mean_name[0]] == 1) |
                                  (nt_m_std[col_mean_name[1]] == 1) |
                                  (nt_m_std[col_mean_name[2]] == 1) |
                                  (nt_m_std[col_mean_name[3]] == 1))]

    idx_min = ~((nt_m_std_no_ones[col_mean_name[0]] < MEAN_MIN_THRESHOLD) |
                (nt_m_std_no_ones[col_mean_name[1]] < MEAN_MIN_THRESHOLD) |
                (nt_m_std_no_ones[col_mean_name[2]] < MEAN_MIN_THRESHOLD) |
                (nt_m_std_no_ones[col_mean_name[3]] < MEAN_MIN_THRESHOLD))

    nt_m_std_threshold_min = nt_m_std_no_ones[idx_min]

    idx_max = ~((nt_m_std_no_ones[col_mean_name[0]] > MEAN_MAX_THRESHOLD) |
                (nt_m_std_no_ones[col_mean_name[1]] > MEAN_MAX_THRESHOLD) |
                (nt_m_std_no_ones[col_mean_name[2]] > MEAN_MAX_THRESHOLD) |
                (nt_m_std_no_ones[col_mean_name[3]] > MEAN_MAX_THRESHOLD))

    nt_m_std_filtered_mean = nt_m_std_threshold_min.loc[idx_max]

    # Filter out positions with higher than 0.25 std
    nt_m_std_filtered = nt_m_std_filtered_mean[(nt_m_std_filtered_mean[col_std_name[0]] < MAX_STD) |
                                               (nt_m_std_filtered_mean[col_std_name[1]] < MAX_STD) |
                                               (nt_m_std_filtered_mean[col_std_name[2]] < MAX_STD) |
                                               (nt_m_std_filtered_mean[col_std_name[3]] < MAX_STD)]
    nt_m_std_filtered[col_std_name] = nt_m_std_filtered[col_std_name].replace(0, MIN_STD)

    return nt_m_std_filtered


def prepare_runkos(main_dir, discard_file=None):
    """Identify the positions with variations between 0.2 to 0.8 in the training population
    and calculate the mean and std for the variation.

    """
    THRESHOLD_DROPNA = 32  # more than 40 columns should have a value not a nan.
    file_count = 0
    list_of_dfs = []
    list_dfs_var = []

    file_csv_list = []
    for file_csv in glob.glob(os.path.join(main_dir, '*.csv')):
        # Ignore the file that is given as the validation dataset
        if discard_file and file_csv == discard_file:
            continue
        dataframe = []
        dataframe = pd.read_csv(file_csv)
        # Ignore Insertions by getting only the first entity of the nucleotide at each row
        dataframe['mut'] = dataframe.apply(lambda row: row['mut'][0], axis=1)
        # Ignore insertion and deletions by merging the rows with equal pos and mut
        dataframe = dataframe.groupby(['pos', 'mut']).sum().reset_index()

        dataframe_pivot = []
        dataframe_pivot = pd.pivot_table(dataframe, index='pos', columns='mut', values='freq', aggfunc='sum')
        # Rename columns
        dataframe_pivot.rename(lambda x: "%s_%s" % (x, file_count), axis='columns', inplace=True)
        list_dfs_var.append(dataframe_pivot.copy())
        file_count += 1
        file_csv_list.append(file_csv)

    df_merged = reduce(lambda x, y: pd.merge(x, y, how='outer', right_index=True, left_index=True),
                       list_dfs_var)
    df_merged_filtered = df_merged.dropna(thresh=THRESHOLD_DROPNA)
    df_mean_std = calc_mean_std_population(df_merged_filtered, discard_file)

    return df_mean_std, list_of_dfs


def analyze_sample_variable_sites(input_file, nt_m_std_filtered):
    """Identify the sites with high variation (e.g. [0.2,0.8]) using information in nt_m_std_filtered
        for the sites with min X datapoints.
    """
    sample_pd = pd.read_csv(input_file)
    # Remove insertions and deletions
    sample_pd.apply(lambda row: row['mut'][0], axis=1)
    sample_pd = sample_pd.groupby(['pos', 'mut']).sum().reset_index()

    # Select the training sites
    sample_filtered = sample_pd[sample_pd['pos'].isin(nt_m_std_filtered.index)]

    sample_df_filtered = pd.pivot_table(sample_filtered, index='pos', columns='mut', values='freq')
    # Count number of rows as the number of training rows
    training_points = sample_df_filtered.shape[0]
    # Filter the training sites for mean in [0.2,0.8] and std < 0.25 then count
    # number of rows
    sample_df_filtered = sample_df_filtered[((sample_df_filtered['A'] > MEAN_MIN_THRESHOLD) |
                                             sample_df_filtered['A'].isnull())]
    sample_df_filtered = sample_df_filtered[((sample_df_filtered['C'] > MEAN_MIN_THRESHOLD) |
                                             sample_df_filtered['C'].isnull())]
    sample_df_filtered = sample_df_filtered[((sample_df_filtered['G'] > MEAN_MIN_THRESHOLD) |
                                             sample_df_filtered['G'].isnull())]
    sample_df_filtered = sample_df_filtered[((sample_df_filtered['T'] > MEAN_MIN_THRESHOLD) |
                                             sample_df_filtered['T'].isnull())]
    sample_df_filtered = sample_df_filtered[((sample_df_filtered['A'] < MEAN_MAX_THRESHOLD) |
                                             sample_df_filtered['A'].isnull())]
    sample_df_filtered = sample_df_filtered[((sample_df_filtered['C'] < MEAN_MAX_THRESHOLD) |
                                             sample_df_filtered['C'].isnull())]
    sample_df_filtered = sample_df_filtered[((sample_df_filtered['G'] < MEAN_MAX_THRESHOLD) |
                                             sample_df_filtered['G'].isnull())]
    sample_df_filtered = sample_df_filtered[((sample_df_filtered['T'] < MEAN_MAX_THRESHOLD) |
                                             sample_df_filtered['T'].isnull())]

    # Divide filtered positions by training positions
    valid_sample_points = sample_df_filtered.shape[0]
    score_sample = float(valid_sample_points)/training_points
    
    df_towrite = pd.merge(nt_m_std_filtered, sample_df_filtered, how='outer', right_index=True, left_index=True)
    report_file = "score_report.txt"
    with open(report_file, 'w') as rep_f:
        rep_f.write("The run control score is %f. \n" % score_sample)
        rep_f.write("%d sample positions out of %d selected positions are present.\n"
                    % (valid_sample_points, training_points))
        if score_sample < 1:
            rep_f.write("WARNING, run control score is below 1. \n")
        
        rep_f.write("\nMean and standard deviation of mutation frequencies in the training dataset followed by the mutation frequencies in the current sample in the last four columns.\n\n")
        dfAsString = df_towrite.to_string(header=True, index=True)
        rep_f.write(dfAsString)

    return score_sample


def main(filein):

    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    # Provide the training datasets directory.
    #training_dir = "runko_nt_csv_667037_training"
    training_dir = "runko_nt_csv_424736_training"
    training_full_path = os.path.join(base_dir, training_dir) #"%s/%s" % (base_dir, training_dir)
    if not os.path.exists(training_full_path):
        logging.error("The training dir %s does not exists." % training_full_path)
        sys.exit("The training dir %s does not exists." % training_full_path)

    # Identify the positions with variations between 0.2 to 0.8 in the training population
    nt_mean_std, population_df = prepare_runkos(training_full_path)
    # Estimate the score for the given sample
    analyze_sample_variable_sites(filein, nt_mean_std)
