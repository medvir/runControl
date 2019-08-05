# -*- coding: utf-8 -*-
"""Test for score calculation function.
"""

import os
import pytest
import pandas as pd
from src.runControl.run import (analyze_sample_variable_sites)

mean_std_population_csv = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                       'nt_m_std_filtered.csv')
sample_file = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           'mutations_nt_pos_ref_aa_G27-P81-runKo-582-424736-E_S12_20.csv')

assert os.path.exists(mean_std_population_csv)
assert os.path.exists(sample_file)


def test_analyze_sample_variable_sites():
    """Calculate and compare the score of a rejected sample with the expected
    score.
    """
    mean_std_population_df = pd.read_csv(mean_std_population_csv, sep='\t')
    mean_std_population_df.set_index('pos', inplace=True)
    score = analyze_sample_variable_sites(sample_file, mean_std_population_df)
    assert score == pytest.approx(0.058824, rel=1e-4)
