import pytest
import pandas as pd
import sys
import os
from tcrregex import mappers

from tcrregex.cdr3_motif import TCRMotif

@pytest.mark.skipif(sys.platform.startswith("win"), reason = "LINUX/MAC ONLY TEST WITH FILES INSTALL")
def test_TCRMotif_generates_all_tcrs():
    fn = os.path.join("tcrregex","test_files", "mouse_pairseqs_v1_parsed_seqs_probs_mq20_clones.tsv")
    clones_df_test = pd.read_csv(fn, sep="\t")
    motif = TCRMotif(clones_df = clones_df_test, organism = "mouse", chains = ["A","B"], epitopes =["PA"])
    assert isinstance(motif.all_tcrs, dict)

@pytest.mark.skipif(sys.platform.startswith("win"), reason = "LINUX/MAC ONLY TEST WITH FILES INSTALL")
def test_TCRMotif_generates_ng_tcrs():
    fn = os.path.join("tcrregex","test_files", "mouse_pairseqs_v1_parsed_seqs_probs_mq20_clones.tsv")
    clones_df_test = pd.read_csv(fn, sep="\t")
    motif = TCRMotif(clones_df = clones_df_test, organism = "mouse", chains = ["A","B"], epitopes =["PA"])
    assert isinstance(motif.ng_tcrs, dict)
    assert len(motif.ng_tcrs['B'].keys()) > 0
    assert len(motif.ng_tcrs['A'].keys()) > 0

# @pytest.mark.skipif(sys.platform.startswith("win"), reason = "LINUX/MAC ONLY TEST WITH FILES INSTALL")
# def test_integration_TCRrep_with_TCRMotif():
#     import pandas as pd
#     import tcrregex as td
#     from tcrregex import mappers
#     from tcrregex.repertoire import TCRrep
#     from tcrregex.cdr3_motif import TCRMotif
#     fn = os.path.join("tcrregex","test_files", "vdjDB_PMID28636592.tsv")
#     pd_df = pd.read_csv(fn, sep = "\t")        # 1
#     t_df = td.mappers.vdjdb_to_tcrregex2(pd_df = pd_df)               # 2
#     t_df.organism.value_counts                                       # 3
#     index_mus = t_df.organism == "MusMusculus"                       # 4
#     t_df_mus = t_df.loc[index_mus,:].copy()                          # 5

#     tr = TCRrep(cell_df = t_df_mus, organism = "mouse")              # 6

#     tr.infer_cdrs_from_v_gene(chain = 'alpha')                       # 7
#     tr.infer_cdrs_from_v_gene(chain = 'beta')                        # 8

#     tr.index_cols = ['subject', 'epitope',                           # subject and epitope
#                      'v_a_gene',  'j_a_gene', 'v_b_gene', 'j_b_gene',# gene usage
#                      'cdr3_a_aa', 'cdr3_b_aa',                       # CDR 3
#                      'cdr1_a_aa', 'cdr2_a_aa', 'pmhc_a_aa',          # alpha CDR 1, 2, and 2.5
#                      'cdr1_b_aa', 'cdr2_b_aa', 'pmhc_b_aa']          # beta CDR 1, 2, and 2.5

#     tr.deduplicate()                                                 # 10

#     motif = TCRMotif(clones_df = tr.tcr_motif_clones_df(), organism = "mouse", chains = ["A","B"], epitopes = ["PA"]) # 11
#     assert isinstance(motif.clones_df, pd.DataFrame)
#     # WAY TO SLOW TO ACTUALLY TEST MOTIF FINDING:
#     #fm = motif.find_cdr3_motifs()
#     #motif.motifs_df.head()
