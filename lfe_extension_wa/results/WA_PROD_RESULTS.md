# Washington production run (W0-W13) -- results

Auto-written by `scripts/wa_prod_report.py` at 2026-10-06 16:08Z. Outputs under `results/wa_prod/` and `results/cut/wa_ext/`; queue `scripts/wa_prod_queue.sh` (log `logs/wa_prod_queue.log`). No downloads, no deletions, no git, nothing existing overwritten. 'Scrambled control' = the same detector with each station site's data shifted to a different time of the same day; 'time-reversed control' = every template reversed in time (the main false-alarm reference, user 2026-10-03).

## In short

- **Fixed rate networks** (rule written before any run, `setup/core_sets.json`): headline E3 2015-2026 = 19 stations (B001, B003, B004, B005, B006, B007, B009, B010, B011, B013, B014, B017, B018, B019, B926, B941, B943, GMW, GNW); E2 2012-2014 = 37; E1 2007-2011 = 12.
- **Confidence (WAVE, own family, out-of-fold on fit days)**: E3 ArmB AUC 0.836 (net_mad 0.612), slope 0.26, matched/ΣP 0.32, control P>0.5 scrambled 0.004 / reversed 0.001; E3 NW AUC 0.833 (net_mad 0.658), slope 1.00, matched/ΣP 1.05, control P>0.5 scrambled 0.005 / reversed 0.001; E3 OlyA AUC 0.872 (net_mad 0.806), slope 1.05, matched/ΣP 0.83, control P>0.5 scrambled 0.001 / reversed 0.000; E3 SW AUC 0.861 (net_mad 0.739), slope 2.80, matched/ΣP 1.83, control P>0.5 scrambled 0.003 / reversed 0.000; all ArmB AUC 0.923 (net_mad 0.675), slope nan, matched/ΣP 0.56, control P>0.5 scrambled 0.003 / reversed 0.001; all NW AUC 0.838 (net_mad 0.622), slope 1.02, matched/ΣP 1.02, control P>0.5 scrambled 0.006 / reversed 0.003; all OlyA AUC 0.825 (net_mad 0.781), slope 0.91, matched/ΣP 0.84, control P>0.5 scrambled 0.001 / reversed 0.000; all SW AUC 0.844 (net_mad 0.800), slope 1.36, matched/ΣP 1.32, control P>0.5 scrambled 0.003 / reversed 0.001 (`conf/eval_dev.csv`).
- **Sealed one-look check** (2007-01 headline, 2009-05 semi-independent; `conf/sealed_bars.json`): sealed_headline|E1|own: AUC 0.856 vs net_mad 0.684; sealed_headline|E2|own: AUC 0.861 vs net_mad 0.665; sealed_headline|E3|own: AUC 0.835 vs net_mad 0.657; sealed_headline|all|own: AUC 0.859 vs net_mad 0.663; sealed_semi|E1|own: AUC 0.765 vs net_mad 0.681; sealed_semi|E2|own: AUC 0.766 vs net_mad 0.708; sealed_semi|E3|own: AUC 0.784 vs net_mad 0.705; sealed_semi|all|own: AUC 0.771 vs net_mad 0.708.
- **PNET gate**: evaluated; adopt = True (dAUC [0.03246201858616095, 0.0244808781223639, 0.04044027687317215], held-out-station block [0.006228014029290867, 0.002224937452258241, 0.01011771373835492]; `conf/pnet_gate.json`).
- **False share at P >= 0.75 (core network)**, time-reversed / scrambled control vs data: ArmB quiet: 4/463 of 641 (0.006 / 0.722); ArmB tremor60: 0/31 of 565 (0.000 / 0.055); NW quiet: 1012/634 of 2956 (0.342 / 0.214); NW tremor60: 130/311 of 29605 (0.004 / 0.011); OlyA quiet: 370/1402 of 1813 (0.204 / 0.773); OlyA tremor60: 33/191 of 1306 (0.025 / 0.146); SW quiet: 6/2 of 0 (nan / nan); SW tremor60: 1/1 of 226 (0.004 / 0.004) (`cutoff_check/false_share_by_cutoff.csv`).
- **Refit risk** (P_own_WAVE >= 0.75 but the older model < 0.2): 0.0391 -> window floor P_max >= 0.1 (`conf/refit_risk.json`).
- Windows: pilot only ({'pilot_month': '2015-10', 'cut_GB': 10.03, 'raw_GB': 23.248, 'cut_over_raw': 0.43143496214728144, 'raw_total_TB': 8.415315553792, 'extrapolated_TB': 3.6306613474076803, 'n_windows': 31308568, 'minutes': 14.8, 'limit_TB': 3.0, 'note': 'the pilot month is the busiest production month; full-store estimate = pilot cut/raw ratio x all raw in the roots', 'pass': False}); full cut not run.
- **Delete-candidate list** (nothing deleted): WA candidates 1363 entries, 7.759 TB by kind {'cache_svi': 1399.0, 'cache_svi_r1': 2.8, 'cache_wa_prod': 5912.6, 'cache_wa_trial': 444.4} GB (`delete/`).
- Caveat: P = probability of matching the array catalogue on its scanned days. Between ETS windows it is a lower bound on 'real', not a measure of it; inter-ETS rates are a control-subtracted bracket.

## W0-W3 freeze, completeness, health, caches

- Frozen inputs: `w0/frozen_inputs.csv`; re-hash at the end: ['/home/jovyan/tremorferometry/lfe_extension_wa/results/wa_prod/w0/verify_20261006T1608.csv']
- Cache determinism: 20/20 station-days rebuilt bit-identical (`cache/determinism.csv`).
- guard cache_p1: [guard:cache_p1] station-days ok 195398, error/missing 0 (0.00%); all {'ok': 195398, 'nodata': 109}
- guard cache_p2: [guard:cache_p2] station-days ok 93032, error/missing 0 (0.00%); all {'ok': 93032, 'nodata': 1777}
- guard hf: [guard:hf] station-days ok 75996, error/missing 0 (0.00%); all {'ok': 75996}
- guard mf_prod: [guard:mf_prod] station-days ok 208249, error/missing 0 (0.00%); all {'nodata': 379155, 'ok': 208249, 'notemplate': 42784, 'dead': 8062, 'excluded': 414}
- guard mf_rich: [guard:mf_rich] station-days ok 125973, error/missing 0 (0.00%); all {'nodata': 253749, 'ok': 125973, 'notemplate': 29749, 'dead': 3524, 'excluded': 197}
- guard mf_truth: [guard:mf_truth] station-days ok 10435, error/missing 0 (0.00%); all {'nodata': 24564, 'ok': 10435, 'dead': 262, 'notemplate': 53, 'excluded': 19}
- guard mf_truth_e12: [guard:mf_truth_e12] station-days ok 10435, error/missing 0 (0.00%); all {'nodata': 24564, 'ok': 10435, 'dead': 262, 'notemplate': 53, 'excluded': 19}
- guard pnet_prod: [guard:pnet_prod] station-days ok 208611, error/missing 0 (0.00%); all {'ok': 208611}
- guard pnet_truth: [guard:pnet_truth] station-days ok 10435, error/missing 0 (0.00%); all {'ok': 10435}
- B014 2023 check: {'pre_tested': 10, 'pre_ok': 3, 'post_tested': 8, 'post_ok': 4, 'rule': "re-template only if the post-outage share of 'ok' tested months is below the pre-outage share minus 0.25", 'pre_share': 0.3, 'post_share': 0.5, 'retemplate': False}

## W4 templates

- W4a library: 19622 (family, node) templates on 116 nodes; by tag {'build': 11358, 'v10b': 7636, 't2012new': 258, 'step4b_t2013': 214, 'step4b_t2012': 110, 'step4b_t2014': 46} (`setup/tpl_library.csv`; 'v10b' = the step-4b templates rebuilt from the V2010 build half so A2011 stays leak-free).
- Node map (W2 verdicts):

| netsta   | node                       | start      | kind          | src_node                  | w2_verdict                    |   cc_after |   ctl95_after |
|:---------|:---------------------------|:-----------|:--------------|:--------------------------|:------------------------------|-----------:|--------------:|
| CN.PGC   | CN.PGC..BH.40@2015-01-12   | 2015-01-12 | pending_retpl | CN.PGC..BH.40@2000-06-28  | AT CONTROL LEVEL              |      0.074 |         0.088 |
| CN.PGC   | CN.PGC..HH.100@2017-08-10  | 2017-08-10 | carry         | CN.PGC..BH.40@2000-06-28  | above control (carry-over OK) |      0.097 |         0.082 |
| CN.VGZ   | CN.VGZ..HH.100@2017-03-03  | 2017-03-03 | pending_retpl | CN.VGZ..HH.100@2009-07-31 | AT CONTROL LEVEL              |      0.108 |         0.119 |
| UW.DOSE  | UW.DOSE..HH.100@2019-05-07 | 2019-05-07 | pending_retpl | UW.DOSE..BH.40@2010-07-21 | AT CONTROL LEVEL              |      0.087 |         0.089 |
| UW.GNW   | UW.GNW..HH.100@2019-05-07  | 2019-05-07 | carry         | UW.GNW..BH.40@2010-12-01  | above control (carry-over OK) |      0.123 |         0.079 |
| UW.HDW   | UW.HDW..EH.100@2008-10-01  | 2008-10-01 | untested      | UW.HDW..EH.100@1998-05-21 | untested                      |      0.162 |       nan     |
| UW.HDW   | UW.HDW..EH.100@2020-12-04  | 2020-12-04 | pending_retpl | UW.HDW..EH.100@1998-05-21 | AT CONTROL LEVEL              |      0.129 |         0.145 |
| UW.HDW   | UW.HDW..HH.100@2022-09-01  | 2022-09-01 | carry         | UW.HDW..EH.100@1998-05-21 | above control (carry-over OK) |      0.165 |         0.132 |
| UW.LRIV  | UW.LRIV..HH.100@2019-05-30 | 2019-05-30 | carry         | UW.LRIV..BH.40@2009-08-27 | above control (carry-over OK) |      0.082 |         0.079 |

- W4b: 132 target nodes; 2409 (family, node) stacks, 156 used on 34 nodes; median cc of a used template vs the scrambled-control-time stack 0.21 (`templates/tpl_qc_w4b.csv`).

| kind                     | netsta   |   n_templates |   med_cc_ctl |
|:-------------------------|:---------|--------------:|-------------:|
| new_node_t2b             | UW.BES3  |             2 |        0.189 |
| new_node_t2b             | UW.BES4  |             2 |        0.176 |
| new_node_t2b             | UW.BES5  |             2 |        0.218 |
| new_node_t2b             | UW.BHW   |             2 |        0.233 |
| new_node_t2b             | UW.BST01 |             3 |        0.189 |
| new_node_t2b             | UW.BST02 |             3 |        0.202 |
| new_node_t2b             | UW.BST03 |             4 |        0.190 |
| new_node_t2b             | UW.BST04 |             6 |        0.187 |
| new_node_t2b             | UW.BST05 |             3 |        0.190 |
| new_node_t2b             | UW.BST07 |             6 |        0.237 |
| new_node_t2b             | UW.BST08 |             6 |        0.213 |
| new_node_t2b             | UW.BST09 |             2 |        0.225 |
| new_node_t2b             | UW.BST10 |             4 |        0.222 |
| new_node_t2b             | UW.BST11 |             5 |        0.219 |
| new_node_t2b             | UW.BST16 |             2 |        0.202 |
| new_node_t2b             | UW.BST20 |            17 |        0.165 |
| new_node_t2b             | UW.BST21 |             9 |        0.218 |
| new_node_t2b             | UW.BST22 |             7 |        0.188 |
| new_node_t2b             | UW.BST23 |            10 |        0.151 |
| new_node_t2b             | UW.CHIMA |             2 |        0.306 |
| new_node_t2b             | UW.COYL  |             3 |        0.295 |
| new_node_t2b             | UW.EQUIL |             9 |        0.153 |
| new_node_t2b             | UW.HANS  |             2 |        0.345 |
| new_node_t2b             | UW.HTW   |             1 |      nan     |
| new_node_t2b             | UW.KNGST |             1 |        0.191 |
| new_node_t2b             | UW.ODEER |             1 |        0.458 |
| new_node_t2b             | UW.ODUC  |             1 |        0.188 |
| new_node_t2b             | UW.OTR   |             2 |        0.244 |
| new_node_t2b             | UW.PRIN  |            10 |        0.352 |
| new_node_t2b             | UW.SP2   |             5 |        0.243 |
| new_node_t2b             | XD.KRES  |            10 |        0.238 |
| retemplate_control_level | CN.VGZ   |             5 |        0.195 |
| retemplate_control_level | UW.DOSE  |             9 |        0.163 |

## W5-W8 confidence

| net   | reg   | label   | model   |      n |   n_pos |   auc |   auc_net_mad |   slope |   matched_over_sumP |   ctl_frac_P_gt05 |   rev_frac_P_gt05 |
|:------|:------|:--------|:--------|-------:|--------:|------:|--------------:|--------:|--------------------:|------------------:|------------------:|
| E1    | ArmB  | own     | WAVE    |   2638 |     119 | 0.836 |         0.627 |   0.405 |               0.395 |             0.003 |             0.002 |
| E1    | ArmB  | own     | FULL    |   2638 |     119 | 0.807 |         0.627 |   1.067 |               1.021 |             0.000 |             0.000 |
| E1    | ArmB  | own     | MF_ONLY |   2638 |     119 | 0.571 |         0.627 |   0.249 |               0.307 |             0.023 |             0.001 |
| E1    | ArmB  | nbr3    | WAVE    |   2638 |     133 | 0.774 |         0.620 |   0.362 |               0.265 |             0.005 |             0.003 |
| E1    | ArmB  | nbr3    | FULL    |   2638 |     133 | 0.796 |         0.620 |   1.049 |               0.928 |             0.000 |             0.000 |
| E1    | ArmB  | nbr3    | MF_ONLY |   2638 |     133 | 0.562 |         0.620 |   0.200 |               0.284 |             0.031 |             0.001 |
| E1    | ArmB  | any     | WAVE    |   2638 |     742 | 0.634 |         0.608 |   0.523 |               0.805 |             0.021 |             0.022 |
| E1    | ArmB  | any     | FULL    |   2638 |     742 | 0.688 |         0.608 |   0.718 |               0.902 |             0.040 |             0.014 |
| E1    | ArmB  | any     | MF_ONLY |   2638 |     742 | 0.549 |         0.608 |   0.407 |               0.847 |             0.083 |             0.020 |
| E1    | NW    | own     | WAVE    | 113154 |   23390 | 0.839 |         0.632 |   1.022 |               1.071 |             0.007 |             0.005 |
| E1    | NW    | own     | FULL    | 113154 |   23390 | 0.879 |         0.632 |   1.022 |               1.020 |             0.014 |             0.003 |
| E1    | NW    | own     | MF_ONLY | 113154 |   23390 | 0.778 |         0.632 |   1.245 |               1.376 |             0.006 |             0.005 |
| E1    | NW    | nbr3    | WAVE    | 113154 |   23654 | 0.842 |         0.631 |   1.040 |               1.023 |             0.008 |             0.005 |
| E1    | NW    | nbr3    | FULL    | 113154 |   23654 | 0.879 |         0.631 |   1.035 |               1.030 |             0.014 |             0.003 |
| E1    | NW    | nbr3    | MF_ONLY | 113154 |   23654 | 0.791 |         0.631 |   1.263 |               1.287 |             0.007 |             0.006 |
| E1    | NW    | any     | WAVE    | 113154 |   39025 | 0.789 |         0.558 |   1.058 |               1.024 |             0.041 |             0.031 |
| E1    | NW    | any     | FULL    | 113154 |   39025 | 0.824 |         0.558 |   1.049 |               1.033 |             0.049 |             0.032 |
| E1    | NW    | any     | MF_ONLY | 113154 |   39025 | 0.744 |         0.558 |   1.209 |               1.169 |             0.025 |             0.024 |
| E1    | OlyA  | own     | WAVE    | 104610 |    4111 | 0.851 |         0.778 |   0.908 |               0.751 |             0.000 |             0.000 |
| E1    | OlyA  | own     | FULL    | 104610 |    4111 | 0.896 |         0.778 |   1.043 |               1.057 |             0.001 |             0.000 |
| E1    | OlyA  | own     | MF_ONLY | 104610 |    4111 | 0.648 |         0.778 |   0.173 |               0.415 |             0.006 |             0.007 |
| E1    | OlyA  | nbr3    | WAVE    | 104610 |   12904 | 0.766 |         0.675 |   1.055 |               1.005 |             0.001 |             0.001 |
| E1    | OlyA  | nbr3    | FULL    | 104610 |   12904 | 0.831 |         0.675 |   0.974 |               1.011 |             0.008 |             0.000 |
| E1    | OlyA  | nbr3    | MF_ONLY | 104610 |   12904 | 0.622 |         0.675 |   0.454 |               0.784 |             0.009 |             0.009 |
| E1    | OlyA  | any     | WAVE    | 104610 |   24567 | 0.697 |         0.659 |   1.047 |               0.989 |             0.006 |             0.005 |
| E1    | OlyA  | any     | FULL    | 104610 |   24567 | 0.780 |         0.659 |   0.988 |               0.997 |             0.025 |             0.007 |
| E1    | OlyA  | any     | MF_ONLY | 104610 |   24567 | 0.606 |         0.659 |   0.543 |               0.851 |             0.026 |             0.029 |
| E1    | SW    | own     | WAVE    |    234 |      76 | 0.803 |         0.709 |   2.091 |               1.793 |             0.003 |             0.001 |
| E1    | SW    | own     | FULL    |    234 |      76 | 0.883 |         0.709 |   1.475 |               1.093 |             0.012 |             0.001 |
| E1    | SW    | own     | MF_ONLY |    234 |      76 | 0.757 |         0.709 | nan     |               2.454 |             0.000 |             0.000 |
| E1    | SW    | nbr3    | WAVE    |    234 |      78 | 0.848 |         0.723 |   2.263 |               1.446 |             0.006 |             0.001 |
| E1    | SW    | nbr3    | FULL    |    234 |      78 | 0.875 |         0.723 |   1.418 |               1.079 |             0.015 |             0.002 |
| E1    | SW    | nbr3    | MF_ONLY |    234 |      78 | 0.789 |         0.723 | nan     |               2.257 |             0.001 |             0.000 |
| E1    | SW    | any     | WAVE    |    234 |      99 | 0.762 |         0.684 |   1.896 |               1.324 |             0.006 |             0.003 |
| E1    | SW    | any     | FULL    |    234 |      99 | 0.778 |         0.684 |   1.676 |               1.160 |             0.024 |             0.006 |
| E1    | SW    | any     | MF_ONLY |    234 |      99 | 0.634 |         0.684 |   0.824 |               1.174 |             0.011 |             0.005 |
| E2    | ArmB  | own     | WAVE    |   1653 |      94 | 0.851 |         0.627 |   0.482 |               0.432 |             0.004 |             0.001 |
| E2    | ArmB  | own     | FULL    |   1653 |      94 | 0.835 |         0.627 |   1.755 |               1.329 |             0.000 |             0.000 |
| E2    | ArmB  | own     | MF_ONLY |   1653 |      94 | 0.629 |         0.627 |   0.288 |               0.412 |             0.022 |             0.000 |
| E2    | ArmB  | nbr3    | WAVE    |   1653 |     101 | 0.756 |         0.614 |   0.339 |               0.267 |             0.006 |             0.001 |
| E2    | ArmB  | nbr3    | FULL    |   1653 |     101 | 0.790 |         0.614 |   1.681 |               1.271 |             0.000 |             0.000 |
| E2    | ArmB  | nbr3    | MF_ONLY |   1653 |     101 | 0.661 |         0.614 |   0.497 |               0.335 |             0.025 |             0.000 |
| E2    | ArmB  | any     | WAVE    |   1653 |     560 | 0.633 |         0.578 |   0.618 |               0.823 |             0.031 |             0.010 |
| E2    | ArmB  | any     | FULL    |   1653 |     560 | 0.684 |         0.578 |   0.734 |               0.976 |             0.036 |             0.001 |
| E2    | ArmB  | any     | MF_ONLY |   1653 |     560 | 0.531 |         0.578 |   0.329 |               1.024 |             0.052 |             0.014 |
| E2    | NW    | own     | WAVE    | 133397 |   26709 | 0.842 |         0.654 |   0.999 |               1.040 |             0.007 |             0.002 |
| E2    | NW    | own     | FULL    | 133397 |   26709 | 0.882 |         0.654 |   1.012 |               0.974 |             0.014 |             0.001 |
| E2    | NW    | own     | MF_ONLY | 133397 |   26709 | 0.737 |         0.654 |   1.028 |               1.177 |             0.007 |             0.004 |
| E2    | NW    | nbr3    | WAVE    | 133397 |   27017 | 0.846 |         0.653 |   1.026 |               1.007 |             0.007 |             0.002 |
| E2    | NW    | nbr3    | FULL    | 133397 |   27017 | 0.881 |         0.653 |   1.028 |               0.989 |             0.012 |             0.001 |
| E2    | NW    | nbr3    | MF_ONLY | 133397 |   27017 | 0.758 |         0.653 |   1.112 |               1.089 |             0.009 |             0.005 |
| E2    | NW    | any     | WAVE    | 133397 |   45736 | 0.790 |         0.575 |   1.040 |               1.008 |             0.040 |             0.016 |
| E2    | NW    | any     | FULL    | 133397 |   45736 | 0.824 |         0.575 |   1.037 |               0.986 |             0.051 |             0.018 |
| E2    | NW    | any     | MF_ONLY | 133397 |   45736 | 0.717 |         0.575 |   1.111 |               1.068 |             0.028 |             0.023 |
| E2    | OlyA  | own     | WAVE    | 104975 |    4099 | 0.864 |         0.798 |   0.932 |               0.747 |             0.001 |             0.000 |
| E2    | OlyA  | own     | FULL    | 104975 |    4099 | 0.899 |         0.798 |   1.067 |               1.039 |             0.001 |             0.000 |
| E2    | OlyA  | own     | MF_ONLY | 104975 |    4099 | 0.744 |         0.798 |   0.444 |               0.474 |             0.005 |             0.022 |
| E2    | OlyA  | nbr3    | WAVE    | 104975 |   12934 | 0.779 |         0.688 |   1.014 |               0.969 |             0.002 |             0.001 |
| E2    | OlyA  | nbr3    | FULL    | 104975 |   12934 | 0.838 |         0.688 |   0.956 |               0.982 |             0.011 |             0.000 |
| E2    | OlyA  | nbr3    | MF_ONLY | 104975 |   12934 | 0.683 |         0.688 |   0.843 |               0.806 |             0.008 |             0.025 |
| E2    | OlyA  | any     | WAVE    | 104975 |   24714 | 0.725 |         0.671 |   1.040 |               0.942 |             0.010 |             0.003 |
| E2    | OlyA  | any     | FULL    | 104975 |   24714 | 0.787 |         0.671 |   0.971 |               0.973 |             0.030 |             0.002 |
| E2    | OlyA  | any     | MF_ONLY | 104975 |   24714 | 0.661 |         0.671 |   0.928 |               0.858 |             0.024 |             0.059 |
| E2    | SW    | own     | WAVE    |    331 |     169 | 0.863 |         0.765 |   1.408 |               1.329 |             0.003 |             0.001 |
| E2    | SW    | own     | FULL    |    331 |     169 | 0.875 |         0.765 |   1.232 |               1.157 |             0.005 |             0.001 |
| E2    | SW    | own     | MF_ONLY |    331 |     169 | 0.680 |         0.765 |   1.349 |               1.861 |             0.001 |             0.000 |
| E2    | SW    | nbr3    | WAVE    |    331 |     172 | 0.867 |         0.770 |   1.461 |               1.329 |             0.003 |             0.001 |
| E2    | SW    | nbr3    | FULL    |    331 |     172 | 0.878 |         0.770 |   1.335 |               1.187 |             0.008 |             0.001 |
| E2    | SW    | nbr3    | MF_ONLY |    331 |     172 | 0.726 |         0.770 |   1.456 |               1.787 |             0.001 |             0.001 |
| E2    | SW    | any     | WAVE    |    331 |     197 | 0.831 |         0.747 |   1.450 |               1.318 |             0.009 |             0.007 |
| E2    | SW    | any     | FULL    |    331 |     197 | 0.844 |         0.747 |   1.328 |               1.259 |             0.018 |             0.005 |
| E2    | SW    | any     | MF_ONLY |    331 |     197 | 0.699 |         0.747 |   1.157 |               1.326 |             0.008 |             0.009 |
| E3    | ArmB  | own     | WAVE    |   2317 |      70 | 0.836 |         0.612 |   0.259 |               0.324 |             0.004 |             0.001 |
| E3    | ArmB  | own     | FULL    |   2317 |      70 | 0.831 |         0.612 | nan     |               0.985 |             0.000 |             0.000 |
| E3    | ArmB  | own     | MF_ONLY |   2317 |      70 | 0.641 |         0.612 | nan     |               0.370 |             0.021 |             0.000 |
| E3    | ArmB  | nbr3    | WAVE    |   2317 |      80 | 0.740 |         0.619 |   0.242 |               0.199 |             0.005 |             0.001 |
| E3    | ArmB  | nbr3    | FULL    |   2317 |      80 | 0.775 |         0.619 |   1.104 |               0.936 |             0.000 |             0.000 |
| E3    | ArmB  | nbr3    | MF_ONLY |   2317 |      80 | 0.637 |         0.619 | nan     |               0.292 |             0.024 |             0.000 |
| E3    | ArmB  | any     | WAVE    |   2317 |     706 | 0.623 |         0.579 |   0.632 |               0.869 |             0.028 |             0.008 |
| E3    | ArmB  | any     | FULL    |   2317 |     706 | 0.698 |         0.579 |   0.781 |               0.966 |             0.034 |             0.003 |

By truth_array_frac (dev, WAVE own):

| net   | reg   | truth_array_frac   |     n |   n_pos |   auc |   auc_net_mad |
|:------|:------|:-------------------|------:|--------:|------:|--------------:|
| E1    | NW    | <=0.45             | 24060 |    6754 | 0.827 |         0.713 |
| E1    | NW    | 0.45-0.55          | 38074 |    7711 | 0.828 |         0.665 |
| E1    | NW    | 0.55-0.65          | 38067 |    6413 | 0.848 |         0.557 |
| E1    | NW    | >0.65              | 12900 |    2512 | 0.843 |         0.643 |
| E1    | SW    | >0.65              |   234 |      76 | 0.803 |         0.709 |
| E2    | NW    | <=0.45             | 35313 |    8321 | 0.849 |         0.725 |
| E2    | NW    | 0.45-0.55          | 46202 |    8944 | 0.840 |         0.677 |
| E2    | NW    | 0.55-0.65          | 38896 |    6855 | 0.839 |         0.575 |
| E2    | NW    | >0.65              | 12901 |    2589 | 0.825 |         0.636 |
| E2    | SW    | >0.65              |   331 |     169 | 0.863 |         0.765 |
| E3    | NW    | <=0.45             | 44060 |    8554 | 0.850 |         0.723 |
| E3    | NW    | 0.45-0.55          | 56064 |    9249 | 0.834 |         0.685 |
| E3    | NW    | 0.55-0.65          | 45715 |    7287 | 0.817 |         0.595 |
| E3    | NW    | >0.65              | 17696 |    2686 | 0.805 |         0.587 |
| E3    | SW    | >0.65              |   192 |      73 | 0.861 |         0.739 |
| all   | NW    | <=0.45             | 23756 |    6794 | 0.837 |         0.701 |
| all   | NW    | 0.45-0.55          | 31224 |    6998 | 0.823 |         0.629 |
| all   | NW    | 0.55-0.65          | 31851 |    5691 | 0.838 |         0.565 |
| all   | NW    | >0.65              | 11108 |    2163 | 0.842 |         0.625 |
| all   | SW    | >0.65              |   296 |     149 | 0.844 |         0.800 |

Olympic leave-region-out (credible if within 0.03 AUC of out-of-fold):

| label   | reg   | net   |      n |   n_pos |   auc_oof |   auc_lro |   diff | credible   |
|:--------|:------|:------|-------:|--------:|----------:|----------:|-------:|:-----------|
| own     | OlyA  | E1    | 104610 |    4111 |     0.851 |     0.605 | -0.246 | False      |
| own     | OlyA  | E2    | 104975 |    4099 |     0.864 |     0.635 | -0.228 | False      |
| own     | OlyA  | E3    | 113631 |    4161 |     0.872 |     0.672 | -0.200 | False      |
| own     | OlyA  | all   |  80203 |    3783 |     0.825 |     0.598 | -0.228 | False      |
| own     | ArmB  | E1    |   2638 |     119 |     0.836 |     0.823 | -0.013 | True       |
| own     | ArmB  | E2    |   1653 |      94 |     0.851 |     0.863 |  0.012 | True       |
| own     | ArmB  | E3    |   2317 |      70 |     0.836 |     0.845 |  0.009 | True       |
| own     | ArmB  | all   |    331 |      30 |     0.923 |     0.931 |  0.009 | True       |
| nbr3    | OlyA  | E1    | 104610 |   12904 |     0.766 |     0.573 | -0.193 | False      |
| nbr3    | OlyA  | E2    | 104975 |   12934 |     0.779 |     0.604 | -0.175 | False      |
| nbr3    | OlyA  | E3    | 113631 |   13397 |     0.790 |     0.628 | -0.163 | False      |
| nbr3    | OlyA  | all   |  80203 |   11143 |     0.753 |     0.561 | -0.192 | False      |
| nbr3    | ArmB  | E1    |   2638 |     133 |     0.774 |     0.804 |  0.031 | False      |
| nbr3    | ArmB  | E2    |   1653 |     101 |     0.756 |     0.855 |  0.098 | False      |
| nbr3    | ArmB  | E3    |   2317 |      80 |     0.740 |     0.831 |  0.091 | False      |
| nbr3    | ArmB  | all   |    331 |      33 |     0.858 |     0.902 |  0.044 | False      |
| any     | OlyA  | E1    | 104610 |   24567 |     0.697 |     0.525 | -0.171 | False      |
| any     | OlyA  | E2    | 104975 |   24714 |     0.725 |     0.590 | -0.134 | False      |
| any     | OlyA  | E3    | 113631 |   25751 |     0.722 |     0.608 | -0.114 | False      |
| any     | OlyA  | all   |  80203 |   21009 |     0.723 |     0.585 | -0.137 | False      |
| any     | ArmB  | E1    |   2638 |     742 |     0.634 |     0.615 | -0.019 | True       |
| any     | ArmB  | E2    |   1653 |     560 |     0.633 |     0.615 | -0.018 | True       |
| any     | ArmB  | E3    |   2317 |     706 |     0.623 |     0.604 | -0.019 | True       |
| any     | ArmB  | all   |    331 |     120 |     0.719 |     0.624 | -0.095 | False      |

Sealed (one look):

| eval            | net   | reg   | label   | model   |     n |   n_pos |   auc |   auc_net_mad |   slope |   matched_over_sumP |   ctl_frac_P_gt05 |   rev_frac_P_gt05 |
|:----------------|:------|:------|:--------|:--------|------:|--------:|------:|--------------:|--------:|--------------------:|------------------:|------------------:|
| sealed_headline | E1    | NW    | own     | WAVE    | 16736 |    3703 | 0.866 |         0.685 |   1.073 |               0.924 |             0.013 |             0.009 |
| sealed_headline | E1    | NW    | own     | FULL    | 16736 |    3703 | 0.893 |         0.685 |   1.066 |               0.877 |             0.021 |             0.004 |
| sealed_headline | E1    | NW    | own     | MF_ONLY | 16736 |    3703 | 0.774 |         0.685 |   1.280 |               1.271 |             0.009 |             0.002 |
| sealed_headline | E1    | NW    | nbr3    | WAVE    | 16736 |    3753 | 0.871 |         0.684 |   1.082 |               0.929 |             0.013 |             0.011 |
| sealed_headline | E1    | NW    | nbr3    | FULL    | 16736 |    3753 | 0.895 |         0.684 |   1.073 |               0.864 |             0.020 |             0.004 |
| sealed_headline | E1    | NW    | nbr3    | MF_ONLY | 16736 |    3753 | 0.778 |         0.684 |   1.299 |               1.268 |             0.011 |             0.003 |
| sealed_headline | E1    | NW    | any     | WAVE    | 16736 |    5550 | 0.800 |         0.627 |   1.039 |               0.897 |             0.070 |             0.050 |
| sealed_headline | E1    | NW    | any     | FULL    | 16736 |    5550 | 0.821 |         0.627 |   1.034 |               0.875 |             0.086 |             0.029 |
| sealed_headline | E1    | NW    | any     | MF_ONLY | 16736 |    5550 | 0.729 |         0.627 |   1.325 |               1.046 |             0.036 |             0.019 |
| sealed_headline | E1    | SW    | own     | WAVE    |   707 |     215 | 0.760 |         0.752 |   1.976 |               2.194 |             0.004 |             0.000 |
| sealed_headline | E1    | SW    | own     | FULL    |   707 |     215 | 0.844 |         0.752 |   1.617 |               1.324 |             0.061 |             0.002 |
| sealed_headline | E1    | SW    | own     | MF_ONLY |   707 |     215 | 0.784 |         0.752 |   1.867 |               2.547 |             0.012 |             0.000 |
| sealed_headline | E1    | SW    | nbr3    | WAVE    |   707 |     229 | 0.745 |         0.750 |   1.991 |               2.457 |             0.001 |             0.002 |
| sealed_headline | E1    | SW    | nbr3    | FULL    |   707 |     229 | 0.810 |         0.750 |   1.587 |               1.377 |             0.045 |             0.007 |
| sealed_headline | E1    | SW    | nbr3    | MF_ONLY |   707 |     229 | 0.792 |         0.750 |   2.249 |               2.637 |             0.015 |             0.000 |
| sealed_headline | E1    | SW    | any     | WAVE    |   707 |     268 | 0.743 |         0.738 |   1.793 |               1.489 |             0.016 |             0.007 |
| sealed_headline | E1    | SW    | any     | FULL    |   707 |     268 | 0.763 |         0.738 |   1.738 |               1.295 |             0.071 |             0.009 |
| sealed_headline | E1    | SW    | any     | MF_ONLY |   707 |     268 | 0.764 |         0.738 |   1.974 |               1.294 |             0.036 |             0.000 |
| sealed_headline | E2    | NW    | own     | WAVE    | 14623 |    3305 | 0.875 |         0.667 |   1.087 |               0.993 |             0.010 |             0.010 |
| sealed_headline | E2    | NW    | own     | FULL    | 14623 |    3305 | 0.894 |         0.667 |   1.081 |               0.886 |             0.022 |             0.004 |
| sealed_headline | E2    | NW    | own     | MF_ONLY | 14623 |    3305 | 0.701 |         0.667 |   1.741 |               2.044 |             0.000 |             0.001 |
| sealed_headline | E2    | NW    | nbr3    | WAVE    | 14623 |    3347 | 0.877 |         0.665 |   1.108 |               0.946 |             0.012 |             0.013 |
| sealed_headline | E2    | NW    | nbr3    | FULL    | 14623 |    3347 | 0.898 |         0.665 |   1.084 |               0.868 |             0.022 |             0.005 |
| sealed_headline | E2    | NW    | nbr3    | MF_ONLY | 14623 |    3347 | 0.731 |         0.665 |   1.979 |               1.717 |             0.000 |             0.001 |
| sealed_headline | E2    | NW    | any     | WAVE    | 14623 |    4953 | 0.799 |         0.611 |   1.072 |               0.947 |             0.048 |             0.044 |
| sealed_headline | E2    | NW    | any     | FULL    | 14623 |    4953 | 0.823 |         0.611 |   1.053 |               0.926 |             0.071 |             0.026 |
| sealed_headline | E2    | NW    | any     | MF_ONLY | 14623 |    4953 | 0.663 |         0.611 |   1.516 |               1.286 |             0.008 |             0.006 |
| sealed_headline | E2    | SW    | own     | WAVE    |   684 |     231 | 0.790 |         0.761 |   2.326 |               3.068 |             0.000 |             0.000 |
| sealed_headline | E2    | SW    | own     | FULL    |   684 |     231 | 0.893 |         0.761 |   1.825 |               1.481 |             0.052 |             0.005 |
| sealed_headline | E2    | SW    | own     | MF_ONLY |   684 |     231 | 0.769 |         0.761 | nan     |               3.132 |             0.003 |             0.000 |
| sealed_headline | E2    | SW    | nbr3    | WAVE    |   684 |     243 | 0.840 |         0.762 |   2.594 |               2.519 |             0.001 |             0.000 |
| sealed_headline | E2    | SW    | nbr3    | FULL    |   684 |     243 | 0.874 |         0.762 |   1.652 |               1.528 |             0.049 |             0.010 |
| sealed_headline | E2    | SW    | nbr3    | MF_ONLY |   684 |     243 | 0.765 |         0.762 | nan     |               3.008 |             0.001 |             0.000 |
| sealed_headline | E2    | SW    | any     | WAVE    |   684 |     289 | 0.792 |         0.739 |   2.506 |               1.682 |             0.013 |             0.008 |
| sealed_headline | E2    | SW    | any     | FULL    |   684 |     289 | 0.838 |         0.739 |   2.113 |               1.547 |             0.055 |             0.008 |
| sealed_headline | E2    | SW    | any     | MF_ONLY |   684 |     289 | 0.651 |         0.739 |   1.706 |               1.411 |             0.023 |             0.000 |
| sealed_headline | E3    | NW    | own     | WAVE    | 17872 |    3616 | 0.846 |         0.660 |   1.047 |               0.951 |             0.011 |             0.008 |
| sealed_headline | E3    | NW    | own     | FULL    | 17872 |    3616 | 0.886 |         0.660 |   1.049 |               0.905 |             0.021 |             0.002 |
| sealed_headline | E3    | NW    | own     | MF_ONLY | 17872 |    3616 | 0.753 |         0.660 |   1.025 |               1.227 |             0.017 |             0.022 |
| sealed_headline | E3    | NW    | nbr3    | WAVE    | 17872 |    3675 | 0.850 |         0.658 |   1.048 |               0.921 |             0.012 |             0.011 |
| sealed_headline | E3    | NW    | nbr3    | FULL    | 17872 |    3675 | 0.890 |         0.658 |   1.061 |               0.890 |             0.021 |             0.003 |
| sealed_headline | E3    | NW    | nbr3    | MF_ONLY | 17872 |    3675 | 0.752 |         0.658 |   1.012 |               1.223 |             0.021 |             0.028 |
| sealed_headline | E3    | NW    | any     | WAVE    | 17872 |    5928 | 0.757 |         0.588 |   0.982 |               0.899 |             0.072 |             0.059 |
| sealed_headline | E3    | NW    | any     | FULL    | 17872 |    5928 | 0.794 |         0.588 |   0.980 |               0.918 |             0.095 |             0.040 |
| sealed_headline | E3    | NW    | any     | MF_ONLY | 17872 |    5928 | 0.693 |         0.588 |   0.958 |               1.044 |             0.078 |             0.072 |
| sealed_headline | E3    | SW    | own     | WAVE    |   649 |     251 | 0.681 |         0.731 |   1.441 |               2.481 |             0.003 |             0.000 |
| sealed_headline | E3    | SW    | own     | FULL    |   649 |     251 | 0.808 |         0.731 |   1.509 |               1.561 |             0.084 |             0.000 |
| sealed_headline | E3    | SW    | own     | MF_ONLY |   649 |     251 | 0.785 |         0.731 |   1.870 |               2.914 |             0.033 |             0.000 |
| sealed_headline | E3    | SW    | nbr3    | WAVE    |   649 |     267 | 0.670 |         0.723 | nan     |               3.240 |             0.003 |             0.000 |
| sealed_headline | E3    | SW    | nbr3    | FULL    |   649 |     267 | 0.797 |         0.723 |   1.716 |               1.666 |             0.064 |             0.000 |
| sealed_headline | E3    | SW    | nbr3    | MF_ONLY |   649 |     267 | 0.786 |         0.723 |   1.944 |               2.997 |             0.037 |             0.000 |
| sealed_headline | E3    | SW    | any     | WAVE    |   649 |     332 | 0.628 |         0.695 |   1.243 |               2.019 |             0.011 |             0.009 |
| sealed_headline | E3    | SW    | any     | FULL    |   649 |     332 | 0.713 |         0.695 |   1.510 |               1.692 |             0.080 |             0.000 |
| sealed_headline | E3    | SW    | any     | MF_ONLY |   649 |     332 | 0.749 |         0.695 |   1.461 |               1.707 |             0.121 |             0.004 |
| sealed_headline | all   | NW    | own     | WAVE    | 14152 |    3239 | 0.868 |         0.664 |   1.072 |               0.920 |             0.013 |             0.011 |
| sealed_headline | all   | NW    | own     | FULL    | 14152 |    3239 | 0.893 |         0.664 |   1.058 |               0.830 |             0.027 |             0.004 |
| sealed_headline | all   | NW    | own     | MF_ONLY | 14152 |    3239 | 0.706 |         0.664 |   1.285 |               1.585 |             0.004 |             0.001 |
| sealed_headline | all   | NW    | nbr3    | WAVE    | 14152 |    3279 | 0.871 |         0.662 |   1.098 |               0.899 |             0.014 |             0.014 |
| sealed_headline | all   | NW    | nbr3    | FULL    | 14152 |    3279 | 0.896 |         0.662 |   1.057 |               0.809 |             0.031 |             0.006 |
| sealed_headline | all   | NW    | nbr3    | MF_ONLY | 14152 |    3279 | 0.712 |         0.662 |   1.363 |               1.433 |             0.004 |             0.002 |
| sealed_headline | all   | NW    | any     | WAVE    | 14152 |    4815 | 0.792 |         0.612 |   1.059 |               0.898 |             0.069 |             0.052 |
| sealed_headline | all   | NW    | any     | FULL    | 14152 |    4815 | 0.819 |         0.612 |   1.019 |               0.865 |             0.102 |             0.033 |
| sealed_headline | all   | NW    | any     | MF_ONLY | 14152 |    4815 | 0.646 |         0.612 |   0.991 |               1.146 |             0.029 |             0.035 |
| sealed_headline | all   | SW    | own     | WAVE    |   659 |     188 | 0.818 |         0.716 |   2.230 |               2.435 |             0.001 |             0.000 |
| sealed_headline | all   | SW    | own     | FULL    |   659 |     188 | 0.908 |         0.716 |   1.600 |               1.241 |             0.077 |             0.002 |
| sealed_headline | all   | SW    | own     | MF_ONLY |   659 |     188 | 0.773 |         0.716 | nan     |               2.587 |             0.010 |             0.000 |
| sealed_headline | all   | SW    | nbr3    | WAVE    |   659 |     199 | 0.853 |         0.715 |   2.227 |               2.114 |             0.005 |             0.000 |
| sealed_headline | all   | SW    | nbr3    | FULL    |   659 |     199 | 0.891 |         0.715 |   1.495 |               1.283 |             0.064 |             0.004 |
| sealed_headline | all   | SW    | nbr3    | MF_ONLY |   659 |     199 | 0.768 |         0.715 |   2.201 |               2.600 |             0.012 |             0.000 |
| sealed_headline | all   | SW    | any     | WAVE    |   659 |     252 | 0.803 |         0.673 |   2.405 |               1.522 |             0.012 |             0.006 |
| sealed_headline | all   | SW    | any     | FULL    |   659 |     252 | 0.845 |         0.673 |   1.989 |               1.385 |             0.061 |             0.006 |
| sealed_headline | all   | SW    | any     | MF_ONLY |   659 |     252 | 0.641 |         0.673 |   1.042 |               1.257 |             0.037 |             0.019 |
| sealed_semi     | E1    | NW    | own     | WAVE    | 18613 |    5779 | 0.765 |         0.681 |   0.928 |               1.401 |             0.016 |             0.019 |
| sealed_semi     | E1    | NW    | own     | FULL    | 18613 |    5779 | 0.802 |         0.681 |   0.902 |               1.239 |             0.040 |             0.016 |
| sealed_semi     | E1    | NW    | own     | MF_ONLY | 18613 |    5779 | 0.733 |         0.681 |   1.398 |               2.008 |             0.005 |             0.007 |
| sealed_semi     | E1    | NW    | nbr3    | WAVE    | 18613 |    5799 | 0.762 |         0.680 |   0.935 |               1.301 |             0.021 |             0.023 |
| sealed_semi     | E1    | NW    | nbr3    | FULL    | 18613 |    5799 | 0.803 |         0.680 |   0.906 |               1.216 |             0.038 |             0.015 |
| sealed_semi     | E1    | NW    | nbr3    | MF_ONLY | 18613 |    5799 | 0.740 |         0.680 |   1.411 |               1.801 |             0.008 |             0.011 |
| sealed_semi     | E1    | NW    | any     | WAVE    | 18613 |    8233 | 0.732 |         0.623 |   0.896 |               1.151 |             0.125 |             0.144 |
| sealed_semi     | E1    | NW    | any     | FULL    | 18613 |    8233 | 0.773 |         0.623 |   0.889 |               1.144 |             0.134 |             0.134 |

- W9 member probability: beta 0.25 from 7665 truth groups (`post/member_beta.json`).

## W10 false share by cut-off

| region   | label   | day_class   |   cut |    data |   rev_control |   scrambled_control |   false_share_rev |   false_share_scrambled |   families_covered |
|:---------|:--------|:------------|------:|--------:|--------------:|--------------------:|------------------:|------------------------:|-------------------:|
| ArmB     | any     | all_days    | 0.500 |  146829 |          8168 |              140182 |             0.056 |                   0.955 |                 41 |
| ArmB     | any     | all_days    | 0.750 |   11484 |           183 |                9240 |             0.016 |                   0.805 |                 41 |
| ArmB     | any     | all_days    | 0.900 |     852 |             1 |                 371 |             0.001 |                   0.435 |                 41 |
| ArmB     | any     | quiet       | 0.500 |  125524 |          7485 |              121815 |             0.060 |                   0.970 |                 41 |
| ArmB     | any     | quiet       | 0.750 |    9307 |           167 |                8489 |             0.018 |                   0.912 |                 33 |
| ArmB     | any     | quiet       | 0.900 |     434 |             1 |                 348 |             0.002 |                   0.802 |                 24 |
| ArmB     | any     | tremor30    | 0.500 |   14025 |           510 |               11554 |             0.036 |                   0.824 |                 41 |
| ArmB     | any     | tremor30    | 0.750 |    1705 |            11 |                 478 |             0.006 |                   0.280 |                 41 |
| ArmB     | any     | tremor30    | 0.900 |     362 |             0 |                  18 |             0.000 |                   0.050 |                 40 |
| ArmB     | any     | tremor60    | 0.500 |   21305 |           683 |               18367 |             0.032 |                   0.862 |                 41 |
| ArmB     | any     | tremor60    | 0.750 |    2177 |            16 |                 751 |             0.007 |                   0.345 |                 41 |
| ArmB     | any     | tremor60    | 0.900 |     418 |             0 |                  23 |             0.000 |                   0.055 |                 41 |
| ArmB     | nbr3    | all_days    | 0.500 |   31095 |           717 |               26377 |             0.023 |                   0.848 |                 41 |
| ArmB     | nbr3    | all_days    | 0.750 |    3125 |            15 |                1952 |             0.005 |                   0.625 |                 41 |
| ArmB     | nbr3    | all_days    | 0.900 |     322 |             0 |                 113 |             0.000 |                   0.351 |                 37 |
| ArmB     | nbr3    | quiet       | 0.500 |   26052 |           644 |               23480 |             0.025 |                   0.901 |                 37 |
| ArmB     | nbr3    | quiet       | 0.750 |    2190 |            15 |                1815 |             0.007 |                   0.829 |                 30 |
| ArmB     | nbr3    | quiet       | 0.900 |     122 |             0 |                 110 |             0.000 |                   0.902 |                 13 |
| ArmB     | nbr3    | tremor30    | 0.500 |    3749 |            57 |                1941 |             0.015 |                   0.518 |                 41 |
| ArmB     | nbr3    | tremor30    | 0.750 |     796 |             0 |                  96 |             0.000 |                   0.121 |                 41 |
| ArmB     | nbr3    | tremor30    | 0.900 |     176 |             0 |                   3 |             0.000 |                   0.017 |                 36 |
| ArmB     | nbr3    | tremor60    | 0.500 |    5043 |            73 |                2897 |             0.014 |                   0.574 |                 41 |
| ArmB     | nbr3    | tremor60    | 0.750 |     935 |             0 |                 137 |             0.000 |                   0.147 |                 41 |
| ArmB     | nbr3    | tremor60    | 0.900 |     200 |             0 |                   3 |             0.000 |                   0.015 |                 37 |
| ArmB     | own     | all_days    | 0.500 |   11856 |           412 |                8391 |             0.035 |                   0.708 |                 41 |
| ArmB     | own     | all_days    | 0.750 |    1206 |             4 |                 494 |             0.003 |                   0.410 |                 40 |
| ArmB     | own     | all_days    | 0.900 |      60 |             0 |                   4 |             0.000 |                   0.067 |                 17 |
| ArmB     | own     | quiet       | 0.500 |    9307 |           358 |                7512 |             0.038 |                   0.807 |                 36 |
| ArmB     | own     | quiet       | 0.750 |     641 |             4 |                 463 |             0.006 |                   0.722 |                 25 |
| ArmB     | own     | quiet       | 0.900 |       4 |             0 |                   4 |             0.000 |                   1.000 |                  4 |
| ArmB     | own     | tremor30    | 0.500 |    2057 |            38 |                 585 |             0.018 |                   0.284 |                 41 |
| ArmB     | own     | tremor30    | 0.750 |     491 |             0 |                  26 |             0.000 |                   0.053 |                 40 |
| ArmB     | own     | tremor30    | 0.900 |      50 |             0 |                   0 |             0.000 |                   0.000 |                 14 |
| ArmB     | own     | tremor60    | 0.500 |    2549 |            54 |                 879 |             0.021 |                   0.345 |                 41 |
| ArmB     | own     | tremor60    | 0.750 |     565 |             0 |                  31 |             0.000 |                   0.055 |                 40 |
| ArmB     | own     | tremor60    | 0.900 |      56 |             0 |                   0 |             0.000 |                   0.000 |                 16 |
| NW       | any     | all_days    | 0.500 |  673507 |        163874 |              559019 |             0.243 |                   0.830 |                132 |
| NW       | any     | all_days    | 0.750 |   76619 |         11008 |               16623 |             0.144 |                   0.217 |                132 |
| NW       | any     | all_days    | 0.900 |   19563 |           357 |                 154 |             0.018 |                   0.008 |                131 |
| NW       | any     | quiet       | 0.500 |  457709 |        142899 |              434874 |             0.312 |                   0.950 |                132 |
| NW       | any     | quiet       | 0.750 |   19377 |          9633 |               12779 |             0.497 |                   0.659 |                128 |
| NW       | any     | quiet       | 0.900 |    1311 |           320 |                  81 |             0.244 |                   0.062 |                 98 |
| NW       | any     | tremor30    | 0.500 |  173641 |         13874 |               90718 |             0.080 |                   0.522 |                132 |
| NW       | any     | tremor30    | 0.750 |   51719 |           935 |                2939 |             0.018 |                   0.057 |                132 |
| NW       | any     | tremor30    | 0.900 |   16569 |            26 |                  61 |             0.002 |                   0.004 |                128 |
| NW       | any     | tremor60    | 0.500 |  215798 |         20975 |              124145 |             0.097 |                   0.575 |                132 |
| NW       | any     | tremor60    | 0.750 |   57242 |          1375 |                3844 |             0.024 |                   0.067 |                132 |
| NW       | any     | tremor60    | 0.900 |   18252 |            37 |                  73 |             0.002 |                   0.004 |                128 |
| NW       | nbr3    | all_days    | 0.500 |  150290 |         27066 |               69186 |             0.180 |                   0.460 |                132 |
| NW       | nbr3    | all_days    | 0.750 |   38830 |          1598 |                1900 |             0.041 |                   0.049 |                132 |
| NW       | nbr3    | all_days    | 0.900 |   11085 |            45 |                  23 |             0.004 |                   0.002 |                123 |
| NW       | nbr3    | quiet       | 0.500 |   63601 |         23918 |               54371 |             0.376 |                   0.855 |                132 |
| NW       | nbr3    | quiet       | 0.750 |    4213 |          1421 |                1352 |             0.337 |                   0.321 |                115 |
| NW       | nbr3    | quiet       | 0.900 |     578 |            40 |                  13 |             0.069 |                   0.022 |                 73 |
| NW       | nbr3    | tremor30    | 0.500 |   76920 |          2133 |               11284 |             0.028 |                   0.147 |                132 |
| NW       | nbr3    | tremor30    | 0.750 |   31561 |           121 |                 433 |             0.004 |                   0.014 |                132 |
| NW       | nbr3    | tremor30    | 0.900 |    9449 |             3 |                  10 |             0.000 |                   0.001 |                122 |
| NW       | nbr3    | tremor60    | 0.500 |   86689 |          3148 |               14815 |             0.036 |                   0.171 |                132 |
| NW       | nbr3    | tremor60    | 0.750 |   34617 |           177 |                 548 |             0.005 |                   0.016 |                132 |
| NW       | nbr3    | tremor60    | 0.900 |   10507 |             5 |                  10 |             0.000 |                   0.001 |                122 |
| NW       | own     | all_days    | 0.500 |  102544 |         20951 |               28943 |             0.204 |                   0.282 |                132 |
| NW       | own     | all_days    | 0.750 |   32561 |          1142 |                 945 |             0.035 |                   0.029 |                132 |
| NW       | own     | all_days    | 0.900 |    7811 |            20 |                  11 |             0.003 |                   0.001 |                109 |
| NW       | own     | quiet       | 0.500 |   31163 |         18370 |               22666 |             0.589 |                   0.727 |                130 |
| NW       | own     | quiet       | 0.750 |    2956 |          1012 |                 634 |             0.342 |                   0.214 |                109 |
| NW       | own     | quiet       | 0.900 |     345 |            19 |                   6 |             0.055 |                   0.017 |                 58 |
| NW       | own     | tremor30    | 0.500 |   64279 |          1807 |                4970 |             0.028 |                   0.077 |                132 |
| NW       | own     | tremor30    | 0.750 |   26866 |            89 |                 257 |             0.003 |                   0.010 |                131 |
| NW       | own     | tremor30    | 0.900 |    6597 |             0 |                   5 |             0.000 |                   0.001 |                107 |
| NW       | own     | tremor60    | 0.500 |   71381 |          2581 |                6277 |             0.036 |                   0.088 |                132 |
| NW       | own     | tremor60    | 0.750 |   29605 |           130 |                 311 |             0.004 |                   0.011 |                131 |
| NW       | own     | tremor60    | 0.900 |    7466 |             1 |                   5 |             0.000 |                   0.001 |                107 |
| OlyA     | any     | all_days    | 0.500 | 1683802 |        134748 |             1654422 |             0.080 |                   0.983 |                 28 |
| OlyA     | any     | all_days    | 0.750 |   14000 |          2666 |               10439 |             0.190 |                   0.746 |                 28 |
| OlyA     | any     | all_days    | 0.900 |     768 |            20 |                  45 |             0.026 |                   0.059 |                 23 |
| OlyA     | any     | quiet       | 0.500 | 1442589 |        120943 |             1419429 |             0.084 |                   0.984 |                 28 |
| OlyA     | any     | quiet       | 0.750 |   11118 |          2453 |                9243 |             0.221 |                   0.831 |                 27 |
| OlyA     | any     | quiet       | 0.900 |     100 |            19 |                  36 |             0.190 |                   0.360 |                 19 |
| OlyA     | any     | tremor30    | 0.500 |  149273 |          8775 |              144722 |             0.059 |                   0.970 |                 28 |
| OlyA     | any     | tremor30    | 0.750 |    2305 |           147 |                 763 |             0.064 |                   0.331 |                 27 |
| OlyA     | any     | tremor30    | 0.900 |     652 |             1 |                   7 |             0.002 |                   0.011 |                 15 |
| OlyA     | any     | tremor60    | 0.500 |  241213 |         13805 |              234993 |             0.057 |                   0.974 |                 28 |
| OlyA     | any     | tremor60    | 0.750 |    2882 |           213 |                1196 |             0.074 |                   0.415 |                 28 |
| OlyA     | any     | tremor60    | 0.900 |     668 |             1 |                   9 |             0.001 |                   0.013 |                 19 |
| OlyA     | nbr3    | all_days    | 0.500 |  194713 |         16800 |              179768 |             0.086 |                   0.923 |                 28 |
| OlyA     | nbr3    | all_days    | 0.750 |    3119 |           403 |                1593 |             0.129 |                   0.511 |                 28 |
| OlyA     | nbr3    | all_days    | 0.900 |     391 |             2 |                  19 |             0.005 |                   0.049 |                 19 |
| OlyA     | nbr3    | quiet       | 0.500 |  167194 |         15243 |              156141 |             0.091 |                   0.934 |                 28 |
| OlyA     | nbr3    | quiet       | 0.750 |    1813 |           370 |                1402 |             0.204 |                   0.773 |                 26 |
| OlyA     | nbr3    | quiet       | 0.900 |      44 |             2 |                  16 |             0.045 |                   0.364 |                 15 |
| OlyA     | nbr3    | tremor30    | 0.500 |   17492 |          1039 |               14519 |             0.059 |                   0.830 |                 28 |
| OlyA     | nbr3    | tremor30    | 0.750 |    1168 |            21 |                 117 |             0.018 |                   0.100 |                 23 |
| OlyA     | nbr3    | tremor30    | 0.900 |     342 |             0 |                   2 |             0.000 |                   0.006 |                 12 |
| OlyA     | nbr3    | tremor60    | 0.500 |   27519 |          1557 |               23627 |             0.057 |                   0.859 |                 28 |
| OlyA     | nbr3    | tremor60    | 0.750 |    1306 |            33 |                 191 |             0.025 |                   0.146 |                 26 |
| OlyA     | nbr3    | tremor60    | 0.900 |     347 |             0 |                   3 |             0.000 |                   0.009 |                 13 |
| OlyA     | own     | all_days    | 0.500 |   26448 |          6601 |               20874 |             0.250 |                   0.789 |                 28 |
| OlyA     | own     | all_days    | 0.750 |    1518 |           159 |                 410 |             0.105 |                   0.270 |                 28 |
| OlyA     | own     | all_days    | 0.900 |     255 |             1 |                   4 |             0.004 |                   0.016 |                 13 |
| OlyA     | own     | quiet       | 0.500 |   21907 |          6059 |               18544 |             0.277 |                   0.846 |                 27 |
| OlyA     | own     | quiet       | 0.750 |     541 |           151 |                 368 |             0.279 |                   0.680 |                 24 |
| OlyA     | own     | quiet       | 0.900 |      24 |             1 |                   3 |             0.042 |                   0.125 |                 10 |
| OlyA     | own     | tremor30    | 0.500 |    3475 |           398 |                1505 |             0.115 |                   0.433 |                 27 |
| OlyA     | own     | tremor30    | 0.750 |     922 |             6 |                  28 |             0.007 |                   0.030 |                 22 |
| OlyA     | own     | tremor30    | 0.900 |     227 |             0 |                   1 |             0.000 |                   0.004 |                 10 |
| OlyA     | own     | tremor60    | 0.500 |    4541 |           542 |                2330 |             0.119 |                   0.513 |                 28 |
| OlyA     | own     | tremor60    | 0.750 |     977 |             8 |                  42 |             0.008 |                   0.043 |                 26 |
| OlyA     | own     | tremor60    | 0.900 |     231 |             0 |                   1 |             0.000 |                   0.004 |                 10 |
| SW       | any     | all_days    | 0.500 |    3012 |          1893 |                 983 |             0.628 |                   0.326 |                 73 |
| SW       | any     | all_days    | 0.750 |     603 |            42 |                  19 |             0.070 |                   0.032 |                 72 |
| SW       | any     | all_days    | 0.900 |      92 |             0 |                   0 |             0.000 |                   0.000 |                 37 |
| SW       | any     | quiet       | 0.500 |     745 |          1645 |                 687 |             2.208 |                   0.922 |                 69 |
| SW       | any     | quiet       | 0.750 |      12 |            38 |                  13 |             3.167 |                   1.083 |                 12 |
| SW       | any     | quiet       | 0.900 |       0 |             0 |                   0 |           nan     |                 nan     |                  0 |
| SW       | any     | tremor30    | 0.500 |    2233 |           182 |                 263 |             0.082 |                   0.118 |                 73 |
| SW       | any     | tremor30    | 0.750 |     589 |             2 |                   6 |             0.003 |                   0.010 |                 71 |
| SW       | any     | tremor30    | 0.900 |      92 |             0 |                   0 |             0.000 |                   0.000 |                 37 |
| SW       | any     | tremor60    | 0.500 |    2267 |           248 |                 296 |             0.109 |                   0.131 |                 73 |
| SW       | any     | tremor60    | 0.750 |     591 |             4 |                   6 |             0.007 |                   0.010 |                 71 |
| SW       | any     | tremor60    | 0.900 |      92 |             0 |                   0 |             0.000 |                   0.000 |                 37 |

## W11 rates

- k factors: 903 (unit, era set) rows; median k by era set and region:

| net   | unit_kind   | reg   |   n |   k_med |   k_raw_med |
|:------|:------------|:------|----:|--------:|------------:|
| E1    | family      | ArmB  |  43 | 161.657 |     163.250 |
| E1    | family      | NW    | 132 |  50.977 |      50.977 |
| E1    | family      | SW    |  79 | nan     |     nan     |
| E1    | olyA_group  | OlyA  |  47 |  82.500 |      82.500 |
| E2    | family      | ArmB  |  43 | 115.860 |     115.875 |
| E2    | family      | NW    | 132 |  39.200 |      39.200 |
| E2    | family      | SW    |  79 | 216.345 |     219.000 |
| E2    | olyA_group  | OlyA  |  47 |  86.646 |      86.750 |
| E3    | family      | ArmB  |  43 | 186.000 |     186.000 |
| E3    | family      | NW    | 132 |  54.000 |      54.000 |
| E3    | family      | SW    |  79 | 352.681 |     353.000 |
| E3    | olyA_group  | OlyA  |  47 | 440.365 |    1111.000 |

- Rate tables `rates/rates_daily.parquet`, `rates/rates_monthly.csv`:

| unit_kind   | label   |   units |   months |   net_count |   rate_scaled |
|:------------|:--------|--------:|---------:|------------:|--------------:|
| family      | nbr3    |     246 |      237 |   42952.000 |         0.000 |
| family      | own     |     247 |      237 |   34886.000 |    954827.005 |
| olyA_group  | nbr3    |      10 |      237 |    2742.000 |    628106.880 |

- Step checks: 1624 of 11552 flagged (|z| > 2 vs the same calendar windows of other years); a step next to an ETS cannot be told from activity (`rates/step_checks.csv`).

| region   | date       | netsta   | change   |   window_days |   d_quiet_rate |   z_quiet |   d_rate_per_tremor_hour |   z_tremor_hour | flag   |
|:---------|:-----------|:---------|:---------|--------------:|---------------:|----------:|-------------------------:|----------------:|:-------|
| ArmB     | 2007-03-12 | UW.SQM   | joins    |            21 |          0.571 |    10.037 |                      nan |             nan | True   |
| ArmB     | 2007-03-12 | UW.SQM   | joins    |            45 |          0.233 |     2.380 |                      nan |             nan | True   |
| ArmB     | 2007-03-13 | PB.B003  | joins    |            21 |          0.762 |    12.729 |                      nan |             nan | True   |
| ArmB     | 2007-03-13 | PB.B003  | joins    |            45 |          0.238 |     2.408 |                      nan |             nan | True   |
| ArmB     | 2007-03-13 | PB.B004  | joins    |            21 |          0.762 |    12.729 |                      nan |             nan | True   |
| ArmB     | 2007-03-13 | PB.B004  | joins    |            45 |          0.238 |     2.408 |                      nan |             nan | True   |
| ArmB     | 2007-03-26 | PB.B003  | joins    |            21 |          0.762 |     2.805 |                      nan |             nan | True   |
| ArmB     | 2007-04-09 | UW.HDW   | joins    |            21 |         -0.810 |    -8.664 |                      nan |             nan | True   |
| ArmB     | 2007-04-12 | PB.B001  | joins    |            21 |         -0.857 |   -10.522 |                      nan |             nan | True   |
| ArmB     | 2007-04-12 | PB.B006  | joins    |            21 |         -0.857 |   -10.522 |                      nan |             nan | True   |
| ArmB     | 2007-04-12 | PB.B003  | joins    |            21 |         -0.857 |   -10.522 |                      nan |             nan | True   |
| ArmB     | 2007-04-12 | PB.B018  | joins    |            21 |         -0.857 |   -10.522 |                      nan |             nan | True   |
| ArmB     | 2007-04-12 | PB.B004  | joins    |            21 |         -0.857 |   -10.522 |                      nan |             nan | True   |
| ArmB     | 2007-04-12 | PB.B005  | joins    |            21 |         -0.857 |   -10.522 |                      nan |             nan | True   |
| ArmB     | 2007-04-12 | PB.B007  | joins    |            21 |         -0.857 |   -10.522 |                      nan |             nan | True   |
| ArmB     | 2007-04-13 | UW.SQM   | joins    |            21 |         -0.905 |   -11.121 |                      nan |             nan | True   |
| ArmB     | 2007-04-19 | UW.GNW   | joins    |            21 |         -0.714 |    -9.321 |                      nan |             nan | True   |
| ArmB     | 2007-05-11 | PB.B007  | leaves   |            21 |          0.476 |     2.977 |                      nan |             nan | True   |
| ArmB     | 2007-05-11 | PB.B007  | leaves   |            45 |          0.333 |     2.181 |                      nan |             nan | True   |
| ArmB     | 2007-05-17 | UW.GNW   | joins    |            21 |          0.905 |     3.871 |                      nan |             nan | True   |
| ArmB     | 2007-05-17 | UW.GNW   | joins    |            45 |          0.644 |     5.265 |                      nan |             nan | True   |
| ArmB     | 2007-05-18 | PB.B007  | joins    |            21 |          0.857 |     3.808 |                      nan |             nan | True   |
| ArmB     | 2007-05-18 | PB.B007  | joins    |            45 |          0.622 |     5.259 |                      nan |             nan | True   |
| ArmB     | 2007-05-20 | UW.SQM   | joins    |            21 |          0.762 |     3.159 |                      nan |             nan | True   |
| ArmB     | 2007-05-20 | UW.SQM   | joins    |            45 |          0.578 |     4.764 |                      nan |             nan | True   |
| ArmB     | 2007-05-20 | PB.B013  | joins    |            21 |          0.762 |     3.159 |                      nan |             nan | True   |
| ArmB     | 2007-05-20 | PB.B013  | joins    |            45 |          0.578 |     4.764 |                      nan |             nan | True   |
| ArmB     | 2007-06-21 | PB.B004  | leaves   |            21 |         -0.714 |    -4.138 |                      nan |             nan | True   |
| ArmB     | 2007-06-21 | PB.B004  | leaves   |            45 |         -0.444 |    -3.111 |                      nan |             nan | True   |
| ArmB     | 2007-06-28 | UW.GNW   | joins    |            21 |         -0.476 |    -5.247 |                      nan |             nan | True   |

## W12 windows

- pilot: {'verify': True, 'extrapolation': {'pilot_month': '2015-10', 'cut_GB': 10.03, 'raw_GB': 23.248, 'cut_over_raw': 0.43143496214728144, 'raw_total_TB': 8.415315553792, 'extrapolated_TB': 3.6306613474076803, 'n_windows': 31308568, 'minutes': 14.8, 'limit_TB': 3.0, 'note': 'the pilot month is the busiest production month; full-store estimate = pilot cut/raw ratio x all raw in the roots', 'pass': False}}
- verify: None
- stack check: None

## W13 delete-candidate list (nothing deleted)

```
WA PRODUCTION delete-CANDIDATE list (NOTHING DELETED; joint SVI + Washington decision by the user):
  WA candidates 1363 entries, 7.759 TB by kind {'cache_svi': 1399.0, 'cache_svi_r1': 2.8, 'cache_wa_prod': 5912.6, 'cache_wa_trial': 444.4} GB
  SVI R9 candidates 9025 files, 0.193 TB (lfe_extension/results/refine/r9/)
  held 5.386 TB (+ other projects' raw, listed only); reasons {'cut not verified': 4.428, "lfe_array_discovery/traces: listed only -- that project's decision": 0.0, "lfe_catalog_repair/traces: listed only -- that project's decision": 0.0, "lfe_picker/traces: listed only -- that project's decision": 0.0, "not an SVI station held for WA (SVI's own decision, R9)": 0.1, "pre-2014 SVI raw: SVI's truth years / not part of the joint list": 0.043, 'symlink (frees nothing; the target belongs to another project)': 0.0, 'truth-era cache: kept': 0.814}
  cut verified: False; symlink readers found: 38379
```

## Deviations from the plan (with reasons)

Written at launch (2026-10-03), before any production result existed. Each item: what differs from WA_PROD_PLAN.md + the
user's decisions of 2026-10-03, and why.

1. **D1, A2011 leakage: rebuild instead of flag.** The step-4b templates (`tpl_step4b_a`) were built from A2011 array
   events. They are replaced by `v10b` = the same frozen recipe on the V2010 BUILD-half array detections at the same nodes
   (7,636 templates, 41 nodes, 200 families; step4b_a had 9,246). A2011 is then out-of-sample for every station (plan D1
   recommendation). The V2010 build half stays in-sample.
2. **Fixed rate networks by a rule written before any run** (`setup/core_sets.json`, `core_set_table.csv`): permanent
   station, raw on >= 80 % of the era's days, a W4a template source (own node with >= 30 families, or a W2-approved
   carry-over) on >= 95 % of its raw days, LRIV out after 2019-05-30 (plan D10 / SVI policy), never XU.W060. Headline E3
   2015-2026 = 19 stations (B001, B003-B007, B009-B011, B013, B014, B017, B018, B019, B926, B941, B943, GMW, GNW), not the
   21 of the plan: PGC, VGZ, DOSE and HDW have no usable W4a template for part of 2015-2026 (their control-level epochs are
   re-templated, user decision), so they cannot be in a FIXED set; B019, B941 and GMW are now on disk (tiers B/C). E1
   2007-2011 = 12 stations, E2 2012-2014 = 37.
3. **Re-templated epochs and other new nodes enter only the rich pass** (W4b after the confidence model; T2B rule P_own_WAVE
   >= 0.3 from the network without the node, families <= 50 km, the node's first two WA ETS else its first 24 months, with
   the DUC_FIN per-channel scrambled-control check: N >= 20, >= 10 day-shifted controls, coherence above all). Node changes
   of non-core stations have no W2 swap check, so they are re-templated rather than carried over. The sealed months are
   never used to build a template.
4. **B014 2023 check** is the W2 monthly stack-test record before / after its 2023-08 .. 2024-08 outage (no metadata change
   exists): ok share 0.30 before, 0.50 after -> no change detected, old template kept (rule in `setup/b014_check.json`).
5. **Health exclusions**: only W2 dead months and dead days (969 station-days). W2 'flat' windows (>= 1 % repeated samples)
   are NOT exclusions: they flag years of normal PB borehole data. Flat / dead channels are caught per window by the
   frozen energy-floor guard and the drop-dead rule (> 5 % of the day guarded -> station absent that day). Non-core
   stations have no W2 health table; the same per-day guard applies to them.
6. **Scrambled-control shifts** are deterministic per (date, site) (the SVI_REF fix) instead of wa_mf's sequential draws,
   so a site's shift does not depend on which other sites are present. Data rows reproduce wa_mf.py exactly (smoke: 10,698 /
   10,698 detections of the go/no-go core-21 run on 2010-08-20, max |d net_mad| = 0); control rows differ by construction.
7. **Truth runs**: networks E3 + all with 2 scrambled controls + the time-reversed control (`truth`), E1 + E2 with 1 scrambled
   + time-reversed (`truth_e12`) -- GPU memory. XU / XG station-days before 2012 are never used (they made the truth).
8. **W8 PNET gate** uses the frozen fine-tune `ft1.pt`, PNET block on the `truth` run (E3 + all). Rule as the user fixed it
   (dAUC own >= 0.01 with day-bootstrap CI > 0, slope 0.8-1.2); the plan's extra condition on the never-trained stations
   (B017, B018, GNW, HDW) is reported as a column, not applied.
9. **Fit rows** are a uniform random subsample (<= 600,000 data rows of the fit days; base rate kept) to bound compute.
10. **Refit risk 'older model'** = FULL features fitted on the V2008 + V2010 test halves only (the truth the earlier WA
    confidence used); the earlier model itself was never saved as a file.
11. **W9 member probability** is a simplified column: softmax(beta x net_mad) over the <= 3-km neighbours' candidates within
    +-2 s (rows with P_nbr3_WAVE >= 0.3), beta by maximum likelihood on the truth fit rows -- not the step-2 station-level
    attribution, which needs a separate group matched filter.
12. **W11**: k factors for the headline labels (own per family; nbr3 per Olympic A position group); the nbr3 count of
    NW / SW / Arm B families is a column without k. The daily-stack activity column (stack_activity.py, a separate GPU
    engine) and the SVI withheld-station inter-ETS validation were NOT run: inter-ETS rates are a control-subtracted bracket.
13. **W12 windows**: S times at non-template stations from a straight-ray harmonic-mean S speed through the USGS CVM v1.7
    1-D profile (royer_NW_WA / royer_SW_WA) + the family's median template residual. 'Every recording station' is limited
    to stations <= 150 km from a Washington family (the California stations of lfe_south/traces are otherwise in the roots);
    stations with no known coordinates are listed in `cut_roots/skipped_stations.csv`. Pre-registered size limit: the
    pilot month extrapolated to all raw must be < 3.0 TB, else the full cut stops for a user decision.
14. **Span** stays 2007-01-01 .. 2026-09-15 (the W1 / W2 span), although tier A now holds data to 2026-09-30.
15. **Rich pass network** = every station with a W4a or W4b template; it runs only on days with a W4b node live.

16. (main session, 2026-10-04) The first fit crashed: fam_code (264 families) exceeds the 255-category limit of
    HistGradientBoosting. fam_code is now a plain numeric feature (only in the FULL / NO_DENSITY sensitivity models;
    the headline WAVE model never uses it); reg_code stays categorical. Queue relaunched; finished stages resume.
16b/main 10-04 09:10Z: production MF stalled 05:00Z→09:10Z waiting for cache/done_p2_2013-09; that month has no phase-2 (traces_dl2) jobs, so the cache stage never wrote its marker. Markers written for all months with 0 phase-2 jobs. Code fix for future runs: write done_p2 for empty months.
17. (main session, 2026-10-06) w4b_targets crashed 10-05 19:43Z: cache manifest_p1_2013-09.csv has only 'empty' rows and was written without a 'node' column. node_days() now skips such files (backup .bak_pre_nodecol); relaunched, finished stages resume. ~11 h idle before it was noticed.
18. (main session, 2026-10-06) User decisions after the cut pilot and the quiet-day waveform test: window floor P_max >= 0.2 (conf/floor_user.json overrides the automatic 0.1; pilot at 0.1 kept in cut_pilot_floor01/); NW quiet-day cut-off 0.85; days without a tremor catalogue labelled 'unknown'. Queue relaunched for pilot -> cut -> verify -> stack check -> delete list -> report.

<!-- post_v2 section -->
## Post-processing v2 (user decisions 2026-10-06)

Full results: [`wa_prod/post_v2/POST_V2_RESULTS.md`](wa_prod/post_v2/POST_V2_RESULTS.md). Day classes per family: 'tremor' = PNSN tremor within 30 km (before 2009-08-06: inside the WECC ETS windows); 'quiet' = the PNSN catalogue covers the day (2009-08-06 .. 2026-05-27) and shows no tremor within 30 km; 'unknown' = every other day. Headline selection: P_own_WAVE_PNET >= 0.75 for NW and Arm B on every day class; SW P_own_WAVE_PNET >= 0.5 on tremor days only; Olympic A position-group rates on P_nbr3_WAVE >= 0.75. The WAVE-model columns and the earlier class rules (tremor 0.75 / SW 0.5, NW quiet 0.85) are kept as secondary columns.

- Headline 2012-2026 (core network; data / time-reversed control / scrambled control; `post_v2/summary_headline_false_share_2012_2026.csv`): ArmB quiet P >= 0.75: 37 / 0 / 7; ArmB tremor P >= 0.75: 276 / 0 / 1; ArmB unknown P >= 0.75: 1 / 0 / 0; NW quiet P >= 0.75: 2,887 / 22 / 16; NW tremor P >= 0.75: 21,637 / 4 / 118; NW unknown P >= 0.75: 312 / 0 / 0; OlyA quiet P >= 0.75: 959 / 36 / 781; OlyA tremor P >= 0.75: 767 / 3 / 80; OlyA unknown P >= 0.75: 11 / 0 / 6; SW tremor P >= 0.5: 770 / 31 / 167.
- Files: `post_v2/catalogue_v2.parquet` (every core-network candidate; `selected_v2` = headline, `selected_v2_wave_classrules` / `selected_v2_pnet_classrules` = secondary), `post_v2/rates_*`.

