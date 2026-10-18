# PNET as independent per-repeat evidence — WA extension (Merlin step 6)

Written by `scripts/pnet_report.py`; queue `scripts/pnet_queue.sh`; outputs `results/pnet/`. Numbers only; interpretation after review.

## Summary (plain language)

PNET is a neural-network S-wave picker (PhaseNet-style) trained on Bostock's southern Vancouver Island LFEs. We asked whether its S-wave probability at the times a Washington family's template says the S should arrive gives independent evidence that a matched-filter repeat is real, and whether it can make candidate times on its own. The scrambled control is the same matched filter with each site's data shifted to a different time of the same day (so nothing lines up across sites except by chance).

- **Off-the-shelf PNET (A)** — AUC of `pnet_max`, held-out days pooled: data vs scrambled control 0.50 (NW 0.53, SW 0.59, Olympic 0.49); truth-matched vs unmatched 0.68 (NW 0.72). For comparison the step-1 P gives 0.51 / 0.66. On 1-component EHZ stations the SVI model returns ~0 (it was trained on 3-C only).
- **Fine-tuned on WA (B)** — 104022 labelled positive windows; held-out stations (6): median per-station AUC (labelled arrival vs quiet-day random time) 0.71 → 0.80; vs episode-day random time 0.60 → 0.65; training stations on held-out days vs episode 0.54 → 0.69; quiet-day peaks > 0.5 per hour 29.0 → 3.3. Pre-set rule 'works': **True**.
- Fine-tuned `pnet_max` AUC data vs scrambled control 0.50 (NW 0.53, SW 0.61); truth-matched vs unmatched 0.77 (NW 0.78).
- **Step-1 refit (C) baseline (step 1)**: calibration slope own V2010 test 0.65, A2011 0.54; P > 0.5 data 0.021 vs scrambled control 0.004 (V2010 test); NW data 0.106 vs control 0.011; expected-count Spearman V2010 0.49, A2011 0.48.
- **Step-1 refit (C) + PNET off-the-shelf**: calibration slope own V2010 test 0.63, A2011 0.54; P > 0.5 data 0.019 vs scrambled control 0.003 (V2010 test); NW data 0.113 vs control 0.012; expected-count Spearman V2010 0.50, A2011 0.47.
- **Step-1 refit (C) + PNET fine-tuned**: calibration slope own V2010 test 0.62, A2011 0.55; P > 0.5 data 0.026 vs scrambled control 0.005 (V2010 test); NW data 0.135 vs control 0.012; expected-count Spearman V2010 0.50, A2011 0.48.
- **PNET candidate times for the withheld-station test (D)** — median recovered-template cc: NW 0.27 (MF) / 0.15 (PNET off-the-shelf) / 0.15 (PNET fine-tuned), scrambled control 0.14; SW 0.15 (MF) / 0.14 (PNET off-the-shelf) / 0.15 (PNET fine-tuned), scrambled control 0.14; Olympic 0.86 (MF) / 0.11 (PNET off-the-shelf) / 0.23 (PNET fine-tuned), scrambled control 0.13.

## Input handling

- Weights: `lfe_picker/results/p1/models/pnet_w05.pt` (P1b selected model: PhaseNet-style U-Net, 1.82 M parameters, 60-s windows at 20 Hz, 3-C Z/N/E, predicted-label weight 0.5, trained 2008–2010 + 2012 on SVI stations).
- Inference exactly as `lfe_picker/scripts/p_scan.cf_day` (step 20 s, central third kept, per-window 3-C peak normalisation, bf16).
- Band/rate: the WA cache is the same frozen 20 Hz, 1–8 Hz Butterworth-4 zero-phase chain the model was trained on — no re-filtering.
- Rotation to ZNE (SVI store was rotated with StationXML azimuths): PB EH1/EH2 rotated (B001 335.5/65.5, B003 239/329, B004 167.5/257.5, B005 90.5/180.5, B006 201/291, B007 36.5/126.5, B009 0/90, B010 228.5/318.5, B011 11/101, B013 31.5/301.5, B014 170/260, B017 206/296, B018 21.5/111.5, B019 305/35, B941 173.5/263.5, B943 301.5/31.5); US.NLWA 4/94; UW.HDW has one horizontal (EH2 at 45°) → Z only; 1-C EHZ stations enter as (Z, 0, 0). Invalid samples (cache mask) set to 0.
- Expected S at a station = candidate origin + the family's template offset at that station (template start = stack envelope peak − 0.5 s ≈ onset); value = max CF within ±0.5 s; scrambled-control candidates read at their site's rolled position (as the MF built them).

## A. Off-the-shelf check (and the same features from the fine-tune)

**off-the-shelf (SVI pnet_w05) — AUC, data vs scrambled control** (reference columns: step-1 `n_z2` and step-1 `P`)

| eval_set      | region   |   n_pos |   n_neg |   pnet_max |   pnet_2nd |   n_sites_pnet03 |   n_sites_pnet05 |   pnet_rel_max |   n_z2 |     P |
|:--------------|:---------|--------:|--------:|-----------:|-----------:|-----------------:|-----------------:|---------------:|-------:|------:|
| fit_in_sample | ALL      |  297796 | 1031035 |      0.516 |      0.520 |            0.523 |            0.518 |          0.510 |  0.567 | 0.648 |
| fit_in_sample | NW       |   55588 |  131546 |      0.561 |      0.581 |            0.588 |            0.569 |          0.543 |  0.714 | 0.761 |
| fit_in_sample | SW       |    7337 |    1397 |      0.629 |      0.678 |            0.688 |            0.629 |          0.617 |  0.924 | 0.977 |
| fit_in_sample | Olympic  |  226905 |  878640 |      0.491 |      0.489 |            0.490 |            0.494 |          0.490 |  0.513 | 0.608 |
| fit_in_sample | Arm B    |    7966 |   19452 |      0.550 |      0.568 |            0.577 |            0.552 |          0.533 |  0.682 | 0.719 |
| V2010_test    | ALL      |  172854 |  675877 |      0.500 |      0.502 |            0.503 |            0.504 |          0.499 |  0.518 | 0.509 |
| V2010_test    | NW       |   21800 |   76010 |      0.535 |      0.546 |            0.550 |            0.539 |          0.529 |  0.567 | 0.578 |
| V2010_test    | SW       |     381 |     671 |      0.621 |      0.657 |            0.673 |            0.620 |          0.605 |  0.552 | 0.766 |
| V2010_test    | Olympic  |  148866 |  591669 |      0.494 |      0.494 |            0.495 |            0.498 |          0.493 |  0.510 | 0.499 |
| V2010_test    | Arm B    |    1807 |    7527 |      0.497 |      0.495 |            0.498 |            0.506 |          0.487 |  0.528 | 0.526 |
| V2008_test    | ALL      |   93754 |  363318 |      0.507 |      0.511 |            0.514 |            0.511 |          0.508 |  0.516 | 0.490 |
| V2008_test    | NW       |   16322 |   54185 |      0.538 |      0.557 |            0.561 |            0.550 |          0.536 |  0.581 | 0.598 |
| V2008_test    | SW       |     305 |     581 |      0.605 |      0.637 |            0.660 |            0.612 |          0.598 |  0.516 | 0.792 |
| V2008_test    | Olympic  |   73842 |  295583 |      0.499 |      0.499 |            0.501 |            0.501 |          0.500 |  0.500 | 0.464 |
| V2008_test    | Arm B    |    3285 |   12969 |      0.490 |      0.497 |            0.497 |            0.495 |          0.488 |  0.513 | 0.511 |
| A2011         | ALL      |  294824 | 1149490 |      0.496 |      0.498 |            0.500 |            0.505 |          0.492 |  0.518 | 0.514 |
| A2011         | NW       |   61999 |  222561 |      0.523 |      0.530 |            0.530 |            0.527 |          0.519 |  0.558 | 0.562 |
| A2011         | SW       |     709 |    2799 |      0.494 |      0.503 |            0.500 |            0.504 |          0.486 |  0.511 | 0.482 |
| A2011         | Olympic  |  215645 |  858740 |      0.486 |      0.487 |            0.490 |            0.497 |          0.482 |  0.507 | 0.501 |
| A2011         | Arm B    |   16471 |   65390 |      0.488 |      0.493 |            0.492 |            0.500 |          0.483 |  0.508 | 0.508 |

**off-the-shelf (SVI pnet_w05) — AUC, own-family truth-matched vs unmatched data** (reference columns: step-1 `n_z2` and step-1 `P`)

| eval_set      | region   |   n_pos |   n_neg |   pnet_max |   pnet_2nd |   n_sites_pnet03 |   n_sites_pnet05 |   pnet_rel_max |   n_z2 |     P |
|:--------------|:---------|--------:|--------:|-----------:|-----------:|-----------------:|-----------------:|---------------:|-------:|------:|
| fit_in_sample | ALL      |   39611 |  120580 |      0.658 |      0.689 |            0.695 |            0.658 |          0.635 |  0.936 | 0.947 |
| fit_in_sample | NW       |   24903 |   18305 |      0.619 |      0.643 |            0.644 |            0.620 |          0.601 |  0.920 | 0.948 |
| fit_in_sample | SW       |    6685 |     168 |      0.647 |      0.632 |            0.633 |            0.612 |          0.623 |  0.948 | 0.989 |
| fit_in_sample | Olympic  |    5181 |  100803 |      0.517 |      0.516 |            0.516 |            0.515 |          0.512 |  0.786 | 0.806 |
| fit_in_sample | Arm B    |    2842 |    1304 |      0.626 |      0.662 |            0.691 |            0.633 |          0.600 |  0.975 | 0.992 |
| V2010_test    | ALL      |    7761 |  107540 |      0.614 |      0.634 |            0.637 |            0.614 |          0.597 |  0.639 | 0.626 |
| V2010_test    | NW       |    3662 |    8865 |      0.657 |      0.689 |            0.682 |            0.663 |          0.632 |  0.677 | 0.710 |
| V2010_test    | Olympic  |    4013 |   97343 |      0.505 |      0.506 |            0.505 |            0.503 |          0.502 |  0.531 | 0.519 |
| V2010_test    | Arm B    |      86 |    1332 |      0.563 |      0.616 |            0.649 |            0.557 |          0.535 |  0.761 | 0.785 |
| V2008_test    | ALL      |    3857 |   10086 |      0.667 |      0.692 |            0.679 |            0.669 |          0.652 |  0.669 | 0.711 |
| V2008_test    | NW       |    3708 |    9939 |      0.667 |      0.690 |            0.675 |            0.666 |          0.652 |  0.665 | 0.705 |
| V2008_test    | SW       |     149 |     147 |      0.682 |      0.759 |            0.736 |            0.724 |          0.667 |  0.558 | 0.819 |
| A2011         | ALL      |    8865 |   53134 |      0.729 |      0.719 |            0.698 |            0.692 |          0.706 |  0.694 | 0.640 |
| A2011         | NW       |    8865 |   53134 |      0.729 |      0.719 |            0.698 |            0.692 |          0.706 |  0.694 | 0.640 |

**off-the-shelf (SVI pnet_w05) — AUC, ≤ 3 km truth-matched vs unmatched data** (reference columns: step-1 `n_z2` and step-1 `P`)

| eval_set      | region   |   n_pos |   n_neg |   pnet_max |   pnet_2nd |   n_sites_pnet03 |   n_sites_pnet05 |   pnet_rel_max |   n_z2 |     P |
|:--------------|:---------|--------:|--------:|-----------:|-----------:|-----------------:|-----------------:|---------------:|-------:|------:|
| fit_in_sample | ALL      |   49027 |  111164 |      0.624 |      0.650 |            0.656 |            0.624 |          0.605 |  0.850 | 0.858 |
| fit_in_sample | NW       |   24945 |   18263 |      0.619 |      0.644 |            0.644 |            0.621 |          0.601 |  0.920 | 0.948 |
| fit_in_sample | SW       |    6685 |     168 |      0.647 |      0.632 |            0.633 |            0.612 |          0.623 |  0.948 | 0.989 |
| fit_in_sample | Olympic  |   14544 |   91440 |      0.506 |      0.510 |            0.515 |            0.506 |          0.502 |  0.609 | 0.596 |
| fit_in_sample | Arm B    |    2853 |    1293 |      0.624 |      0.662 |            0.692 |            0.632 |          0.599 |  0.974 | 0.991 |
| V2010_test    | ALL      |   16425 |   98876 |      0.550 |      0.560 |            0.561 |            0.551 |          0.543 |  0.566 | 0.542 |
| V2010_test    | NW       |    3686 |    8841 |      0.657 |      0.690 |            0.682 |            0.663 |          0.632 |  0.677 | 0.711 |
| V2010_test    | Olympic  |   12647 |   88709 |      0.504 |      0.504 |            0.505 |            0.504 |          0.503 |  0.520 | 0.484 |
| V2010_test    | Arm B    |      92 |    1326 |      0.550 |      0.604 |            0.636 |            0.547 |          0.523 |  0.755 | 0.787 |
| V2008_test    | ALL      |    3892 |   10051 |      0.668 |      0.693 |            0.678 |            0.669 |          0.653 |  0.669 | 0.713 |
| V2008_test    | NW       |    3739 |    9908 |      0.667 |      0.690 |            0.675 |            0.666 |          0.653 |  0.666 | 0.707 |
| V2008_test    | SW       |     153 |     143 |      0.689 |      0.768 |            0.744 |            0.734 |          0.671 |  0.546 | 0.829 |
| A2011         | ALL      |    8959 |   53040 |      0.729 |      0.718 |            0.698 |            0.691 |          0.706 |  0.695 | 0.641 |
| A2011         | NW       |    8959 |   53040 |      0.729 |      0.718 |            0.698 |            0.691 |          0.706 |  0.695 | 0.641 |

**WA fine-tune (ft1) — AUC, data vs scrambled control** (reference columns: step-1 `n_z2` and step-1 `P`)

| eval_set      | region   |   n_pos |   n_neg |   pnet_max |   pnet_2nd |   n_sites_pnet03 |   n_sites_pnet05 |   pnet_rel_max |   n_z2 |     P |
|:--------------|:---------|--------:|--------:|-----------:|-----------:|-----------------:|-----------------:|---------------:|-------:|------:|
| fit_in_sample | ALL      |  297796 | 1031035 |      0.509 |      0.517 |            0.520 |            0.514 |          0.504 |  0.567 | 0.648 |
| fit_in_sample | NW       |   55588 |  131546 |      0.568 |      0.621 |            0.622 |            0.588 |          0.557 |  0.714 | 0.761 |
| fit_in_sample | SW       |    7337 |    1397 |      0.652 |      0.667 |            0.643 |            0.589 |          0.653 |  0.924 | 0.977 |
| fit_in_sample | Olympic  |  226905 |  878640 |      0.496 |      0.485 |            0.490 |            0.497 |          0.494 |  0.513 | 0.608 |
| fit_in_sample | Arm B    |    7966 |   19452 |      0.543 |      0.556 |            0.561 |            0.538 |          0.533 |  0.682 | 0.719 |
| V2010_test    | ALL      |  172854 |  675877 |      0.506 |      0.501 |            0.505 |            0.508 |          0.505 |  0.518 | 0.509 |
| V2010_test    | NW       |   21800 |   76010 |      0.541 |      0.561 |            0.565 |            0.550 |          0.537 |  0.567 | 0.578 |
| V2010_test    | SW       |     381 |     671 |      0.696 |      0.731 |            0.712 |            0.680 |          0.696 |  0.552 | 0.766 |
| V2010_test    | Olympic  |  148866 |  591669 |      0.500 |      0.492 |            0.496 |            0.501 |          0.500 |  0.510 | 0.499 |
| V2010_test    | Arm B    |    1807 |    7527 |      0.514 |      0.502 |            0.508 |            0.517 |          0.513 |  0.528 | 0.526 |
| V2008_test    | ALL      |   93754 |  363318 |      0.507 |      0.509 |            0.513 |            0.511 |          0.507 |  0.516 | 0.490 |
| V2008_test    | NW       |   16322 |   54185 |      0.546 |      0.585 |            0.585 |            0.560 |          0.544 |  0.581 | 0.598 |
| V2008_test    | SW       |     305 |     581 |      0.665 |      0.672 |            0.654 |            0.662 |          0.683 |  0.516 | 0.792 |
| V2008_test    | Olympic  |   73842 |  295583 |      0.497 |      0.492 |            0.497 |            0.499 |          0.498 |  0.500 | 0.464 |
| V2008_test    | Arm B    |    3285 |   12969 |      0.495 |      0.490 |            0.493 |            0.502 |          0.491 |  0.513 | 0.511 |
| A2011         | ALL      |  294824 | 1149490 |      0.503 |      0.499 |            0.507 |            0.507 |          0.501 |  0.518 | 0.514 |
| A2011         | NW       |   61999 |  222561 |      0.522 |      0.558 |            0.557 |            0.537 |          0.520 |  0.558 | 0.562 |
| A2011         | SW       |     709 |    2799 |      0.501 |      0.493 |            0.505 |            0.500 |          0.490 |  0.511 | 0.482 |
| A2011         | Olympic  |  215645 |  858740 |      0.496 |      0.482 |            0.492 |            0.498 |          0.496 |  0.507 | 0.501 |
| A2011         | Arm B    |   16471 |   65390 |      0.494 |      0.487 |            0.497 |            0.502 |          0.489 |  0.508 | 0.508 |

**WA fine-tune (ft1) — AUC, own-family truth-matched vs unmatched data** (reference columns: step-1 `n_z2` and step-1 `P`)

| eval_set      | region   |   n_pos |   n_neg |   pnet_max |   pnet_2nd |   n_sites_pnet03 |   n_sites_pnet05 |   pnet_rel_max |   n_z2 |     P |
|:--------------|:---------|--------:|--------:|-----------:|-----------:|-----------------:|-----------------:|---------------:|-------:|------:|
| fit_in_sample | ALL      |   39611 |  120580 |      0.627 |      0.697 |            0.692 |            0.631 |          0.614 |  0.936 | 0.947 |
| fit_in_sample | NW       |   24903 |   18305 |      0.675 |      0.708 |            0.700 |            0.677 |          0.664 |  0.920 | 0.948 |
| fit_in_sample | SW       |    6685 |     168 |      0.634 |      0.600 |            0.587 |            0.586 |          0.642 |  0.948 | 0.989 |
| fit_in_sample | Olympic  |    5181 |  100803 |      0.681 |      0.562 |            0.574 |            0.617 |          0.669 |  0.786 | 0.806 |
| fit_in_sample | Arm B    |    2842 |    1304 |      0.598 |      0.659 |            0.652 |            0.591 |          0.588 |  0.975 | 0.992 |
| V2010_test    | ALL      |    7761 |  107540 |      0.745 |      0.692 |            0.698 |            0.705 |          0.727 |  0.639 | 0.626 |
| V2010_test    | NW       |    3662 |    8865 |      0.753 |      0.764 |            0.749 |            0.753 |          0.741 |  0.677 | 0.710 |
| V2010_test    | Olympic  |    4013 |   97343 |      0.709 |      0.567 |            0.582 |            0.635 |          0.697 |  0.531 | 0.519 |
| V2010_test    | Arm B    |      86 |    1332 |      0.674 |      0.706 |            0.696 |            0.677 |          0.654 |  0.761 | 0.785 |
| V2008_test    | ALL      |    3857 |   10086 |      0.777 |      0.775 |            0.766 |            0.763 |          0.766 |  0.669 | 0.711 |
| V2008_test    | NW       |    3708 |    9939 |      0.781 |      0.777 |            0.769 |            0.765 |          0.768 |  0.665 | 0.705 |
| V2008_test    | SW       |     149 |     147 |      0.745 |      0.761 |            0.723 |            0.741 |          0.748 |  0.558 | 0.819 |
| A2011         | ALL      |    8865 |   53134 |      0.780 |      0.790 |            0.771 |            0.772 |          0.763 |  0.694 | 0.640 |
| A2011         | NW       |    8865 |   53134 |      0.780 |      0.790 |            0.771 |            0.772 |          0.763 |  0.694 | 0.640 |

**WA fine-tune (ft1) — AUC, ≤ 3 km truth-matched vs unmatched data** (reference columns: step-1 `n_z2` and step-1 `P`)

| eval_set      | region   |   n_pos |   n_neg |   pnet_max |   pnet_2nd |   n_sites_pnet03 |   n_sites_pnet05 |   pnet_rel_max |   n_z2 |     P |
|:--------------|:---------|--------:|--------:|-----------:|-----------:|-----------------:|-----------------:|---------------:|-------:|------:|
| fit_in_sample | ALL      |   49027 |  111164 |      0.627 |      0.668 |            0.666 |            0.618 |          0.614 |  0.850 | 0.858 |
| fit_in_sample | NW       |   24945 |   18263 |      0.676 |      0.708 |            0.700 |            0.678 |          0.665 |  0.920 | 0.948 |
| fit_in_sample | SW       |    6685 |     168 |      0.634 |      0.600 |            0.587 |            0.586 |          0.642 |  0.948 | 0.989 |
| fit_in_sample | Olympic  |   14544 |   91440 |      0.630 |      0.556 |            0.566 |            0.576 |          0.617 |  0.609 | 0.596 |
| fit_in_sample | Arm B    |    2853 |    1293 |      0.597 |      0.659 |            0.652 |            0.589 |          0.587 |  0.974 | 0.991 |
| V2010_test    | ALL      |   16425 |   98876 |      0.666 |      0.613 |            0.620 |            0.624 |          0.653 |  0.566 | 0.542 |
| V2010_test    | NW       |    3686 |    8841 |      0.754 |      0.765 |            0.750 |            0.754 |          0.742 |  0.677 | 0.711 |
| V2010_test    | Olympic  |   12647 |   88709 |      0.631 |      0.556 |            0.566 |            0.578 |          0.619 |  0.520 | 0.484 |
| V2010_test    | Arm B    |      92 |    1326 |      0.653 |      0.683 |            0.669 |            0.668 |          0.634 |  0.755 | 0.787 |
| V2008_test    | ALL      |    3892 |   10051 |      0.778 |      0.776 |            0.768 |            0.764 |          0.767 |  0.669 | 0.713 |
| V2008_test    | NW       |    3739 |    9908 |      0.782 |      0.778 |            0.771 |            0.766 |          0.769 |  0.666 | 0.707 |
| V2008_test    | SW       |     153 |     143 |      0.741 |      0.770 |            0.733 |            0.749 |          0.748 |  0.546 | 0.829 |
| A2011         | ALL      |    8959 |   53040 |      0.780 |      0.791 |            0.772 |            0.773 |          0.763 |  0.695 | 0.641 |
| A2011         | NW       |    8959 |   53040 |      0.780 |      0.791 |            0.772 |            0.773 |          0.763 |  0.695 | 0.641 |

## B. Fine-tune on Washington

Label rule: (family, station node) pairs whose build-half stack (≥ 30 events) has a clear-ish arrival (imp ≥ 3 and dur < 3 s; catalogue clear = imp > 5); label weight 1.0 if clear, 0.5 otherwise; every build-half array detection of the family gives a Gaussian (σ 0.5 s) target at t + template offset; loss only within ±3 s of a labelled arrival in positive windows. Negatives: quiet-day windows and windows at the stations' scrambled-control candidate positions (build half), masked ±3 s around every known arrival. Held-out stations: UW.HDW, UW.OSD, PB.B018, PB.B017, US.NLWA, UW.GNW (each its own site); held-out days: V2008/V2010 test half and A2011.

Training windows (BUILD half of V2008 + V2010, 36 training stations, 68 days):

| kind                       |   n_horiz |   windows |
|:---------------------------|----------:|----------:|
| positive                   |         0 |     25848 |
| positive                   |         2 |     78174 |
| quiet-day negative         |         0 |     41454 |
| quiet-day negative         |         2 |     65096 |
| scrambled-control negative |         0 |     88094 |
| scrambled-control negative |         2 |    137805 |

Fine-tune history (internal validation = every 10th date of the BUILD half; step 0 = SVI weights):

|      step |   val_loss |   val_auc |   pos_med |   neg_p99 |   train_loss |
|----------:|-----------:|----------:|----------:|----------:|-------------:|
|     0.000 |      0.269 |     0.621 |     0.045 |     0.551 |      nan     |
|  1000.000 |      0.177 |     0.843 |     0.222 |     0.381 |        0.311 |
|  2000.000 |      0.171 |     0.851 |     0.211 |     0.359 |        0.306 |
|  3000.000 |      0.172 |     0.853 |     0.218 |     0.371 |        0.303 |
|  4000.000 |      0.173 |     0.848 |     0.218 |     0.391 |        0.297 |
|  5000.000 |      0.169 |     0.849 |     0.207 |     0.359 |        0.294 |
|  6000.000 |      0.173 |     0.848 |     0.229 |     0.427 |        0.291 |
|  7000.000 |      0.176 |     0.846 |     0.241 |     0.453 |        0.287 |
|  8000.000 |      0.174 |     0.846 |     0.224 |     0.439 |        0.283 |
|  9000.000 |      0.176 |     0.848 |     0.243 |     0.475 |        0.279 |
| 10000.000 |      0.176 |     0.846 |     0.245 |     0.500 |        0.276 |
| 11000.000 |      0.176 |     0.839 |     0.220 |     0.527 |        0.275 |
| 12000.000 |      0.178 |     0.839 |     0.229 |     0.520 |        0.270 |
| 13000.000 |      0.176 |     0.844 |     0.220 |     0.464 |        0.268 |
| 14000.000 |      0.182 |     0.836 |     0.242 |     0.578 |        0.265 |
| 15000.000 |      0.181 |     0.836 |     0.220 |     0.516 |        0.262 |
| 16000.000 |      0.179 |     0.841 |     0.218 |     0.508 |        0.259 |
| 17000.000 |      0.182 |     0.835 |     0.224 |     0.551 |        0.257 |
| 18000.000 |      0.187 |     0.829 |     0.236 |     0.605 |        0.254 |
| 19000.000 |      0.188 |     0.833 |     0.239 |     0.613 |        0.252 |
| 20000.000 |      0.191 |     0.823 |     0.242 |     0.621 |        0.250 |

Kept checkpoint: step 3000 (internal val AUC 0.853).

Held-out days (V2008/V2010 TEST half, A2011): per-station AUC summary, stations with ≥ 50 labelled arrivals (`held_out` = the 6 stations never used for fine-tuning):

| held_out   | tag   |   n_stations |   med_auc_vs_quiet |   med_auc_vs_episode |   med_fp03_per_h |   med_fp05_per_h |   med_pos_frac03 |
|:-----------|:------|-------------:|-------------------:|---------------------:|-----------------:|-----------------:|-----------------:|
| False      | ft1   |           30 |              0.845 |                0.690 |            7.488 |            1.014 |            0.134 |
| False      | w05   |           30 |              0.698 |                0.545 |           29.099 |            8.741 |            0.181 |
| True       | ft1   |            6 |              0.802 |                0.655 |           17.985 |            3.292 |            0.174 |
| True       | w05   |            6 |              0.706 |                0.602 |           93.265 |           29.038 |            0.278 |

Per station (≥ 50 labelled arrivals on held-out days):

| netsta   | tag   | held_out   |   n_pos |   auc_vs_quiet |   auc_vs_episode |   auc_vs_episode_V2008 |   auc_vs_episode_V2010 |   auc_vs_episode_A2011 |   pos_frac03 |   fp03_per_h |   fp05_per_h |
|:---------|:------|:-----------|--------:|---------------:|-----------------:|-----------------------:|-----------------------:|-----------------------:|-------------:|-------------:|-------------:|
| PB.B017  | ft1   | True       |     659 |          0.575 |            0.576 |                  0.587 |                nan     |                nan     |        0.082 |       22.725 |        5.489 |
| PB.B017  | w05   | True       |     659 |          0.655 |            0.586 |                  0.615 |                nan     |                nan     |        0.205 |      162.808 |       45.886 |
| PB.B018  | ft1   | True       |    6573 |          0.809 |            0.620 |                  0.763 |                  0.426 |                nan     |        0.247 |       25.997 |        7.784 |
| PB.B018  | w05   | True       |    6573 |          0.817 |            0.619 |                  0.715 |                  0.483 |                nan     |        0.389 |      132.020 |       40.710 |
| US.NLWA  | ft1   | True       |    3574 |          0.809 |            0.690 |                  0.697 |                  0.674 |                nan     |        0.382 |       99.048 |       30.489 |
| US.NLWA  | w05   | True       |    3574 |          0.879 |            0.681 |                  0.713 |                  0.668 |                nan     |        0.352 |       54.511 |       17.366 |
| UW.GNW   | ft1   | True       |    6885 |          0.795 |            0.743 |                  0.793 |                  0.754 |                nan     |        0.389 |        8.312 |        0.878 |
| UW.GNW   | w05   | True       |    6885 |          0.740 |            0.718 |                  0.745 |                  0.722 |                nan     |        0.516 |      385.228 |       58.014 |
| UW.HDW   | ft1   | True       |   27931 |          0.522 |            0.563 |                  0.593 |                  0.623 |                  0.567 |        0.062 |       13.245 |        0.926 |
| UW.HDW   | w05   | True       |   27931 |          0.460 |            0.464 |                  0.430 |                  0.470 |                  0.503 |        0.007 |       22.251 |       12.156 |
| UW.OSD   | ft1   | True       |    3031 |          0.823 |            0.696 |                  0.658 |                  0.726 |                nan     |        0.102 |        9.175 |        1.094 |
| UW.OSD   | w05   | True       |    3031 |          0.671 |            0.570 |                  0.582 |                  0.560 |                nan     |        0.012 |       14.162 |        7.746 |
| CN.PGC   | ft1   | False      |    6147 |          0.881 |            0.692 |                  0.677 |                  0.683 |                  0.713 |        0.224 |       11.423 |        1.678 |
| CN.PGC   | w05   | False      |    6147 |          0.803 |            0.581 |                  0.557 |                  0.589 |                  0.629 |        0.223 |       25.934 |        6.302 |
| CN.VGZ   | ft1   | False      |    4865 |          0.905 |            0.734 |                  0.732 |                  0.738 |                  0.729 |        0.227 |        5.991 |        0.939 |
| CN.VGZ   | w05   | False      |    4865 |          0.821 |            0.643 |                  0.639 |                  0.672 |                  0.633 |        0.240 |       22.959 |        4.681 |
| PB.B001  | ft1   | False      |   32386 |          0.910 |            0.774 |                  0.784 |                  0.772 |                  0.773 |        0.510 |       29.498 |        2.718 |
| PB.B001  | w05   | False      |   32386 |          0.843 |            0.642 |                  0.623 |                  0.670 |                  0.639 |        0.282 |       65.775 |       23.885 |
| PB.B003  | ft1   | False      |    1081 |          0.802 |            0.608 |                  0.664 |                  0.561 |                  0.601 |        0.109 |        4.915 |        0.756 |
| PB.B003  | w05   | False      |    1081 |          0.805 |            0.487 |                  0.428 |                  0.492 |                  0.543 |        0.112 |       36.421 |       13.592 |
| PB.B005  | ft1   | False      |    4087 |          0.959 |            0.823 |                  0.806 |                  0.823 |                  0.847 |        0.267 |        4.975 |        1.051 |
| PB.B005  | w05   | False      |    4087 |          0.923 |            0.734 |                  0.713 |                  0.734 |                  0.761 |        0.401 |       43.483 |       12.119 |
| PB.B006  | ft1   | False      |    3370 |          0.927 |            0.792 |                  0.752 |                  0.798 |                  0.814 |        0.431 |       11.669 |        2.093 |
| PB.B006  | w05   | False      |    3370 |          0.886 |            0.711 |                  0.685 |                  0.716 |                  0.729 |        0.442 |       62.920 |       17.529 |
| PB.B007  | ft1   | False      |    6956 |          0.932 |            0.788 |                  0.780 |                  0.805 |                  0.799 |        0.324 |        8.579 |        1.167 |
| PB.B007  | w05   | False      |    6956 |          0.877 |            0.708 |                  0.678 |                  0.727 |                  0.738 |        0.409 |       66.975 |       18.774 |
| PB.B009  | ft1   | False      |    2282 |          0.904 |            0.741 |                  0.703 |                  0.711 |                  0.797 |        0.265 |        6.623 |        1.486 |
| PB.B009  | w05   | False      |    2282 |          0.884 |            0.705 |                  0.669 |                  0.710 |                  0.769 |        0.409 |       70.052 |       23.846 |
| PB.B010  | ft1   | False      |    4726 |          0.819 |            0.608 |                  0.603 |                  0.602 |                  0.636 |        0.183 |        3.455 |        0.828 |
| PB.B010  | w05   | False      |    4726 |          0.813 |            0.544 |                  0.499 |                  0.537 |                  0.613 |        0.263 |       25.934 |        9.108 |
| PB.B011  | ft1   | False      |    2082 |          0.918 |            0.764 |                  0.738 |                  0.765 |                  0.787 |        0.290 |        4.127 |        0.978 |
| PB.B011  | w05   | False      |    2082 |          0.909 |            0.696 |                  0.665 |                  0.679 |                  0.756 |        0.402 |       36.637 |       12.262 |
| PB.B013  | ft1   | False      |   58490 |          0.909 |            0.832 |                  0.849 |                  0.869 |                  0.804 |        0.629 |       67.994 |       20.757 |
| PB.B013  | w05   | False      |   58490 |          0.784 |            0.674 |                  0.738 |                  0.606 |                  0.693 |        0.418 |       41.103 |        7.086 |
| PB.B019  | ft1   | False      |     108 |          0.932 |            0.841 |                  0.477 |                nan     |                nan     |        0.046 |        2.086 |        0.482 |
| PB.B019  | w05   | False      |     108 |          0.527 |            0.552 |                  0.603 |                nan     |                nan     |        0.324 |      340.910 |       36.005 |
| PB.B941  | ft1   | False      |     179 |          0.573 |            0.509 |                  0.470 |                  0.541 |                nan     |        0.011 |        8.918 |        1.917 |
| PB.B941  | w05   | False      |     179 |          0.592 |            0.538 |                  0.454 |                  0.599 |                nan     |        0.140 |      116.327 |       34.264 |
| PB.B943  | ft1   | False      |   59301 |          0.909 |            0.838 |                  0.859 |                  0.873 |                  0.805 |        0.614 |       55.882 |       15.926 |
| PB.B943  | w05   | False      |   59301 |          0.793 |            0.677 |                  0.747 |                  0.607 |                  0.694 |        0.429 |       32.263 |        6.545 |
| UW.BHW   | ft1   | False      |     100 |          0.380 |            0.414 |                  0.459 |                nan     |                nan     |        0.000 |       19.179 |        1.676 |
| UW.BHW   | w05   | False      |     100 |          0.379 |            0.452 |                  0.483 |                nan     |                nan     |        0.000 |       25.270 |       13.277 |
| UW.BLN   | ft1   | False      |   23636 |          0.848 |            0.720 |                  0.688 |                  0.745 |                  0.721 |        0.171 |       14.594 |        2.028 |
| UW.BLN   | w05   | False      |   23636 |          0.411 |            0.512 |                  0.520 |                  0.507 |                  0.504 |        0.013 |       70.298 |       17.998 |
| UW.CPW   | ft1   | False      |     906 |          0.721 |            0.591 |                  0.550 |                  0.511 |                nan     |        0.007 |        1.187 |        0.098 |
| UW.CPW   | w05   | False      |     906 |          0.475 |            0.493 |                  0.501 |                  0.526 |                nan     |        0.003 |        3.958 |        1.898 |
| UW.DOSE  | ft1   | False      |     649 |          0.842 |            0.741 |                nan     |                  0.929 |                  0.732 |        0.159 |        1.074 |        0.307 |
| UW.DOSE  | w05   | False      |     649 |          0.829 |            0.712 |                nan     |                  0.865 |                  0.712 |        0.313 |       21.274 |        1.244 |
| UW.GMW   | ft1   | False      |    1942 |          0.814 |            0.687 |                  0.617 |                  0.737 |                nan     |        0.066 |        0.489 |        0.013 |
| UW.GMW   | w05   | False      |    1942 |          0.406 |            0.474 |                  0.494 |                  0.455 |                nan     |        0.002 |        4.571 |        2.119 |
| UW.HTW   | ft1   | False      |     928 |          0.566 |            0.520 |                  0.507 |                  0.531 |                nan     |        0.009 |        5.302 |        0.385 |
| UW.HTW   | w05   | False      |     928 |          0.502 |            0.503 |                  0.507 |                  0.501 |                nan     |        0.004 |       11.954 |        6.841 |
| UW.LRIV  | ft1   | False      |    1980 |          0.949 |            0.816 |                nan     |                  0.831 |                  0.807 |        0.378 |        4.222 |        0.650 |
| UW.LRIV  | w05   | False      |    1980 |          0.862 |            0.745 |                nan     |                  0.757 |                  0.740 |        0.368 |       52.982 |        5.889 |
| UW.OTR   | ft1   | False      |    1848 |          0.707 |            0.588 |                  0.633 |                  0.571 |                nan     |        0.034 |        6.979 |        0.288 |
| UW.OTR   | w05   | False      |    1848 |          0.615 |            0.526 |                  0.588 |                  0.536 |                nan     |        0.017 |       15.997 |        8.375 |
| UW.PCFR  | ft1   | False      |    2563 |          0.528 |            0.515 |                  0.511 |                  0.520 |                nan     |        0.010 |       18.243 |        1.933 |
| UW.PCFR  | w05   | False      |    2563 |          0.493 |            0.490 |                  0.493 |                  0.488 |                nan     |        0.002 |        5.471 |        3.081 |
| UW.RER   | ft1   | False      |     247 |          0.555 |            0.552 |                  0.527 |                  0.535 |                nan     |        0.000 |        8.145 |        0.809 |
| UW.RER   | w05   | False      |     247 |          0.518 |            0.478 |                  0.540 |                  0.466 |                nan     |        0.000 |        3.489 |        1.911 |
| UW.RRHS  | ft1   | False      |     534 |          0.561 |            0.540 |                  0.530 |                  0.542 |                nan     |        0.004 |        9.430 |        1.503 |
| UW.RRHS  | w05   | False      |     534 |          0.562 |            0.545 |                  0.551 |                  0.557 |                nan     |        0.002 |        2.751 |        1.490 |
| UW.SMW   | ft1   | False      |    1995 |          0.806 |            0.664 |                  0.727 |                  0.667 |                  0.621 |        0.023 |        2.524 |        0.215 |
| UW.SMW   | w05   | False      |    1995 |          0.505 |            0.440 |                  0.394 |                  0.421 |                  0.464 |        0.008 |        5.197 |        2.563 |
| UW.SQM   | ft1   | False      |    4772 |          0.894 |            0.842 |                  0.815 |                  0.845 |                nan     |        0.405 |       35.526 |        8.286 |
| UW.SQM   | w05   | False      |    4772 |          0.782 |            0.811 |                  0.745 |                  0.777 |                nan     |        0.228 |       65.168 |       16.142 |
| UW.STW   | ft1   | False      |    2527 |          0.863 |            0.639 |                  0.625 |                  0.568 |                  0.446 |        0.000 |        0.815 |        0.000 |
| UW.STW   | w05   | False      |    2527 |          0.379 |            0.465 |                  0.479 |                  0.523 |                  0.435 |        0.004 |       82.910 |       20.750 |
| UW.SVOH  | ft1   | False      |     534 |          0.413 |            0.475 |                  0.479 |                  0.475 |                nan     |        0.002 |       13.791 |        1.411 |
| UW.SVOH  | w05   | False      |     534 |          0.530 |            0.534 |                  0.595 |                  0.500 |                nan     |        0.011 |       24.924 |       13.765 |
| UW.VVHS  | ft1   | False      |     433 |          0.483 |            0.492 |                  0.455 |                  0.528 |                nan     |        0.000 |        7.998 |        0.817 |
| UW.VVHS  | w05   | False      |     433 |          0.511 |            0.512 |                  0.507 |                  0.516 |                nan     |        0.005 |        9.795 |        5.417 |

False positives on QUIET days (test half): median over stations of declustered CF peaks per hour:

| held_out   | tag   |   fp03 |   fp05 |
|:-----------|:------|-------:|-------:|
| False      | ft1   |  6.979 |  1.051 |
| False      | w05   | 36.421 | 11.022 |
| True       | ft1   | 17.985 |  3.292 |
| True       | w05   | 93.265 | 29.038 |

Verdict rule (fixed before the run): fine-tune 'works' if, over held-out stations, the median AUC gain is > 0.02 vs quiet AND > 0 vs episode random times. Median ΔAUC vs quiet 0.024, vs episode 0.017 → **works = True**.

## C. Step-1 refit with PNET features

Same lfdr design (HistGradientBoosting data-vs-control + isotonic, 4 scrambled controls, π0 from quiet days, fit on the BUILD half); only the feature list differs (+ `pnet_max`, `pnet_2nd`, `n_sites_pnet03`, `n_sites_pnet05`). Held-out sets only. slope = calibration slope of the truth-matched fraction on P (target 0.8–1.2); P05 = fraction of candidates with P > 0.5; Spearman = expected count ΣP per family-day vs the true array count.

| run                  | eval_set   | region   |   slope_own |   slope_3km |   matched_over_sumP_own |   data_P05 |   ctl_P05 |   ctl_P05_never_seen |   data_sumP |   ctl_sumP_per_real |   spearman_own |   spearman_3km |
|:---------------------|:-----------|:---------|------------:|------------:|------------------------:|-----------:|----------:|---------------------:|------------:|--------------------:|---------------:|---------------:|
| baseline (step 1)    | V2010_test | ALL      |       0.650 |       0.578 |                   1.287 |      0.021 |     0.004 |                0.005 |    8795.923 |            5952.251 |          0.491 |          0.551 |
| baseline (step 1)    | V2010_test | NW       |       0.556 |       0.561 |                   1.482 |      0.106 |     0.011 |                0.015 |    2723.345 |             561.451 |          0.514 |          0.469 |
| baseline (step 1)    | V2010_test | SW       |     nan     |     nan     |                 nan     |      0.386 |     0.024 |                0.042 |     158.714 |              12.950 |        nan     |        nan     |
| baseline (step 1)    | V2010_test | Olympic  |       0.208 |       0.193 |                   1.153 |      0.007 |     0.003 |                0.004 |    5827.517 |            5328.555 |          0.373 |          0.312 |
| baseline (step 1)    | V2010_test | Arm B    |       0.596 |       0.602 |                   1.375 |      0.032 |     0.010 |                0.016 |      86.347 |              49.295 |          0.180 |          0.086 |
| baseline (step 1)    | V2008_test | ALL      |       0.556 |       0.566 |                   1.546 |      0.029 |     0.007 |                0.007 |    5110.635 |            3437.870 |          0.452 |          0.404 |
| baseline (step 1)    | V2008_test | NW       |       0.550 |       0.560 |                   1.589 |      0.128 |     0.008 |                0.013 |    2411.336 |             386.823 |          0.519 |          0.488 |
| baseline (step 1)    | V2008_test | SW       |       0.679 |       0.700 |                   0.870 |      0.518 |     0.071 |                0.100 |     162.531 |              20.076 |          0.278 |          0.236 |
| baseline (step 1)    | V2008_test | Olympic  |     nan     |     nan     |                 nan     |      0.006 |     0.007 |                0.006 |    2477.650 |            2978.939 |        nan     |        nan     |
| baseline (step 1)    | V2008_test | Arm B    |     nan     |     nan     |                 nan     |      0.001 |     0.001 |                0.002 |      59.119 |              52.031 |        nan     |        nan     |
| baseline (step 1)    | A2011      | ALL      |       0.544 |       0.552 |                   1.047 |      0.028 |     0.014 |                0.036 |   24821.464 |           20152.425 |          0.476 |          0.453 |
| baseline (step 1)    | A2011      | NW       |       0.544 |       0.552 |                   1.047 |      0.085 |     0.017 |                0.032 |    8464.658 |            3961.714 |          0.476 |          0.453 |
| baseline (step 1)    | A2011      | SW       |     nan     |     nan     |                 nan     |      0.477 |     0.488 |                0.657 |     362.880 |             369.081 |        nan     |        nan     |
| baseline (step 1)    | A2011      | Olympic  |     nan     |     nan     |                 nan     |      0.011 |     0.011 |                0.036 |   14677.801 |           14610.222 |        nan     |        nan     |
| baseline (step 1)    | A2011      | Arm B    |     nan     |     nan     |                 nan     |      0.023 |     0.016 |                0.028 |    1316.125 |            1211.407 |        nan     |        nan     |
| + PNET off-the-shelf | V2010_test | ALL      |       0.635 |       0.564 |                   1.274 |      0.019 |     0.003 |                0.005 |    8711.919 |            5703.909 |          0.498 |          0.551 |
| + PNET off-the-shelf | V2010_test | NW       |       0.548 |       0.553 |                   1.409 |      0.113 |     0.012 |                0.015 |    2838.936 |             559.868 |          0.522 |          0.472 |
| + PNET off-the-shelf | V2010_test | SW       |     nan     |     nan     |                 nan     |      0.417 |     0.030 |                0.042 |     164.144 |              11.911 |        nan     |        nan     |
| + PNET off-the-shelf | V2010_test | Olympic  |       0.182 |       0.176 |                   1.178 |      0.004 |     0.002 |                0.003 |    5621.329 |            5085.482 |          0.379 |          0.329 |
| + PNET off-the-shelf | V2010_test | Arm B    |       0.614 |       0.623 |                   1.288 |      0.033 |     0.011 |                0.017 |      87.510 |              46.648 |          0.175 |          0.082 |
| + PNET off-the-shelf | V2008_test | ALL      |       0.553 |       0.563 |                   1.492 |      0.031 |     0.007 |                0.005 |    5172.381 |            3326.741 |          0.459 |          0.412 |
| + PNET off-the-shelf | V2008_test | NW       |       0.546 |       0.555 |                   1.534 |      0.136 |     0.010 |                0.013 |    2484.876 |             373.830 |          0.526 |          0.495 |
| + PNET off-the-shelf | V2008_test | SW       |       0.737 |       0.752 |                   0.827 |      0.548 |     0.083 |                0.107 |     168.321 |              20.305 |          0.284 |          0.246 |
| + PNET off-the-shelf | V2008_test | Olympic  |     nan     |     nan     |                 nan     |      0.006 |     0.007 |                0.004 |    2460.678 |            2884.566 |        nan     |        nan     |
| + PNET off-the-shelf | V2008_test | Arm B    |     nan     |     nan     |                 nan     |      0.003 |     0.001 |                0.002 |      58.506 |              48.040 |        nan     |        nan     |
| + PNET off-the-shelf | A2011      | ALL      |       0.543 |       0.552 |                   1.060 |      0.033 |     0.020 |                0.044 |   28725.232 |           23912.579 |          0.473 |          0.449 |
| + PNET off-the-shelf | A2011      | NW       |       0.543 |       0.552 |                   1.060 |      0.082 |     0.015 |                0.024 |    8359.841 |            3945.652 |          0.473 |          0.449 |
| + PNET off-the-shelf | A2011      | SW       |     nan     |     nan     |                 nan     |      0.203 |     0.210 |                0.569 |     253.507 |             261.016 |        nan     |        nan     |
| + PNET off-the-shelf | A2011      | Olympic  |     nan     |     nan     |                 nan     |      0.020 |     0.021 |                0.049 |   18846.491 |           18555.223 |        nan     |        nan     |
| + PNET off-the-shelf | A2011      | Arm B    |     nan     |     nan     |                 nan     |      0.020 |     0.013 |                0.018 |    1265.393 |            1150.687 |        nan     |        nan     |
| + PNET fine-tuned    | V2010_test | ALL      |       0.618 |       0.547 |                   1.081 |      0.026 |     0.005 |                0.005 |   10366.146 |            6806.889 |          0.505 |          0.551 |
| + PNET fine-tuned    | V2010_test | NW       |       0.569 |       0.574 |                   1.193 |      0.135 |     0.012 |                0.014 |    3341.371 |             634.919 |          0.526 |          0.472 |
| + PNET fine-tuned    | V2010_test | SW       |     nan     |     nan     |                 nan     |      0.457 |     0.031 |                0.045 |     176.157 |              12.148 |        nan     |        nan     |
| + PNET fine-tuned    | V2010_test | Olympic  |       0.149 |       0.128 |                   1.002 |      0.009 |     0.004 |                0.004 |    6735.418 |            6102.006 |          0.366 |          0.322 |
| + PNET fine-tuned    | V2010_test | Arm B    |       0.651 |       0.661 |                   0.965 |      0.043 |     0.011 |                0.017 |     113.200 |              57.816 |          0.205 |          0.102 |
| + PNET fine-tuned    | V2008_test | ALL      |       0.582 |       0.591 |                   1.331 |      0.034 |     0.007 |                0.006 |    5889.410 |            3768.516 |          0.467 |          0.418 |
| + PNET fine-tuned    | V2008_test | NW       |       0.578 |       0.586 |                   1.357 |      0.153 |     0.010 |                0.012 |    2818.539 |             426.907 |          0.539 |          0.505 |
| + PNET fine-tuned    | V2008_test | SW       |       0.738 |       0.752 |                   0.880 |      0.518 |     0.086 |                0.120 |     166.676 |              20.485 |          0.282 |          0.245 |
| + PNET fine-tuned    | V2008_test | Olympic  |     nan     |     nan     |                 nan     |      0.007 |     0.007 |                0.005 |    2822.154 |            3254.064 |        nan     |        nan     |
| + PNET fine-tuned    | V2008_test | Arm B    |     nan     |     nan     |                 nan     |      0.005 |     0.001 |                0.003 |      82.042 |              67.060 |        nan     |        nan     |
| + PNET fine-tuned    | A2011      | ALL      |       0.547 |       0.556 |                   0.894 |      0.026 |     0.008 |                0.021 |   26879.437 |           21251.527 |          0.485 |          0.460 |
| + PNET fine-tuned    | A2011      | NW       |       0.547 |       0.556 |                   0.894 |      0.091 |     0.012 |                0.016 |    9918.017 |            4616.186 |          0.485 |          0.460 |
| + PNET fine-tuned    | A2011      | SW       |     nan     |     nan     |                 nan     |      0.142 |     0.165 |                0.324 |     247.377 |             253.443 |        nan     |        nan     |
| + PNET fine-tuned    | A2011      | Olympic  |     nan     |     nan     |                 nan     |      0.007 |     0.007 |                0.023 |   15432.163 |           15244.039 |        nan     |        nan     |
| + PNET fine-tuned    | A2011      | Arm B    |     nan     |     nan     |                 nan     |      0.016 |     0.006 |                0.011 |    1281.880 |            1137.859 |        nan     |        nan     |

PNET features inside the refits (permutation log-loss gain, rank among all features):

| run                  | feature        |   perm_logloss_gain |   auc_data_vs_ctl |   rank |   of |
|:---------------------|:---------------|--------------------:|------------------:|-------:|-----:|
| + PNET off-the-shelf | n_sites_pnet03 |              0.0007 |            0.5225 |     14 |   25 |
| + PNET off-the-shelf | pnet_2nd       |              0.0007 |            0.5198 |     15 |   25 |
| + PNET off-the-shelf | pnet_max       |              0.0003 |            0.5157 |     21 |   25 |
| + PNET off-the-shelf | n_sites_pnet05 |              0.0000 |            0.5181 |     24 |   25 |
| + PNET fine-tuned    | pnet_max       |              0.0013 |            0.5089 |     10 |   25 |
| + PNET fine-tuned    | pnet_2nd       |              0.0012 |            0.5168 |     11 |   25 |
| + PNET fine-tuned    | n_sites_pnet03 |              0.0009 |            0.5198 |     15 |   25 |
| + PNET fine-tuned    | n_sites_pnet05 |              0.0001 |            0.5140 |     24 |   25 |

## D. PNET-derived candidate times → T2 withheld-station template recovery

Withheld XU/XG station, template rebuilt from candidate times only; cc (2–8 Hz) vs the station's true template. `data` = beam over the 3 nearest permanent sites, `data_n1` = nearest site only; `control*` = the same on the site-rolled scrambled control; `ceiling` = stack at the true test-half detections.

| source                                  | region   | series     |   n_pairs |   med_cc_al |   frac_cc_gt07 |   frac_cc_gt05 |   own_best |   best_nbr3 |   med_rank |   med_pur_own |   med_n_stacked |
|:----------------------------------------|:---------|:-----------|----------:|------------:|---------------:|---------------:|-----------:|------------:|-----------:|--------------:|----------------:|
| PNET beam, WA fine-tune (ft1)           | ALL      | control    |      1003 |       0.140 |          0.000 |          0.000 |      0.007 |       0.028 |     94.000 |         0.000 |         773.000 |
| PNET beam, WA fine-tune (ft1)           | ALL      | control_n1 |      1003 |       0.138 |          0.000 |          0.000 |      0.006 |       0.025 |     88.000 |         0.000 |         774.000 |
| PNET beam, WA fine-tune (ft1)           | ALL      | data       |      1003 |       0.149 |          0.001 |          0.002 |      0.007 |       0.033 |     81.000 |         0.004 |         771.000 |
| PNET beam, WA fine-tune (ft1)           | ALL      | data_n1    |      1003 |       0.146 |          0.001 |          0.013 |      0.008 |       0.032 |     79.000 |         0.001 |         776.000 |
| PNET beam, off-the-shelf (SVI pnet_w05) | ALL      | control    |      1003 |       0.138 |          0.000 |          0.000 |      0.010 |       0.027 |     90.000 |         0.000 |         774.000 |
| PNET beam, off-the-shelf (SVI pnet_w05) | ALL      | control_n1 |      1003 |       0.137 |          0.000 |          0.000 |      0.007 |       0.022 |     93.000 |         0.000 |         776.000 |
| PNET beam, off-the-shelf (SVI pnet_w05) | ALL      | data       |      1003 |       0.142 |          0.000 |          0.000 |      0.005 |       0.019 |     85.000 |         0.000 |         778.000 |
| PNET beam, off-the-shelf (SVI pnet_w05) | ALL      | data_n1    |      1003 |       0.139 |          0.000 |          0.000 |      0.006 |       0.015 |     86.000 |         0.000 |         776.000 |
| T2 permanent MF (round 0)               | ALL      | ceiling    |      1003 |       0.659 |          0.383 |          0.819 |      0.685 |       0.719 |      1.000 |       nan     |         119.000 |
| T2 permanent MF (round 0)               | ALL      | control    |      1003 |       0.139 |          0.000 |          0.000 |      0.007 |       0.026 |     85.000 |         0.000 |         771.000 |
| T2 permanent MF (round 0)               | ALL      | data       |      1003 |       0.207 |          0.095 |          0.167 |      0.196 |       0.220 |     21.000 |         0.032 |         773.000 |
| PNET beam, WA fine-tune (ft1)           | Arm B    | control    |       123 |       0.138 |          0.000 |          0.000 |      0.008 |       0.049 |     97.000 |         0.000 |         816.000 |
| PNET beam, WA fine-tune (ft1)           | Arm B    | control_n1 |       123 |       0.132 |          0.000 |          0.000 |      0.000 |       0.049 |     89.000 |         0.000 |         816.000 |
| PNET beam, WA fine-tune (ft1)           | Arm B    | data       |       123 |       0.138 |          0.000 |          0.000 |      0.008 |       0.033 |     83.000 |         0.002 |         816.000 |
| PNET beam, WA fine-tune (ft1)           | Arm B    | data_n1    |       123 |       0.139 |          0.000 |          0.000 |      0.016 |       0.057 |    102.000 |         0.001 |         816.000 |
| PNET beam, off-the-shelf (SVI pnet_w05) | Arm B    | control    |       123 |       0.131 |          0.000 |          0.000 |      0.016 |       0.073 |     97.000 |         0.000 |         816.000 |
| PNET beam, off-the-shelf (SVI pnet_w05) | Arm B    | control_n1 |       123 |       0.136 |          0.000 |          0.000 |      0.000 |       0.041 |     90.000 |         0.000 |         816.000 |
| PNET beam, off-the-shelf (SVI pnet_w05) | Arm B    | data       |       123 |       0.131 |          0.000 |          0.000 |      0.000 |       0.041 |     84.000 |         0.000 |         816.000 |
| PNET beam, off-the-shelf (SVI pnet_w05) | Arm B    | data_n1    |       123 |       0.132 |          0.000 |          0.000 |      0.024 |       0.073 |     82.000 |         0.000 |         816.000 |
| T2 permanent MF (round 0)               | Arm B    | ceiling    |       123 |       0.568 |          0.163 |          0.675 |      0.748 |       0.813 |      1.000 |       nan     |          78.000 |
| T2 permanent MF (round 0)               | Arm B    | control    |       123 |       0.138 |          0.000 |          0.000 |      0.008 |       0.073 |     73.000 |         0.000 |         816.000 |
| T2 permanent MF (round 0)               | Arm B    | data       |       123 |       0.147 |          0.000 |          0.000 |      0.008 |       0.033 |     64.000 |         0.004 |         816.000 |
| PNET beam, WA fine-tune (ft1)           | NW       | control    |       627 |       0.138 |          0.000 |          0.000 |      0.005 |       0.014 |     96.000 |         0.000 |         816.000 |
| PNET beam, WA fine-tune (ft1)           | NW       | control_n1 |       627 |       0.138 |          0.000 |          0.000 |      0.005 |       0.013 |     89.000 |         0.000 |         815.000 |
| PNET beam, WA fine-tune (ft1)           | NW       | data       |       627 |       0.150 |          0.000 |          0.000 |      0.003 |       0.008 |     83.000 |         0.007 |         816.000 |
| PNET beam, WA fine-tune (ft1)           | NW       | data_n1    |       627 |       0.147 |          0.000 |          0.000 |      0.006 |       0.018 |     79.000 |         0.003 |         816.000 |
| PNET beam, off-the-shelf (SVI pnet_w05) | NW       | control    |       627 |       0.139 |          0.000 |          0.000 |      0.010 |       0.021 |     84.000 |         0.000 |         815.000 |
| PNET beam, off-the-shelf (SVI pnet_w05) | NW       | control_n1 |       627 |       0.138 |          0.000 |          0.000 |      0.008 |       0.014 |     92.000 |         0.000 |         816.000 |
| PNET beam, off-the-shelf (SVI pnet_w05) | NW       | data       |       627 |       0.145 |          0.000 |          0.000 |      0.005 |       0.010 |     87.000 |         0.001 |         816.000 |
| PNET beam, off-the-shelf (SVI pnet_w05) | NW       | data_n1    |       627 |       0.143 |          0.000 |          0.000 |      0.005 |       0.008 |     81.000 |         0.000 |         816.000 |
| T2 permanent MF (round 0)               | NW       | ceiling    |       627 |       0.675 |          0.424 |          0.834 |      0.612 |       0.628 |      1.000 |       nan     |         129.000 |
| T2 permanent MF (round 0)               | NW       | control    |       627 |       0.138 |          0.000 |          0.000 |      0.008 |       0.014 |     85.000 |         0.000 |         816.000 |
| T2 permanent MF (round 0)               | NW       | data       |       627 |       0.268 |          0.097 |          0.211 |      0.254 |       0.271 |      6.000 |         0.048 |         816.000 |
| PNET beam, WA fine-tune (ft1)           | Olympic  | control    |        45 |       0.117 |          0.000 |          0.000 |      0.000 |       0.156 |     98.000 |         0.000 |         816.000 |
| PNET beam, WA fine-tune (ft1)           | Olympic  | control_n1 |        45 |       0.119 |          0.000 |          0.000 |      0.044 |       0.178 |    101.000 |         0.000 |         816.000 |
| PNET beam, WA fine-tune (ft1)           | Olympic  | data       |        45 |       0.234 |          0.022 |          0.044 |      0.044 |       0.467 |     17.000 |         0.005 |         816.000 |
| PNET beam, WA fine-tune (ft1)           | Olympic  | data_n1    |        45 |       0.405 |          0.022 |          0.289 |      0.000 |       0.200 |     20.000 |         0.015 |         816.000 |
| PNET beam, off-the-shelf (SVI pnet_w05) | Olympic  | control    |        45 |       0.129 |          0.000 |          0.000 |      0.000 |       0.044 |    126.000 |         0.000 |         816.000 |
| PNET beam, off-the-shelf (SVI pnet_w05) | Olympic  | control_n1 |        45 |       0.115 |          0.000 |          0.000 |      0.000 |       0.111 |    137.000 |         0.001 |         816.000 |
| PNET beam, off-the-shelf (SVI pnet_w05) | Olympic  | data       |        45 |       0.110 |          0.000 |          0.000 |      0.000 |       0.089 |    150.000 |         0.000 |         816.000 |
| PNET beam, off-the-shelf (SVI pnet_w05) | Olympic  | data_n1    |        45 |       0.112 |          0.000 |          0.000 |      0.000 |       0.022 |    154.000 |         0.000 |         816.000 |
| T2 permanent MF (round 0)               | Olympic  | ceiling    |        45 |       0.873 |          0.733 |          0.911 |      0.889 |       0.956 |      1.000 |       nan     |         100.000 |
| T2 permanent MF (round 0)               | Olympic  | control    |        45 |       0.130 |          0.000 |          0.000 |      0.022 |       0.156 |     93.000 |         0.000 |         816.000 |
| T2 permanent MF (round 0)               | Olympic  | data       |        45 |       0.856 |          0.756 |          0.800 |      0.689 |       0.800 |      1.000 |         0.055 |         816.000 |
| PNET beam, WA fine-tune (ft1)           | SW       | control    |       208 |       0.145 |          0.000 |          0.000 |      0.014 |       0.029 |     87.500 |         0.000 |         686.000 |
| PNET beam, WA fine-tune (ft1)           | SW       | control_n1 |       208 |       0.145 |          0.000 |          0.000 |      0.005 |       0.014 |     82.000 |         0.000 |         684.000 |
| PNET beam, WA fine-tune (ft1)           | SW       | data       |       208 |       0.147 |          0.000 |          0.000 |      0.010 |       0.014 |     85.000 |         0.001 |         685.500 |
| PNET beam, WA fine-tune (ft1)           | SW       | data_n1    |       208 |       0.144 |          0.000 |          0.000 |      0.010 |       0.024 |     80.000 |         0.001 |         685.000 |
| PNET beam, off-the-shelf (SVI pnet_w05) | SW       | control    |       208 |       0.140 |          0.000 |          0.000 |      0.010 |       0.014 |     95.000 |         0.000 |         686.000 |
| PNET beam, off-the-shelf (SVI pnet_w05) | SW       | control_n1 |       208 |       0.138 |          0.000 |          0.000 |      0.010 |       0.014 |     95.000 |         0.000 |         685.000 |
| PNET beam, off-the-shelf (SVI pnet_w05) | SW       | data       |       208 |       0.144 |          0.000 |          0.000 |      0.010 |       0.019 |     73.000 |         0.000 |         686.000 |
| PNET beam, off-the-shelf (SVI pnet_w05) | SW       | data_n1    |       208 |       0.141 |          0.000 |          0.000 |      0.000 |       0.000 |     83.500 |         0.000 |         685.000 |
| T2 permanent MF (round 0)               | SW       | ceiling    |       208 |       0.644 |          0.312 |          0.837 |      0.822 |       0.885 |      1.000 |       nan     |         109.000 |
| T2 permanent MF (round 0)               | SW       | control    |       208 |       0.141 |          0.000 |          0.000 |      0.000 |       0.005 |     89.500 |         0.000 |         686.000 |
| T2 permanent MF (round 0)               | SW       | data       |       208 |       0.148 |          0.000 |          0.000 |      0.029 |       0.053 |     90.500 |         0.011 |         685.000 |

