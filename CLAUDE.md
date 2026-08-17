# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this repository is

ANSYMB2024 is a mobile-brain/body-imaging (MoBI) research project, not a software product. Participants (`sub-5` … `sub-18`) perform seated knee flexion/extension angle-tracking while wearing a knee exoskeleton driven by a pneumatic artificial muscle (PAM) at three assistance pressures (1, 3, 6 bar → `P1`/`P3`/`P6`). Each trial is also rated by the participant on a subjective difficulty scale of 1–10 (`Score`). Recorded modalities: 64-ch EEG (LiveAmp, 500 Hz), 6-ch EMG (Delsys Trigno, 2000 Hz), and an "experiment" stream (knee encoder angle, reference angle, force sensor, beep/trigger channel) — all captured together via LSL into `.xdf`.

The code is **MATLAB scripts, run interactively**. There is no build, no linter config, no test suite, and no package manager. Scripts are written as `%%`-section notebooks meant to be stepped through cell by cell in the MATLAB desktop; large blocks are commented out on purpose because that stage was "run once already" (this is a deliberate convention throughout — read the comments before uncommenting anything). Several stages block on GUI apps (`.mlapp`) or `questdlg`/`inputdlg` prompts.

`Code/Python/main_EMG_to_LSL.py` is standalone acquisition tooling (streams Delsys Trigno EMG to LSL); it is not part of the analysis pipeline.

## Running code

- Everything runs from MATLAB on Windows. External toolboxes live **outside** the repo and are added by absolute path inside each script:
  - EEGLAB (`D:\Morteza\Toolboxes\EEGLAB\eeglab2026.0.0`, older scripts point at `eeglab2025.1.0` / `eeglab2025.0.0` / `eeglab2024.2.1`)
  - FieldTrip (`D:\Morteza\Toolboxes\Fieldtrip\fieldtrip-20231127`)
  - xdf-Matlab (`D:\Morteza\LSL\xdf-Matlab-master`)
  - BeMoBIL pipeline — git submodule at `Code/Matlab/data_processing/BeMoBIL_Pipeline` (currently uninitialized; also listed in `.gitignore`)
  - EEGLAB plugins in use: AMICA, ICLabel, dipfit, Zapline-Plus, SIFT (connectivity), and MATLAB's Statistics & ML, Signal Processing, and Wavelet toolboxes (`fitlme`, `cwt`, `kmeans`, Classification Learner exports).
- The MATLAB MCP server is available (`mcp__matlab__*`). Prefer `check_matlab_code` for static checks. Be careful with `run_matlab_file` / `evaluate_matlab_code`: most top-level scripts start with `clc; clear`, load multi-GB `.mat`/`.set` files, take minutes to hours, open blocking GUIs, or overwrite saved results — ask before running one.
- There is no automated test. "Verifying" a change here means running the relevant `%%` section on one subject and eyeballing the figure or the saved structure.

## Data flow

Raw and derived data live in `data/` (gitignored, not in the repo). Folder numbers *are* the pipeline stages; scripts address them by these literal names:

| Stage | Produced by | Contents |
| --- | --- | --- |
| `0_source_data/sub-N/ses-S00X/eeg/*.xdf` | recording | raw LSL streams (3 sessions per subject), plus `Subjects.xlsx` metadata |
| `1_BIDS_data` | `Main.m` → `bemobil_xdf2bids` | BIDS export |
| `2_raw-EEGLAB` | `Main.m` → `bemobil_bids2set` | `sub-N_merged_EEG.set` |
| `3_EEG-preprocessing`, `4_spatial-filters/4-1_AMICA` | `bemobil_process_all_EEG_preprocessing`, `bemobil_process_all_AMICA` | cleaning, ICA, dipfit |
| `5_single-subject-EEG-analysis` | BeMoBIL wrappers | `sub-N_cleaned_with_ICA.set`; `timewarp_test/Epoched_data` holds the newer EEGLAB-native epoched/time-warped sets + `.icatimef` |
| `6_0_Trials_Info_and_Events` | `Main_add_events` + `find_flexion_extension_events.mlapp` | `sub-N_Trials_encoder_events.mat` |
| `6_Trials_Info_and_Epoched_data` | `Main_evevnt_and_epoch_selection.m` | `Trials_Info.mat`, `Epochs_{FlextoFlex,Flexion,Extension,Trial}_based.mat` |
| `7_STUDY` | `Group_Level_PostProcessing` | EEGLAB `.study` files, `multiple_clustering/<ROI>/`, ERSP results |
| `8_Classification/ROIs_features` | `Classification_analysis` | `ROIs_*.mat` feature tables |
| `9_EXP_Analysis`, `10_Time_Frequency_Analysis` | `EXP_analysis`, `Time_Frequency_Analysis` | torque/force structures, TF content |

Some derived artifacts are stored **next to the code** rather than in `data/`, e.g. `EMG_processing/structured_EMG_data/sub-N/`, `Classification_analysis/Classifier_training/sub-N/` (Classification Learner exports, accuracy encoded in the filename), `Events/sub-N/events_*.txt`, and the `_NHB/*/*.mat` result tables.

## Pipeline entry points (in order)

1. **`data_processing/Main.m`** — one subject at a time: `runs_concatenated` (loads/concatenates all `.xdf` sessions) → xdf2bids → bids2set → drop ACC channels 65–67 → load `chanlocs.ced` → `Main_add_events` → reject non-experimental periods → BeMoBIL preprocessing + AMICA. Parameters come from **`BeMoBIL_Configuration.m`**, the single source of truth for folder names, filenames, channel cleaning, AMICA, dipfit, and ICLabel settings.
2. **`Main_add_events.m`** — derives trial events from the experiment stream and imports them into `EEG.event`. Event types: `SB_Start_Beep`, `PC_Pressure_Change`, `SM_Start_Move`, `FB_Finish_Beep`, `SP_Score_Press`, and the movement events `FlxS`/`FlxE`/`ExtS`/`ExtE`. Each event's `desc` field ends in `_<trialnumber>`, which is how downstream code re-associates events with trials.
3. **`Main_evevnt_and_epoch_selection.m`** — `Main_event_selection` builds `Trials_Info`, then `Main_epoch_selection` cuts epochs. Only the FlexToFlex variant is enabled by default.
4. **`Group_Level_PostProcessing/main_group_level_postprocessing.m`** — builds the STUDY, preclusters (dipole weight 3, scalp 1), then runs BeMoBIL **repeated k-means clustering** (200–1000 iterations) once per anatomical ROI. Each ROI is defined by an MNI coordinate + cluster count + a 6-element `quality_measure_weights` vector; results land in `7_STUDY/.../multiple_clustering/<ROI name>/`. The commented-out block at the bottom sweeps N clusters to pick the best.
5. **Analyses** branching off the epoched data:
   - `EMG_processing/main_EMG_processing.m` — builds `structured_EMG_data` (RMS/iEMG features, time-warping, outlier flags).
   - `Time_Frequency_Analysis/main_TF_analysis.m` — Morlet/`cwt` TF per IC per ROI.
   - `Classification_analysis/` — per-subject decoding of pressure (`P1P3P6`, `P1P6`) or score cluster (`S1S2S3`). Two pipelines exist, both programmatic: **Pipeline B (`regardless_of_GroupLevel_clustering/`, June 2025) is the most recent complete run** — PSD band integrals over all brain ICs, no group-level clustering, no SMOTE, via `runClassificationWithoutSMOTE.m`. Pipeline A (`main_classification.m`, March 2025) uses band-RMS features over group-clustered ICs with SMOTE, via `runClassificationWithSMOTE.m`. Both runners use ANOVA/eta-squared feature selection, 10 repeats of 5-fold CV, and `fitcauto` (SVM / naive Bayes / neural net, 30 evaluations, inner 5-fold on the training fold) — so there are 50 model fits per subject per feature set and no single "chosen" architecture behind a reported accuracy. `Based_on_Scores_Clustering/Score_Clustering.m` k-means-clusters each subject's 1–10 scores into S1/S2/S3, with hardcoded per-subject initial centers. All 71 saved result logs come from these two runners. The MATLAB Classification Learner app branch (`main_classifier_training*.m`, `Classifier_training/`) predates both and produced none of the current results — it only assembles tables into the base workspace for manual use.
   - `EXP_analysis/` — force/torque and kinematics, incl. OpenSim.
   - `Cortico_Muscular_Coherence/`, `Detailed_per_subject_processing/` — per-subject CMC and diagnostic plots.
6. **`Group_Level_PostProcessing/Final_paper_plot_generation/`** — figure/stat generation for the manuscript. **`_NHB/` is the current paper's analysis set** and the most actively developed area: `behavioural_results/` (builds `behavior_table.mat`), `manual_TF_outlier_removal/` (→ `EEG_features.mat`), `full_LMM_Analysis_Err_EMG_EEG/FULL_LMM_PIPELINE.m` (merges behavior + EEG, `fitlme` mixed models), `Corticomuscular_Coherence/`, `connectivity_analysis/` (SIFT), `sPCA_denoising_of_TF_data/`, `BCI_classifier/`, `Supplementary_section/`.

## Core data structures

- **`Trials_Info{1,i}`** (one cell per trial): `.General` = `Description` (`'Experiment'` / contains `'Reject'`), `Session`, `Pressure` (1/3/6), `Score` (1–10, `0` = not rated), `Case` (3 or 4 — see gotchas). `.Events` mirrors the same index fields (`Trial_start_indx`, `Pressure_Change_indx`, `Movement_start/end_indx`, `flexion_*`, `extension_*`, `flextoflex_*`) **three times over**, once per stream: `EEG_stream.Raw` (indices into the concatenated raw stream), `EEG_stream.Preprocessed` (latencies in the cleaned `.set`, after non-experimental segments were cut), `EMG_stream`, and `EXP_stream`. Mixing these index spaces up is the most common source of misalignment bugs.
- **`Epochs_FlextoFlex_based{1,i}`**: `.General`, `.EEG_stream.{Raw,Preprocessed}.{Times,Channels,Sources}`, `.EMG_stream.{Times,Sensors_Raw,Sensors_Preprocessed}`, `.EXP_stream.{Times,Forces,Encoder_angle,Ref_angle}`. Each field is a cell array with one entry per movement cycle within the trial. EMG preprocessing is fixed at bandpass 20–450 Hz → rectify → 4 Hz lowpass (Butterworth, `filtfilt`).
- **Muscle order** is positional everywhere: `{Vastus_med_R, Rectus_femoris_R, Gastrocnemius_R, Biceps_femoris_R, Trapezius_R, Trapezius_L}`, from Delsys sensor ids `[2 3 4 5 6 7]`.
- **`ROIs_*.mat`**: struct keyed by brain-region name; each value is an N×3+ cell array of `{subject_id, IC_index, per-epoch feature struct}`.

## Conventions and gotchas

- **Absolute paths are hardcoded in every script.** Two roots coexist: newer code uses `D:\Morteza\MyProjects\ANSYMB2024`, older code still says `C:\Morteza\...` (≈25 files). When editing a script, follow the root already in that file; when a script fails at `addpath`, that mismatch is usually why.
- **`subject_id > 9` branches everywhere.** The experimental protocol changed after sub-9: the beep encoding differs (single vs. double `diff` values on `All_Experiment(6,:)`), pressure change moves from 2 s *before* to 2 s *after* the start beep, and epoch structures are saved differently. Any new code touching raw streams needs both branches.
- **`Case` 3 vs 4** describes whether a trial ends on a high or low encoder peak; it decides whether the last `Flexion_Start` is dropped when building flex-to-flex cycles. Always handle both.
- **Subject exclusions are per-analysis and inconsistent.** `5:18` is the nominal list; EMG analyses drop sub-10, CMC drops sub-8 and sub-10, and force-sensor/PAM analyses only have `[11 12 15 16 17 18]`. Copy the list from the analysis you are extending rather than assuming.
- **`Main.m` is stale**: it calls `Main_add_events` with 8 arguments, but the function was reduced to 5 (`EEG, output, subject_id, processing_path, study_path`) on 01/04/2026.
- **A folder rename is mid-flight**: `Classifiction_analysis` → `Classification_analysis` is staged in git, but 6 scripts still build paths from the old misspelling and will fail at runtime.
- Data files (`*.mat`, `*.set`, `*.fdt`, `*.xdf`) and the `data/` tree are gitignored — only code, small text/figure artifacts, and trained-model `.mat` files under the code tree are tracked. Backup-by-copy is the norm (`old_codes/`, `Old files/`, `temp_codes/`, `Events - Backup at 01.04.2026/`, `*_old.mat`); leave those alone unless asked.
- Figure colors for the three pressure conditions are fixed across the paper: P1 `[1,115,178]/255`, P3 `[222,143,5]/255`, P6 `[148,73,92]/255`.
- Paper-level EEG features use alpha/mu `[8 14]` Hz and beta `[14 30]` Hz over a movement cycle time-normalized to 0–100 %, where 0–50 % is flexion and 50–100 % is extension.


## Known methodological issue (priority)

Feature selection (ANOVA, eta-squared > 0.05) runs ONCE over the full dataset
before the CV split, in both runClassificationWithSMOTE.m and
runClassificationWithoutSMOTE.m (selection fixed at :60-63; fold loop starts
at :88). Every "Selected" accuracy in every saved log is optimistically biased.
The "Full" column is unaffected. Fixing this is the current work.

Pipeline A additionally z-scores before the split
(createInputTableforClassification.m:74).

Full audit: ANSYMB2024_Classification_Report.docx at the repo root.


- NEVER modify or delete anything under *_results/ or
  subjects_classification_results/. Those are the historical outputs the
  report documents; they are the only record of the March and June 2025 runs.