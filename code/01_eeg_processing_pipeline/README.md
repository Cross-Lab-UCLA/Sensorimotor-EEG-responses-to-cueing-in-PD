# DoD_Gait_Pipeline
- [DoD\_Gait\_Pipeline](#dod_gait_pipeline)
		- [Step 1. A1\_Merge\_Sets\_withStanding\_run1only](#step-1-a1_merge_sets_withstanding_run1only)
		- [Step 2. A2\_Run\_AMICA\_EEG\_repeated](#step-2-a2_run_amica_eeg_repeated)
		- [Step 3. A3\_Cluster\_EEG](#step-3-a3_cluster_eeg)
		- [Step 4. A4\_DipRemove\_EEG.m](#step-4-a4_dipremove_eegm)
		- [Step 5. A5\_Epoch\_EEG](#step-5-a5_epoch_eeg)
		- [Step 6. A6\_TFDecomp\_EEG](#step-6-a6_tfdecomp_eeg)
		- [Step 7. A7\_Target\_Clusters.m](#step-7-a7_target_clustersm)
		- [Step 8. A8\_MakeFigure\_allersp\_percent.m](#step-8-a8_makefigure_allersp_percentm)
		- [Step 9. A9\_Extract2CSV\_mean](#step-9-a9_extract2csv_mean)
		- [Step 10. A10\_Extract2CSV\_affectedSide\_mean](#step-10-a10_extract2csv_affectedside_mean)
		- [Step 11. A11\_Extract\_PSD](#step-11-a11_extract_psd)
		- [Step 12. A12\_PSD\_Stats\_and\_Plot](#step-12-a12_psd_stats_and_plot)


Notes:
- scripts are written for the PD gait cueing study
- run the scripts in the following order
- some scripts requires functions from  \code\01_eeg_processing_pipeline\func

### Step 1. A1_Merge_Sets_withStanding_run1only
	- prompt and load subject data from subject-->raw
	- currently set to skipping rest trials
	- clean each condition using function `cleanTrial()`
    	- lowpass cutoff at 100Hz with transition band of 30Hz and passband edge of 85Hz
    	- highpass cutoff at 1.5Hz with transition band of 2Hz and passband edge at 2.5Hz
    	- down sample from 1024 to 256
    	- zapline plus cleaning
    	- remove channels whose rms is more than 3 std above or below the mean rms of all channels
  	- assign more and less affected sides labels to gait events (PD only)
  	- merge conditions, removing channels excluded in any single condition.
  	- ICC using notch filter
  	- remove channels using clean_channels
  	- no interpolation or re-referencing here; both are done in A2
	- save the merged .set in `data\01_combined_withStanding_run1only_moreAffected`

### Step 2. A2_Run_AMICA_EEG_repeated
	- create a tmp copy of the data
	- remove bad segments using jointprob on 1s windows and clean_windows with default size
	- interpolate, add ref chan, re-ref, and check rank
	- AMICA on tmp dataset using corrected rank
	- add ICA weights to the orginal data wuith chans added back in

### Step 3. A3_Cluster_EEG
	- warp digitized channels to standard MNI head and perform dipole fitting on each IC across repeated AMICA runs
	- pool all ICs from all runs, build an all-pairs |corr| similarity on scalp maps and activations, and agglomerative-cluster into rank clusters
	- drop clusters recovered in fewer than `minRuns` runs
	- keep each cluster's centrotype IC as carrier; per-run dipole positions and rv are stored in `cltr` for downstream selection
	- save to `data\03_clustered`

### Step 4. A4_DipRemove_EEG.m
	- loads the `03_clustered` folder; dipole fits first if not already done
	- remove "non-brain" dipoles
    	- upward PSD slope between 4-40 Hz
    	- PowPowCat_ICrej flagged from Jacobsen et al, 2024
            - https://github.com/jacobsen-noelle/ExoAdapt-DualEEG-Processing
    	- outside of the brain using depth calc from ft_sourcedepth
    	- RV > 15%
	- consistency across repeated AMICA runs is handled in A3 (`minRuns`), not here

### Step 5. A5_Epoch_EEG
	- epoch on right heel strike (RHS); the affected-side epoching option is disabled (`PerformAffectedSide = false`)
	- epoch data from to [-1 2] seconds
	- remove epochs and gait cycle data that does not contain a complete gait cycle [ RHS LTO LHS RTO RHS]

### Step 6. A6_TFDecomp_EEG
	- jointprob to remove outlier epoches
	- uses newtimef() to extract ERSP
  	- normalized ERSP

### Step 7. A7_Target_Clusters.m
	- select ICs whose cluster dipole falls inside an ellipsoid ROI around the pre-set left/right sensorimotor MNI coordinates
	- remove outlier ICs by scalp-map (topo) correlation
	- when more than one IC remains in an ROI, keep the one with the strongest contralateral-minus-ipsilateral swing beta suppression

### Step 8. A8_MakeFigure_allersp_percent.m
	- create group ERSP plot
	- also performs cluster analysis on the ersp plots to grab significant mask
	- mask can be used for extracting region of interest for stats

### Step 9. A9_Extract2CSV_mean
	- extracts LSM and RSM data, and saves them in csv for statistical modeling

### Step 10. A10_Extract2CSV_affectedSide_mean
	- extracts more- and less-affected side data, and saves them in csv for statistical modeling

### Step 11. A11_Extract_PSD
	- extracts walk PSD, standing PSD, and fooofed 1/f removed PSD from standing and walking.
	- extracts PSD using spectopo(), with 50% overlap

### Step 12. A12_PSD_Stats_and_Plot
	- perform stats and plot fooofed PSD walking - fooofed PSD standing

