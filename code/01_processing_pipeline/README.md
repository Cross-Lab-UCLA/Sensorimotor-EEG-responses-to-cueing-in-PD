# DoD_Gait_Pipeline
- [DoD\_Gait\_Pipeline](#dod_gait_pipeline)
		- [Step 1. A1\_Merge\_Sets\_withStanding\_run1only](#step-1-a1_merge_sets_withstanding_run1only)
		- [Step 2. A2\_Run\_AMICA\_EEG\_repeated](#step-2-a2_run_amica_eeg_repeated)
		- [Step 3. A3\_Bootstrap\_EEG](#step-3-a3_bootstrap_eeg)
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
- some scripts requires functions from  \code\func

### Step 1. A1_Merge_Sets_withStanding_run1only
	- prompt and load subject data from subject-->raw
	- currently set to skipping rest trials
	- clean each condition using function `cleanTrial()`
    	- lowpass cufoff at 100Hz with transition band of 30Hz and passband edge of 85Hz
    	- highpass cutoff at 1.5Hx with transition band of 2Hz and passband edge at 2.5Hz
    	- down sample from 1024 to 256
    	- zapline plus cleaning
    	- remove channels that is 3 std above the median rms of all channels
    	- assign more and less affected sides labels to gait events
  	- merge conditions, removing channels excluded in any single condition.
  	- ICC using notch filter
  	- remove channels using clean_channels
  	- add back ref channel and average re-reference
	- save the merged .set in `data\01_combined_withStanding_run1only_moreAffected`

### Step 2. A2_Run_AMICA_EEG_repeated
	- create a tmp copy of the data
	- remove bad segments using jointprob on 1s windows and clean_windows with default size
	- interpolate, add ref chan, re-ref, and check rank
	- AMICA on tmp dataset using corrected rank
	- add ICA weights to the orginal data wuith chans added back in

### Step 3. A3_Bootstrap_EEG
	- warp digitized channels to standard MNI head and performing dipole fitting
	- compute correlations between IC activations and topographies across AMICA runs.
	- match ICs and identify poorly correlated ICs and mark them as NaN.
	- perform bootstrap resampling to compute robust dipole positions and rv

### Step 4. A4_DipRemove_EEG.m
	- warp digitized channels to standard MNI head (if A3 is skipped)
	- remove "non-brain" dipoles
    	- upward PSD slope between 4-40 Hz
    	- IC with poor correlation in repeated AMICA
    	- PowPowCat_ICrej flagged from Jacobsen et al, 2024
            - https://github.com/jacobsen-noelle/ExoAdapt-DualEEG-Processing
    	- outside of the brain using depth calc from ft_sourcedepth
    	- RV > 15

### Step 5. A5_Epoch_EEG
	- ask user if to perform epoch on RHS
	- epoch data from to [-1 2] seconds
	- remove epochs and gait cycle data that does not contain a complete gait cycle [ RHS LTO LHS RTO RHS]

### Step 6. A6_TFDecomp_EEG
	- jointprob to remove outlier epoches
	- uses newtimef() to extract ERSP
  	- normalized ERSP

### Step 7. A7_Target_Clusters.m
	- grab dipoles from a pre-set MNI corrdinate and make a cluster of dipoles with a pre-set radius from that corrdinate
	- check for outliers within that cluster
	- check for overlapping dipoles
	- check each dipoles for outlier gait cycles and remove them
	- normalized the ERSP data from the dipoles

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

