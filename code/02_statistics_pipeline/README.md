## DOD gait statistic analysis pipeline
- [DOD gait statistic analysis pipeline](#dod-gait-statistic-analysis-pipeline)
  - [01\_MakeTable\_Demographics.r](#01_maketable_demographicsr)
  - [02\_MakeTable\_GaitMeasures.r](#02_maketable_gaitmeasuresr)
  - [03\_MakeFigure\_GaitMeasures\_6panel.r](#03_makefigure_gaitmeasures_6panelr)
  - [04\_MakeTable\_EEGMeasures.r](#04_maketable_eegmeasuresr)
  - [05\_MakeTable\_EEGMeasures\_AffectedSide.r](#05_maketable_eegmeasures_affectedsider)
  - [06\_MakeFigure\_EEGMeasures\_8panel.r](#06_makefigure_eegmeasures_8panelr)
  - [07\_MakeTableAndFigure\_EEGOnGaitDiff\_MoreAffected\_perCue.r](#07_maketableandfigure_eegongaitdiff_moreaffected_percuer)
  - [Helper Functions (`func/`)](#helper-functions-func)

### 01_MakeTable_Demographics.r
- Performs statistical testing on and generate summary tables for demographic data.

### 02_MakeTable_GaitMeasures.r
- Performs linear mixed-effects modeling and planned contrasts on gait kinematic measures
- (dv ~ group * cue + age + (1 | subject)) using cluster-robust standard errors.

### 03_MakeFigure_GaitMeasures_6panel.r
- Generates kinematic plots for the manuscript, including boxplots and effect sizes for gait measures.

### 04_MakeTable_EEGMeasures.r
- Performs linear mixed-effects modeling and planned contrasts on EEG measures
- (dv ~ subgroup * cue + age + (1 | subject)) using cluster-robust standard errors.

### 05_MakeTable_EEGMeasures_AffectedSide.r
- Performs linear mixed-effects modeling and planned contrasts on EEG measures specifically for more/less affected sides
- (dv ~ cue x group + age + pas + updrs3 + (1 | subject)) using cluster-robust standard errors.

### 06_MakeFigure_EEGMeasures_8panel.r
- Plots beta suppression boxplots and effect size for RSM, LSM, more, and less affected sides

### 07_MakeTableAndFigure_EEGOnGaitDiff_MoreAffected_perCue.r
- Performs linear model of changes in contralateal beta on gait parameters per cue in more and affected sides
- (gait change ~ group * beta change + age + pas + updrs3)


### Helper Functions (`func/`)
- `calculate_tstat.r`: Calculates t-statistics from estimates and standard errors.
- `club_contrast.r`: Computes robust standard errors, confidence intervals, and p-values for contrasts using `clubSandwich` and `emmeans`.