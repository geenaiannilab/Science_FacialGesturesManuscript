This repository contains all code necessary to regenerate the figures found in:  Ianni G.R., Vázquez Y., Rouse A.G., Schieber M.H., Prut Y., Freiwald W.A. Facial gestures are enacted via a cortical hierarchy of dynamic and stable codes. 
 
 
 
Figure 1 
Facial gestures are distinguishable during a naturalistic social paradigm
1B : Facial posture prior to movement onset 
1C: TSNE plot of marker positions 
1D-F: Face PCs by weights, by temporal resolution, by trajectory 
Data: markerTrajectories.mat, markerTrajectoriesTSNE.mat
Plotting: plotBehavioralTrajectories.m
Analysis: generateBehavioralTrajectories.m
 
Figure 2
Single-Cell Activity and Selectivity in Cortical Face-Motor Regions. 
2B: Sample population activity of all simultaneously recorded cells and facial markers 
Data: Fig2B.mat 
Plotting: plotNeuralActivityMatrix.m
Analysis: plotRasterAllCells2021.m
 
2C: Peri-event time histograms of four neurons aligned to onset of three gestures
Data: Fig2C.mat
Plotting: plotPSTHs.m
 
2D: Bee plot of gesture preference indices (GPI) distributions for each cortical area 
Data: Fig2D_GPIvalues.mat 
Plotting: plotGPIByRegion.m
Analysis: gesturePreferenceIndex.m
 
2E: Fractions of cells in each region with significant activity modulations by facial gesture type, time, and their interaction 
Data: Fig2E.mat
Plotting: plotBhvTimeANOVAbyRegion.m
Analysis: facePrefANOVA.m, chiSquaredGOF.m
 
Figure 3 
Population encoding of Facial Gestures
3A: Single-event population activity vectors cluster by gesture-type in neural state space 
Data: Fig3A.mat
Plotting: plotSingleTrialsNeuralStateSpace.m
Analysis: calculateSingleTrialPCAallCells2021.m 
 
3B: Decoding categorical gesture type from neural population activity
Data: Fig3B.mat
Plotting: plotSVMClassifierAcrossDays.m
Analysis: analyzeSVMClassifierAcrossDays.m, SVMClassifierAcrossDays.m, SVMClassifierPermutationTestingAcrossDays.m
 
3C, 3D: Trial-averaged, time resolved gesture-specific neural trajectories 
Data: Fig3C.mat
Plotting: plotTrialAvgPCAallCells.m
Analysis: calculateTrialAvgPCAallCellsWithDistances.m
 
Figure 4
Kinematic decoding and unique neural correlations during gestures
4A, 4B : Kinematic decoding of facial gesture components by region
Data Fig4AResults.mat, Fig4ANullResults.mat
Plotting script: plotPSIDResultsVNull.m
Analysis: prepData2023.m, makePseduoPopPSID23.m, analyzePSIDWrapper.m (runs runCorePSIDAnalysis.m underneath)
Dependencies: /PSID/
 
4C: Strength of pairwise cross-gesture neural correlations (R2, threat v. lipsmack, threat v. chew, and lipsmack v. chew, one dot per recording) for each region 
Data:  Fig4C.mat
Plotting: plotCorrStructureAcrossExpressions.m
Analysis: calcR2BetweenPairwiseCorrelations.m, neuralCorrelationStructureAcrossExpressions_wRankOrder.m
 
4D: Population response dissimilarity matrices, all regions 
Data: Fig4D/(region).mat
Plotting: plotCorrStructureAcrossExpressions.m
Analysis: neuralCorrelationStructureAcrossExpressions_wRankOrder.m
 
4E: Spike-triggered movement averages 
Data, plotted: STMA/(pair).fig 
Analysis: STMA2021wrapper.m, nullSTMAwrapper.m, evaluateSTMAwrapper.m 
Dependencies: /STMA/
 
Figure 5 
Stable and dynamic coding of facial gestures across cortex 
Left plots 
Data: Fig5_left.mat
Plotting Scripts: plotCrossTemporalDecodingAcrossDays.m, 
Analysis: analyzeCrossTemporalSVMClassifierAcrossDays.m, 
Dependencies: crossTemporalSVMClassifierAcrossDays.m
 
Right plots
Data: neuralTrajectory_combined.mat
Plotting: neuralTrajectoryCharacteristicsAcrossDays.m
Analysis: calculateNeuralTrajectoryAcrossDays.m
 
SUPPLEMENTAL FIGURES 
 
Figure S1
Handscored gesture onsets compared to automatic detection by continuous facial marker tracking
Data: markerTrajectories.mat
Analysis, Plotting: compareManualScoring.m 
 
Figure S2
Mutual information distributions by region 
Data: FigS2.mat
Plotting: plotMutualInformationStatic.m
 
Figure S3A
Between-region categorical decoding accuracy over time
Data: Fig3B.mat
Plotting: plotSVMClassifierAcrossDays.m
 
Figure S3B
Categorical decoding curves by region 
Data: Fig3B.mat
Plotting: plotSVMClassifierAcrossDays.m
 
Figure S4A, S4B
Euclidean distances between gesture-specific neural trajectories, per day
Data: FigS4.mat
Plotting: summarizeEuclideanDistancesOverDays.m
 
Figure S4C
Euclidean distances between gesture-specific neural trajectories, across pseudopopulations 
Data: /Figure5/matfiles/neuralTrajectory_combined.mat
Plotting: /Figure5/neuralTrajectoryCharacteristicsAcrossDays.m
 
Figure S5
Kinematic Decoding Performance by Region, with statistical comparisons overlaid 
Data: /Figure4/matfiles/Fig4AResults.mat, /Figure4/matfilesFig4ANullResults.mat
Plotting script: plotPSIDResultsVNull.m
 
Figure S6A, S6B
Rank order similarity of gesture-specific neural correlations by region
Data: /Supplemental/matfiles/combined_rankOrderStats.csv
Plotting: summarizeRankOrderStats.m
 
Figure S7 
Diagonal Index, Temporal Generalization Windows of Cross-temporal Generalization Matrices by Region 
Data: /Figure5/matfiles/Fig5_left.mat
Plotting: metricsCrossTemporalDecoding.m
 
Figure S8
Examples of per day, per region behavioral clustering of individual trials in neural state space, separated by region 
Data: /Figure3/matfiles/Fig3A_(region).mat
Plotting: plotSingleTrialsNeuralStateSpace.m
 
Figure S9 
Examples of per day, per region gesture trajectories in 3D neural state space, separated by region 
Data: Fig3C.mat
Plotting: plotTrialAvgPCAallCells.m
 
Figure S10
Data: /matfiles/FigS10.mat
Analysis Scripts:  /Supplemental/calcMutualInformationTimeResolved.m 
Plotting Scripts: /Supplemental/plotMutualInformationTimeResolved.m 
