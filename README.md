# Average Treatment Effect for Misclassified Outcomes in Observational Studies

## 2025 University of Pennsylvania SUIP Project, Mentors: Kristin Linn (Associate Professor of Biostatistics), Dane Isenberg (PhD Student)

## Abstract
Measurement error in Electronic Health Records (EHR) undermines the data quality needed to support valid causal conclusions, increasing the risk of bias. When dealing with errors that may bias results within an observational setting, special considerations must be made when attempting to establish a causal relationship. When a non-random validation subset exists (of gold standard outcomes), we can define an average treatment effect (ATE) estimator for binary outcomes that leverages both error prone outcomes (silver standard) defined on everyone and the gold standard outcomes. We propose a novel estimator that incorporates both a propensity and silver standard classification model that aims to estimate the ATE with desirable statistical properties, such as low bias and nominal coverage. The estimator is tested using Monte Carlo simulation studies to validate our estimator within three scenarios (Independent and identically distributed individuals, Clustered data but individual treatment assignment, Clustered data with cluster level treatment). Our results show that our estimator adequately reaches our desired statistical properties regardless of the scenario, but is sensitive to classification model misspecification. This estimator extends previous work on ATE estimators from cluster randomized trials to observational settings.

## File Description

01_Ranomized_Cluster_Simulation.R: Simulating Dane's original paper's estimator and finding if similar results were achieved

02_Observational_IID_Simulation.R: Extending Dane's estimator for the IID scenario within an observational setting

03_Single_Simulation.R: Code used mainly for debugging purposes

04_Observational_Clustered_Simulation.R: Extending Dane's estimator for the Clustered but individual treatment assignment scenario within an observational setting

05_Observational_CA_Simulation.R: Extending Dane's estimator for the clustered and cluster level treatment assignement within an observational setting

06_Viz.R: Code used for generating boxplots

