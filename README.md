# MSCT
Master's Thesis: Predicting promising treatment response in personalized cancer therapy

The Codes are used for the analysis performed in chapter 5 of the Master's thesis.

To start the scripts for the comparison study (Comparison_study_prediction_performance.R) and the application to Oslo data (Application_oslo_data.R) on abel one an use the following command (here as example of the drug 17-AAG):

R CMD BATCH --no-save --no-restore '--args drugname="AAG"' Comparison_study_prediction_performance.R AAG_script.Rout
