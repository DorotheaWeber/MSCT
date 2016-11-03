# MSCT
Master's Thesis: Predicting promising treatment response in personalized cancer therapy

The Codes are used for the analysis performed in chapter 5 of the Master's thesis.

To start the scripts for the comparison study and the application to Oslo data on abel one an use the following command (here as example of the drug 17-AAG):

R CMD BATCH --no-save --no-restore '--args drugname="AAG"' Abel_script.R AAG_script.Rout
