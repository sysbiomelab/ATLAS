This files describe the folder's content                                
 *"Disease_PangenomeWAS.xml" Preliminary final figure, same layout as slide 4
 *"Enriched_depleated_functional_clusters.xml" Preliminary final figure, replacement for figure d) in slide 5. 
 
 Folder "code" has the scripts   for estimating effect size in disease cohorts.
        "effect_size.py" Python script, takes as input the metadata table, the msp abundance table csv format and a file with the cohorts talble; Outputs "Effect_size.Country_ctrl.tsv", "Effect_size.Matched_ctrl.tsv" and "Volcano_plot.Country.tsv"
        "cohort_colors.ob.txt" csv file with list of cohorts and and their disease category.
        "plot_manhatan.R" Rscript takes and input "Effect_size.Country_ctrl.tsv" and plots in pdf files manhatan-plots for all species and for species in vectors defined in line 118 and 119
        "plot_volcano_MSPdisease.R" Rscript takes as input "Volcano_plot.Country.tsv", "Effect_size.Country_ctrl.tsv", "funcMat.20190808.RData", and "funcModules.20191227.RData" files and ouputs Volcano plots for MSP and functional_clusters
	* Find the names for the functional_clusters commonly eriched/depleated across diseases.
 Folder "output" contains the output text files from "effect_size.py" script    
        "Effect_size.Country_ctrl.tsv" list of effect sizes for enriched and depleated by cohort using healthy samples from the same country                as background
        "Effect_size.Matched_ctrl.tsv" list of effect sizes for enriched and depletead by cohort using healthy samples from the same study as               background
        "Volcano_plot.Country.tsv" list of MSP and the number of cohorts where it was found enriched and the number cohorts where it was found              depleated
 Folders "Figures_PDF" contain all output figures in pdf (some in jpeg format) format(various versions and subsets of final version figures)
 	Filenames describe the type of plot and datatype.
 Folder "Tables" contain the updates for the Supplementary tables S3 and S4 in tsv formant, and tables of enriched/depleated species (func_clusters) by disease 
 	*"S3.csv" New file replacing supplementary table S3
	*"S4.csv" New file replacing suppplementary table S4

Some notes:
The disease cohorts definitions were built from information in  the metadata table.
The contol for effect size estimation were the healthy individials from the same country's cohort.
