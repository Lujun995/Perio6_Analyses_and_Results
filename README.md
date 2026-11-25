# Perio6 Analyses and Results

Analysis code and results for *Machine Learning Enhances Accuracy of Partial Mouth Diagnostics for Periodontitis: A Multi-Cohort Study*. The workflow prepares periodontal examination data, fits the Perio6 generative classifier, searches for partial-mouth tooth combinations, and compares classifier performance across study cohorts.

## Browse the results

Open [`Year-2-Wrap-Ups.html`](Year-2-Wrap-Ups.html) in a browser after cloning or downloading the repository. Keep [`Year-2-Wrap-Ups_files/`](Year-2-Wrap-Ups_files/) next to the HTML for its figures. The [`res/`](res/) directory contains selected figures, fitted NHANES models, and processed data for the analysis.

| File or directory | Purpose |
| --- | --- |
| [`Year-2 Wrap-Ups.Rmd`](Year-2%20Wrap-Ups.Rmd) | Report source for cohort preparation, model fitting, tooth selection, validation, and figures. |
| [`src/Model Fit.R`](src/Model%20Fit.R) and [`src/Model Predict.R`](src/Model%20Predict.R) | Study workflows for fitting the ordinal models and running prediction. |
| [`src/validation.R`](src/validation.R) and [`src/calculate_metrics2.R`](src/calculate_metrics2.R) | Cohort validation and performance metrics. |
| [`src/`](src/) | Classifier code, comparison models, evolutionary search scripts (`EA_*.R`), and plotting helpers. |
| [`tooth numbering/`](tooth%20numbering/) | Tooth-numbering SVGs. `Universal_Numbering_System.svg` is the template used by `paint_tooth.R`. |

For the reusable classifier and prediction entry points, see the companion [Perio6 source repository](https://github.com/Lujun995/Perio6).

## Reproduce the analysis

The rendered report can be viewed as shipped. A full rerun requires additional inputs and adjustments to the study scripts. Use a working copy because the report writes generated files into `res/`.

1. Set up R and Pandoc. The report uses `rmarkdown` and `knitr`, `tidyverse` and `haven` for data preparation, and `ordinalNet`, `ranger`, `ROCR`, `glmnet`, `e1071`, and `pROC` for modelling and evaluation. Tables and figures use packages including `kableExtra`, `tableone`, `Hmisc`, `ggtext`, `ggrepel`, `cowplot`, `ggpubr`, `igraph`, and the Bioconductor package `ggtree`. The tooth plots also use `xml2`, `rsvg`, and `magick`. Individual search scripts declare their own dependencies.
2. Place a copy of [`Year-2 Wrap-Ups.Rmd`](Year-2%20Wrap-Ups.Rmd) in a `research/` subdirectory of the working copy. Its `../src/`, `../res/`, and cohort paths assume that layout. Preserve the filename case shown in the report when setting up inputs on a case-sensitive filesystem.
3. Supply the additional inputs listed below. Check the hard-coded `workdir` and worker counts in the evolutionary search and validation scripts before running them. Before regenerating validation outputs, check cohort selection and how `src/validation.R` assigns output filenames so each file can be traced to the correct cohort, tooth combination, and repetition. Verify fold provenance before interpreting results as out-of-fold estimates.
4. With the inputs in place and the scripts configured, render from the root of the working copy.

   ```r
   rmarkdown::render("research/Year-2 Wrap-Ups.Rmd")
   ```

### Additional inputs for a full rerun

| Expected location | Input |
| --- | --- |
| `NHANES/` | NHANES 2009–2014 dental examination (`OHXPER`), demographics (`DEMO`), oral health questionnaire (`OHQ`), and smoking questionnaire (`SMQ`) files for cycles F, G, and H. |
| `ARIC/AS1996.01 Dental Study V4/` | Authorized ARIC dental and derived demographic SAS datasets. |
| `Temple/` | Authorized Temple summaries and prediction artifacts used by the report. |
| `res/` | Evolutionary search results and validation prediction RDS files referenced by the report. Only selected results and fitted models are included in the checkout. |
| `ORIGINS/` | Authorized ORIGINS SAS data, needed only when `enable_origins` is enabled in the report. The default is `FALSE`. |

## License and data access

The code is released under the [MIT License](LICENSE). Obtain cohort inputs through their respective data providers and follow the applicable access and redistribution terms for participant-level data and derived outputs.
