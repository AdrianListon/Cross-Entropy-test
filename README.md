# FlowCytoScript
A simplified complete workflow in R for analysis of high parameter flow cytometry data including the Crossentropy method.

This simplified version of the flowcytoscript (Crossentropy test) is intended to be usable by people with little to no experience in programming. All inputs are via plain language prompts in an RStudio markdown notebook. All outputs are organized in folders as before, but additionally an HTML document summarizing the results is created with each analysis. Check out "flowcytoscript.nb.html" for an example of the output.

![tsne_plot_all_conditions__cluster](https://github.com/AdrianListon/Cross-Entropy-test/assets/50043996/2bc3aeb5-0f71-4cb6-9763-f6bf64d7cb28)

**Features**  
-Clustering with histograms of expression, barcharts and tables of frequencies  
-Automated cluster identification  
-tSNE and UMAP  
-PCA on MFIs to show sample-level differences  
-Crossentropy testing on tSNE and UMAP  
-Statistical testing on markers and clusters  
-Heatmaps, dendrograms  
-  

**Improvements**  
-Speed. Optimizations throughout should render this approximately 10x faster, although this will vary depending on multithreading.  
-Both FCS and CSV files are accepted as input types.  
-FCS data are automatically transformed as best befits the cytometer used. This avoids potentially serious issues with scaling of the data by inexperienced users.  
-Clustering can be performed using Phenograph or FlowSOM (via EmbedSOM).  
-Clusters are automatically identified and named via matching to a cell type database.  
-tSNE performed in line with OptSNE modifications to learning rate.  
-  
![umap_plot_all_conditions__cluster](https://github.com/AdrianListon/Cross-Entropy-test/assets/50043996/1faccbe6-b2db-422e-ac8e-29d30fbfaf2c)

![histogram_by_sample](https://github.com/AdrianListon/Cross-Entropy-test/assets/50043996/b7fe83ad-bf06-4f30-a695-afa8dde1a804)

![pca_sample_mfi_loadings](https://github.com/AdrianListon/Cross-Entropy-test/assets/50043996/b5dfe9a6-ad3b-45a0-bd0f-2d7d6173ddd6)

![density_cluster](https://github.com/AdrianListon/Cross-Entropy-test/assets/50043996/f872cfba-f65d-4e14-babe-b63074fcc566)

**Using the script:**
Install R, Rstudio and Rtools. For Mac, you’ll need command line tools and OpenMP. The flowcytoscript_setup.r script (in 00_source_files) can be used to facilitate set-up for new users of R.

In your favorite flow cytometry data analysis program (FlowJo, FCS Express), gate on the cells you wish to analyze and export those cells in new fcs or csv files. While exporting, adding group or variable tags to the file names will help you sort the files with the script. 

To export your data in CSV format, preserving the transformations from FlowJo, see the instructions in "Exporting data in csv format.PNG". For more, see https://docs.flowjo.com/flowjo/graphs-and-gating/gw-transform-overview/

Create a folder for your analysis (preferably not in Dropbox or OneDrive). In this folder, put these items:
A copy of flowcytoscript.Rmd
A copy of 00_source_files
Your files, inside a folder called “Data”
Double click on the flowcytoscript file to open it in Rstudio.
Run each code chunk in order by clicking on the green arrow in the upper right corner of the chunk.

Read through the presentation "Simplified flowcytoscript--instructions for use.pptx" for more detail.

To read the publication on the Crossentropy test:
https://www.cell.com/cell-reports-methods/pdfExtended/S2667-2375(22)00295-8
