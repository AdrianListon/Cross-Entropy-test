# Cross-Entropy-test
Unified script for applying the t-SNE diff method on flow cytometry or single cell RNA-seq data

Code is available and free for academic users. Commercial users should contact Adrian Liston to discuss licensing options.

For flow cytometry analysis, we recommend using the CSV-based method. This allows users to set the data scales (transformations) 
in FlowJo or any other standard flow cytometry data analysis software. Scaling of the data is critical for optimal visualization and clustering. 
To understand how to set scales in FlowJo, see the instruction file "Setting axis in FlowJo for Aurora data.pptx". Use the channel values format.
To export your data in CSV format, preserving the transformations from FlowJo, see the instructions in "Exporting data in csv format.PNG". For more 
detail, see https://docs.flowjo.com/flowjo/graphs-and-gating/gw-transform-overview/

The CSV-based version incorporates EmbedSOM for fast parallelized SOM clustering and visualization. There is also an option to use the Irish lab's 
T-REX method as an add-on for visualizing and identifying changed regions in the tSNE, UMAP or EmbedSOM plots.
https://elifesciences.org/articles/64653

Alternatively, you may use the legacy version "flow analysis", which starts with FCS files.

For all flow cytometry analysis, we recommend pre-gating and exporting the populations of cells you are interested in. While exporting, 
adding group or variable tags to the file names will help you sort the files with the script. For a tutorial on using the script, see the lab's website:
http://www.listonlab.uk/cross-entropy-test/

Publication here:
https://www.cell.com/cell-reports-methods/pdfExtended/S2667-2375(22)00295-8

To install all the packages and dependencies for running the script, run flowcytoscript_setup.R

The parameter file (analyze_flow_cytometry_parameter_csv.R) contains all the variables you may wish to modify, 
so no changes should be necessary in the main script (analyze_flow_cytometry_csv.R) or source (00_src) files.
