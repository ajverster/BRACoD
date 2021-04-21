
library('reticulate')
source_python("/mnt/Freya/Analysis/Microbiome/Lib/SVSS.py")
infile <- "/mnt/Freya/Analysis/Microbiome/Data/Obesity_OTUCounts.csv"
Df.Microbiome <- read.table(infile, sep = "\t", header = TRUE)

X <- Df.Microbiome[,which(colnames(Df.Microbiome) == "Otu00001"):which(colnames(Df.Microbiome) == "Otu00821")]
Y <- Df.Microbiome[,"butyric"]

r <- remove_null(X,Y)
X <- r[[1]]
Y <- r[[2]]

X.Norm = scale_counts(X, 1000, 100)
trace = run_svss(X.Norm, Y, n_sample = 5000,njobs=4)

# Test if the important parameters in the model have converged properly
convergence_tests(trace, inclusion_cutoff=0.1)

df_summary = summarize_trace(trace, colnames(X.Norm), inclusion_cutoff=0.1)
