Sweave("PatuS_countsummary.Rnw");
library(tools);

texi2dvi("PatuS_countsummary.tex",pdf=TRUE);

