Sweave("PatuT_countsummary.Rnw");
library(tools);

texi2dvi("PatuT_countsummary.tex",pdf=TRUE);

