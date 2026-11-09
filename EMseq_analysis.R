# H3K27ac ChIP-EM-seq analysis summary
# Final use: Figure 4F-G
# Detailed implementation: ./analysis_code/ (numbered in workflow order).

library(data.table)
library(ggplot2)

PROJECT_ROOT <- "RKO_TAZ_v_VAL_project"

# VAL acEM-seq methylation was averaged across two biological replicates.
# Signal was summarized around final ATAC peak centers. Peaks lacking adequate
# positional representation were retained in recovery summaries but excluded
# from quantitative methylation comparisons.

# Final waterfall display:
# - +/-1.5 kb around each ATAC peak center;
# - 500-bp display bins;
# - rows sorted from highest to lowest mean VAL methylation;
# - methylation scale fixed from 0 to 1;
# - missing values shown with a contrasting neutral color.

# Under the established coverage criterion, H3K27ac-associated fragments were
# classified as methylated when their mean VAL acEM-seq methylation was >=50%.
# The adjacent status bar in Figure 4G shows this binary classification.

