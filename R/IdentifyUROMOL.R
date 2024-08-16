#' UROMOL2021 class of non-muscle-invasive bladder cancer (NMIBC) samples identified by methylation values
#'
#' This package performs a Pearson nearest-centroid classification according to the UROMOL2021 classes of NMIBC based on methylation values from Infinium HumanMethylation450 BeadChip (HM450K)
#'
#' @param x data.frame with CpGs from HM450K or Infinium MethylationEPIC (Epic) in rows and samples to be classified in columns (accepting a dataframe with a single column namely one sample).
#'
#' @param minCor numeric value specifying a minimal threshold for best Pearson's correlation between a sample's methylation profile and centroids profiles.
#' A sample showing no correlation above this threshold will remain unclassifed and prediction results will be set to NA. Default minCor value is 0.2.
#'
#' @return dataframe with classification results. The UROMOL_2021_class column returns the predicted class labels of the samples.
#' The cor_pval column returns the p-value associated with the Pearson's correlation between the sample and the nearest centroid.
#' The separationLevel gives a measure (ranging from 0 to 1) of how a sample is representative of its consensus class,
#' with 0 meaning the sample is too close to the other classes to be confidently assigned to one class label, and 1 meaning the sample is very representative of its class.
#' This separationLevel is measured as follows: (correlation to nearest centroid - correlation to second nearest centroid) / median difference of sample-to-centroid correlation.
#' The remaining four columns return the Pearson's correlation values for each sample and each centroid.
#' UROMOL_2021_class predictions are set to NA if the minCor condition is not verified.
#'
#' @export
#'
#' @examples
#' data(test_data)
#' IdentifyUROMOL(test_data, minCor = .2)
#'
#' @note This classifier was built referenced to the published tool for the consensus classes of NMIBC:
#' Lindskrog, S.V., Prip, F., Lamy, P. et al. An integrated multi-omics analysis identifies prognostic molecular subtypes of non-muscle-invasive bladder cancer. Nat Commun 12, 2301 (2021). https://doi.org/10.1038/s41467-021-22465-w


IdentifyUROMOL <- function(x, minCor = .2) {
  # load centroids_UROMOL
  data(centroids_UROMOL)

  classes <- c("Class_1", "Class_2a", "Class_2b", "Class_3")

  gkeep <- intersect(rownames(centroids_UROMOL), rownames(x))
  if (length(gkeep) == 0) stop("Empty intersection between profiled CpGs and the CpGs used for classification.\n Make sure that CpGs names correspond to the type of identifiers specified by Illumina Methylation Array HM450K/Epic")
  if (length(gkeep) < 0.8 * nrow(centroids_UROMOL)) warning("Input methylation profiles include less than 80% of the CpGs used for classification. Results may not be relevant")

  col_N <- colnames(x)
  x <- as.data.frame(x[gkeep, ])
  colnames(x) <- col_N
  rownames(x) <- gkeep
  cor.dat <- as.data.frame(cor(x, centroids_UROMOL[match(gkeep, rownames(centroids_UROMOL)), classes], use = "complete.obs", method = "pearson"), row.names = colnames(x))

  # Best correlated centroid
  cor.dat$nearestCentroid <- apply(cor.dat, 1, function(y) {
    classes[which.max(y)]
  })
  cor.dat$corToNearest <- apply(cor.dat[, classes], 1, max)
  cor.dat$cor_pval <- sapply(colnames(x), function(smp) {
    cor.test(x[gkeep, smp], centroids_UROMOL[match(gkeep, rownames(centroids_UROMOL)), cor.dat[smp, "nearestCentroid"]])$p.value
  })

  # Separation level metrics
  cor.dat$deltaSecondNearest <- apply(cor.dat$corToNearest - cor.dat[, classes], 1, function(x) {
    sort(x)[2]
  })
  cor.dat$deltaMed <- apply(cor.dat$corToNearest - cor.dat[, classes], 1, median)
  cor.dat$separationLevel <- cor.dat$deltaSecondNearest / cor.dat$deltaMed

  cor.dat$UROMOL_2021_class <- cor.dat$nearestCentroid

  # Set to NA if best correlation < minCor
  for (i in 1:length(rownames(cor.dat))) {
    if (cor.dat[i, "corToNearest"] < minCor) {
      cor.dat[i, "UROMOL_2021_class"] <- NA
      cor.dat[i, "separationLevel"] <- NA
    }
  }

  cor.dat <- cor.dat[, c("UROMOL_2021_class", "cor_pval", "separationLevel", classes)]
  return(cor.dat)
}
