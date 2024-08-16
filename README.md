# IdentifyNMIBC

This package provides two Pearson nearest-centroid classifications based on DNA methylation array data:

1)  *IdentifyMethyl* intended for assigning class labels to single samples according to the three BCPP methylation subgroups of non-muscle-invasive bladder cancer (NMIBC): M1, M2 and M3.

2)  *IdentifyUROMOL* intended for assigning class labels to single samples according to the four UROMOL2021 classes of non-muscle-invasive bladder cancer (NMIBC): Class 1, Class 2a, Class 2b and Class 3.

## Citation

Please cite xxx. DOI: xxx

## Install

You may install this package with devtools:

```{r}
library(devtools)

devtools::install_github("yumh0811/IdentifyNMIBC", build_vignettes = TRUE)

library(IdentifyNMIBC)
```

## Usage

1.  **IdentifyMethyl**

```{r}
IdentifyMethyl(x, minCor = .2)
```

`x`: data.frame with CpGs from HM450K or Infinium MethylationEPIC (Epic) in rows and samples to be classified in columns (accepting a dataframe with a single column namely one sample).

`minCor`: numeric value specifying a minimal threshold for best Pearson's correlation between a sample's methylation profile and 3 centroids profiles. A sample showing no correlation above this threshold will remain unclassifed and prediction results will be set to NA. Default minCor value is 0.2.

2.  **IdentifyUROMOL**

```{r}
IdentifyUROMOL(x, minCor = .2)
```

`x`: data.frame with CpGs from HM450K or Infinium MethylationEPIC (Epic) in rows and samples to be classified in columns (accepting a dataframe with a single column namely one sample).

`minCor`: numeric value specifying a minimal threshold for best Pearson's correlation between a sample's methylation profile and 4 centroids profiles. A sample showing no correlation above this threshold will remain unclassifed and prediction results will be set to NA. Default minCor value is 0.2.

## Example

1.  **IdentifyMethyl**

```{r}
data(test_data)

Methylation_subgroups <- IdentifyMethyl(test_data, minCor = .2)

head(Methylation_subgroups)
```

|         | Methylation_subgroups | cor_pval | separationLevel | M1        | M2        | M3        |
|----------|------------|----------|----------|----------|----------|----------|
| sample1 | M3                    | 0        | 1               | 0.9689645 | 0.9255264 | 0.9773381 |
| sample2 | M2                    | 0        | 1               | 0.9434526 | 0.9652757 | 0.9439114 |
| sample3 | M1                    | 0        | 1               | 0.9378622 | 0.9006545 | 0.9147292 |
| sample4 | M3                    | 0        | 1               | 0.9648775 | 0.9311472 | 0.9806707 |

The classifier returns a data.frame with 6 columns:

`Methylation_subgroups`: the predicted class labels of the samples.

`cor_pval`: the p-value associated with the Pearson's correlation between the sample and the nearest centroid.

`separationLevel`: gives a measure (ranging from 0 to 1) of how a sample is representative of its consensus class, with 0 meaning the sample is too close to the other classes to be confidently assigned to one class label, and 1 meaning the sample is very representative of its class. It is measured as follows: (correlation to nearest centroid - correlation to second nearest centroid) / median difference of sample-to-centroid correlation.

`M1` `M2` `M3`: The remaining three columns return the Pearson's correlation values for each sample and each centroid.

2.  **IdentifyUROMOL**

```{r}
data(test_data)

UROMOL2021_class <- IdentifyUROMOL(test_data, minCor = .2)

head(UROMOL2021_class)
```

|         | UROMOL_2021_class | cor_pval | separationLevel | Class_1   | Class_2a  | Class_2b  | Class_3   |
|---------|---------|---------|---------|---------|---------|---------|---------|
| sample1 | Class_1           | 0        | 0.9511100       | 0.9808939 | 0.9428596 | 0.9711279 | 0.9701239 |
| sample2 | Class_2a          | 0        | 0.9913753       | 0.9492292 | 0.9687753 | 0.9488891 | 0.9433909 |
| sample3 | Class_2b          | 0        | 0.9940678       | 0.9142828 | 0.9159124 | 0.9374620 | 0.9161666 |
| sample4 | Class_3           | 0        | 0.2628464       | 0.9799166 | 0.9481460 | 0.9663972 | 0.9823269 |

The classifier returns a data.frame with 7 columns:

`UROMOL_2021_class`: the predicted class labels of the samples.

`cor_pval`: the p-value associated with the Pearson's correlation between the sample and the nearest centroid.

`separationLevel`: gives a measure (ranging from 0 to 1) of how a sample is representative of its consensus class, with 0 meaning the sample is too close to the other classes to be confidently assigned to one class label, and 1 meaning the sample is very representative of its class. It is measured as follows: (correlation to nearest centroid - correlation to second nearest centroid) / median difference of sample-to-centroid correlation.

`Class_1` `Class_2a` `Class_2b` `Class_3`: The remaining three columns return the Pearson's correlation values for each sample and each centroid.
