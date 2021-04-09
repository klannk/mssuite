# Documentation

## Defaults 

The mssuite package locates the needed columns in the input dataframe by a string/label based approach. It searches the column names for matching strings and uses them in the downstream analysis. This object contains PD specific default names for columns. The other modules in the package will access this object to derive the needed column names.

### Class variables

```python
Defaults.MasterProteinAccession = "Master Protein Accessions"
```

This variable contains the default string that is used to extract the protein identifiers used in all downstream analysis.

```python
Defaults.labelsForLMM = ['Annotated Sequence','Master Protein Accessions','Abundance:']
```

This variable contains the needed column labels for the peptide based models. If you change it, please provide the names in the following order: 1. Peptide based identifier (Sequence or something similar), 2. Identifier on protein level (e.g. UniProt Accession or Gene Name) and 3. A string that is shared by all intensity/abundance containing columns.

```python
Defaults.AbundanceColumn = "Abundance:"
```

This variable contains the string that is shared by all intensity/abundance containing columns i.e. the sample quantifications.

```python
Defaults.file_id = "File ID"
```

If your PSM file contains PSMs and quantifications from multiple files that need to be splitted then this variable contains the column name string for the file ID. Is used by `Preprocessing.psm_splitting`.

### Methods

```python
Defaults.processor(self,list_of_df,function,*args,**kwargs)
```

Processor function for batch processing an array of dataframes.

- `list_of_df`:_list_ or _array_
  Array that contains the dataframes to be processed.
- `function`: _function_
  Function that is applied to dataframes. Needs to return a dataframe
- args and kwargs for the function

- **Returns**: Array of processed dataframes

```python
Defaults.get_channels(self, input_file, custom=None)
```

Function that returns all intensity/abundance containing column names of the input dataframe

- `input_file`: _dataframe_
  Dataframe containing protein/peptide level quantifications
- `custom`: _string_,optional
  If intensity/abundance column is not from default PD output, it can be set here.

- **Returns**: Array of column names

## Preprocessing

### Methods

```python
Preprocessing.filter_peptides(self,input_file)
```

Function that filters peptides/PSMs for belonging to contaminants and removes razor peptides.
Only possible for PD inputs at the current stage.

- `input_file`: _dataframe_
  Dataframe containing the peptides/PSMs.
- **Returns**: _dataframe_

Uses internally the `Defaults.MasterProteinAccession` variable.

```python
Preprocessing.psm_splitting(self, input_file)
```

Only possible for PD inputs at the current stage.

ProteomeDiscoverer outputs contain all PSMs in one large table long-format without splitting it into multiple files, if experiments have been analysed in the same workflow (e.g. replicates). This makes normalisation and preprocessing difficult. Therfore this function splits the large PSM table into multiple dataframes and returns an array of dataframes.

- `input_file`: _dataframe_
  Dataframe containing the PSMs.
- **Returns**: Array of dataframes

Uses internally the `Defaults.file_id` variable.

```python
Preprocessong.psm_joining(self, input_list)
```

Only possible for PD inputs at the current stage.

This function joins different PSM files that have been split by `Preprocessing.psm_splitting` function and processed further. It joins the provided dataframes in a wide format compatible with downstream analysis such as IRS normalisation.

- `input_list`: _array of dataframes_ Dataframes to join
- **Returns**: _dataframe_

```python
Preprocessing.median_normalisation(self, input_file, channels)
```

Performs median normalisation (column-wise) for all `channels` (the columns to use).

- `input_file`:_dataframe_ Dataframe to be normalised
- `channels`:_array_ of column names to use for normalisation. To get channel array use `Defaults.get_channels` function
- **Returns**: _dataframe_

```python
Preprocessing.total_intensity(self, input_file, channels)
```

Performs total intensity normalisation (column-wise) for all `channels` (the columns to use).

- `input_file`:_dataframe_ Dataframe to be normalised
- `channels`:_array_ of column names to use for normalisation. To get channel array use `Defaults.get_channels` function
- **Returns**: _dataframe_

```python
Preprocessing.TMM(self, input_file, channels)
```

Performs TMM normalisation (column-wise) for all `channels` (the columns to use).
Robinson, M.D., Oshlack, A. A scaling normalization method for differential expression analysis of RNA-seq data. Genome Biol 11, R25 (2010). https://doi.org/10.1186/gb-2010-11-3-r25

- `input_file`:_dataframe_ Dataframe to be normalised
- `channels`:_array_ of column names to use for normalisation. To get channel array use `Defaults.get_channels` function
- **Returns**: _dataframe_

```python
Preprocessing.chunks(self,l,n)
```

Yield successive n-sized chunks from. Used internally by `Preprocessing.IRS_normalisation`.

```python
Preprocessing.IRS_normalisation(self, input_file, bridge, plexes)
```

This function performs IRS normalisation for a input pandas df. Bridge channels have to be same TMT channel in all multiplexes and plexes must have same size.
Plubell DL, Wilmarth PA, Zhao Y, Fenton AM, Minnier J, Reddy AP, Klimek J, Yang X, David LL, Pamir N. Extended Multiplexing of Tandem Mass Tags (TMT) Labeling Reveals Age and High Fat Diet Specific Proteome Changes in Mouse Epididymal Adipose Tissue. Mol Cell Proteomics. 2017 May;16(5):873-890. doi: https://10.1074/mcp.M116.065524. Epub 2017 Mar 21. PMID: 28325852; PMCID: PMC5417827

- `input_file`: _dataframe_ Dataframe to be normalised, containing multiple TMT runs.
- `bridge`: _str_ Channel to use as bridge channel between multiplexes. e.g. '127N'.
- `plexes`: _int_ Number of multiplexes in dataframe.
- **Returns**: _dataframe_

## Annotation

### Methods

```python
Annotation.basic_annotation(self, input_file)
```

Performs basic annotation by adding Gene Symbols, Protein Names, Taxonomy and Molecular Weight to dataframe. Intended to use as last step of analysis pipeline. Dataframes index needs to be UniProt Accession. Works currently only for human proteins.

- `input_file`: _dataframe_ Dataframe to be annotated. Index should be Accession.
- **Returns**: _dataframe_

## Rollup

This class contains functions for calculating protein level quantifications from Peptide level data.

### Methods

```python
Rollup.protein_rollup_sum(self,input_file, channels)
```

This function takes Peptide level (or PSM) dataframes and performs a sum based rollup to protein level.
the channels variable takes an array of column names that contain the quantifictions. You can create such an array via `Defaults.get_channels` function.

- `input_file`: _dataframe_ Peptide level dataframe
- `channels`: _array_ of column names to use for rollup. To get channel array use `Defaults.get_channels` function
- **Returns**: _dataframe_

```python
Rollup.protein_rollup_median(self,input_file, channels)
```

This function takes Peptide level (or PSM) dataframes and performs a median based rollup to protein level.
the channels variable takes an array of column names that contain the quantifictions. You can create such an array via `Defaults.get_channels` function.

- `input_file`: _dataframe_ Peptide level dataframe
- `channels`: _array_ of column names to use for rollup. To get channel array use `Defaults.get_channels` function
- **Returns**: _dataframe_

```python
Rollup.protein_rollup_mean(self,input_file, channels)
```

This function takes Peptide level (or PSM) dataframes and performs a mean based rollup to protein level.
the channels variable takes an array of column names that contain the quantifictions. You can create such an array via `Defaults.get_channels` function.

- `input_file`: _dataframe_ Peptide level dataframe
- `channels`: _array_ of column names to use for rollup. To get channel array use `Defaults.get_channels` function
- **Returns**: _dataframe_

## HypothesisTesting

### Class variables

After initialisation contains 2 variables:

- `HypothesisTesting.pair_names`: _array_ contains all performed comparisons after differential expression analysis.
- `HypothesisTesting.pair_names`: _dict_ contains column names for p_values, q_values and fold changes after differential expression analysis.

### Methods

```python
HypothesisTesting.t_test(self, input_data, matrix1, matrix2, name='')
```

Calculates p values and corrected p values (q values, BH-FDR) according to a students t-test (two-sided) for each row/protein in a dataframe. Needs column names for the two matrices as arrays.
Name defines the suffix that is added to the resulting columns e.g. comparison.

- `input_data`: _dataframe_ contains the quantification data
- `matrix1`: _array_ of column names that make the first population
- `matrix2`: _array_ of column names that make the second population
- `name`: _str_, optional. Default = ''. Adds the name as suffix to the resulting columns to distinguish multiple tests-
- **Returns**: _dataframe_

```python
HypothesisTesting.tessa(self,source)
```

Returns all possible pairs from source as array

- `source`: _array_ with elements
- **Returns**: nested array with all pairs

```python
HypothesisTesting.peptide_based_lmm(self, input_file, conditions, norm=Preprocessing.total_intensity,pairs=None)
```

Performs gene-wise linear mixed modelling for differential expression on peptide level.

- `input_file`: _dataframe_ Input data
- `conditions`: _array_ of treatment conditions matching the order of the quantification columns. E.g. `['Control','Treatment','Control','Treatment']`
- `norm`: _function_, normalisation function of the Preprocessing module. If no normalisation should be applied set to `None`
- `pairs`: _array_, optional. Nested array containing pairs that should be tested for differential analysis. E.g.`[['Control','Treatment1'],['Control','Treatment2']]`. If None, performs analysis for all possible pairs.
- **Returns**: _dataframe_

```python
HypothesisTesting.get_comparisons(self)
```

Returns all comparisons that have been performed on that dataframe for further use.

- **Returns**: `HypothesisTesting.pair_names`

```python
HypothesisTesting.export_comparison_strings(self)
```

Returns a nested dictionary (json format) with all comparisons tested during hypothesis testing as keys and the column names for P values, q values and fold changes.
Automatically used when performing t-test or peptide_based_lmm

- **Returns**: _dict_

```python
HypothesisTesting.get_columnnames_for_comparison(self, comparison)
```

Extracts the column names from `HypothesisTesting.comparison_data` and returns column names for fold change, p value and q value.

- `comparison`: _str_ comparison to return. Element of `HypothesisTesting.pair_names` array.
- **Returns**: fold_change: _str_, p_value: _str_, q_value: _str_ (column names)

```python
HypothesisTesting.get_significant_hits(self, input_data, comparison, fc_cutoff=0.5, p_cutoff=0.05, use_q=True)
```

Creates a _dict_ containing upregulated and downregulated genes from comparison.

- `input_data`:_dataframe_ that contains columns with differential expression analysis of comparison pair
- `comparison`: _str_ comparison to return. Element of `HypothesisTesting.pair_names` array.
- `fc_cutoff`: _float_ fold change cutoff used to filter significant hits
- `p_cutoff`:_float_ p value cutoff used to filter significant hits
- `use_q`: _bool_ If True uses q_value instead of p value as cutoff for filtering
- **Returns**: _dict_

## Visualisation

This class contains functions to generate automated plots in a data analysis pipeline.

### Methods

```python
Visualization.volcano_plot(self, input_file, fold_change, pval, comparison, wd)
```

Produces a volcano plot (plotting fold changes against p-values) for a comparison and saves it to the workind directory

- `input_file`:_dataframe_ containing differential expression data
- `fold_change`:_str_ column name of column that contains the fold change values. Extract column names by `HypothesisTesting.get_columnnames_for_comparison`.
- `pval`:_str_ column name of column that contains the p values. Extract column names by `HypothesisTesting.get_columnnames_for_comparison`.
- `comparison`: _str_ Element of `HypothesisTesting.pair_names`.
- `wd`:_str_ Working directory path, where output will be saved.

```python
Visualization.boxplots(self,input_file,channels,wd)
```

Produces a boxplot of all quantification channels for quality control and saves it in working directory.

- `input_file`: _dataframe_ that contains the quantification values.
- `channels`: _array_ of column names that contain the quantification values. Generate by `Defaults.get_channels`.
- `wd`:_str_ Working directory path, where output will be saved.

```python
Visualization.heatmap(self, input_file,channels,conditions,wd)
```

Produces a clustered heatmap with conditions as labels and saves it to the working directory.

- `input_file`:_dataframe_ containing the quantification values.
- `channels`: _array_ of column names that contain the quantification values. Generate by `Defaults.get_channels`.
- `conditions`: _array_ of treatment conditions matching the order of the quantification columns. E.g. `['Control','Treatment','Control','Treatment']`.
- `wd`:_str_ Working directory path, where output will be saved.

## Pathway Enrichment

This class performs Reactome pathway enrichment for gene sets, tested either against a custom background or the genome.

### Class Variables

- `PathwayEnrichment.database` - contains the filtered reactome database after `PathwayEnrichment.get_background_sizes` or `PathwayErichment.get_pathway_sizes`.
- `PathwayEnrichment.counts` - contains the occurances of each pathway in the background
- `PathwayEnrichment.total` - number of background genes

### Methods

```python
PathwayEnrichment.get_background_sizes(self, background=list)
```

Calculates the occurances of pathways in a custom background list and stores them to the `PathwayEnrichment.total` variable.

- `background`:_list_ of genes used as background

```python
PathwayEnrichment.get_pathway_sizes(self, species='Homo sapiens')
```

If no custom background is provided then all Genes belonging to `species` are used as background. Calucluates the occurances of each pathway in the species wide background.

- `species`:_str_ with species name to use.

```python
PathwayEnrichment.get_enrichment(self, genes)
```

After background occurences have been calculated by either `PathwayEnrichment.get_background_sizes` or `PathwayEnrichment.get_pathway_sizes`, enrichment calculation will be performed with a hypergeometric test testing the provided genes list against the background for enrichment.

- `genes`:_list_ of Accessions to test for enrichment. Have to be UniProt Accessions
- **Returns**: _dataframe_ containing pathways and enrichment statistics.

## Pipelines

This class contains prebuilt pipelines, that streamline data analysis starting from a PSM/Peptide file, performing processing and normalisation, differential expression analysis and automated plots.

### Methods

```python
Pipelines.singlefile_lmm(self, psms, conditions, pairs=None, wd=None,filter=True)
```

This pipeline takes a single PSM file as an input and performs complete data analysis including: Filtering, Normalisation, differential expression analysis and plotting.

- `psms`:_dataframe_ containing the PSM level data from PD
- `conditions`: _array_ of treatment conditions matching the order of the quantification columns. E.g. `['Control','Treatment','Control','Treatment']`.
- `pairs`: _array_, optional. Nested array containing pairs that should be tested for differential analysis. E.g.`[['Control','Treatment1'],['Control','Treatment2']]`. If None, performs analysis for all possible pairs.
- `wd`:_str_ Working directory path, where output will be saved.
- `filter`:_bool_ If True performs contaminant filtering. Only for PD input files.

```python
Pipelines.multifile_lmm(self, psms, conditions, bridge, pairs=None, wd=None,filter=True)
```

This pipeline takes a single PSM file with multiple multiplexes analysed together as an input and performs complete data analysis including: Filtering, Normalisation,IRS, differential expression analysis and plotting.

- `psms`:_dataframe_ containing the PSM level data from PD
- `conditions`: _array_ of treatment conditions matching the order of the quantification columns. E.g. `['Control','Treatment','Control','Treatment']`.
- `bridge`: _str_ Channel to use as bridge channel between multiplexes. e.g. '127N'.
- `pairs`: _array_, optional. Nested array containing pairs that should be tested for differential analysis. E.g.`[['Control','Treatment1'],['Control','Treatment2']]`. If None, performs analysis for all possible pairs.
- `wd`:_str_ Working directory path, where output will be saved.
- `filter`:_bool_ If True performs contaminant filtering. Only for PD input files.
