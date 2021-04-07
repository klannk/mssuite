# Documentation

## Defaults (class)
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
* `list_of_df`:*list* or *array* 
Array that contains the dataframes to be processed.
* `function`: *function* 
Function that is applied to dataframes. Needs to return a dataframe
* args and kwargs for the function

* **Returns**: Array of processed dataframes

```python
Defaults.get_channels(self, input_file, custom=None)
```
Function that returns all intensity/abundance containing column names of the input dataframe
* `input_file`: *dataframe*
Dataframe containing protein/peptide level quantifications
* `custom`: *string*,optional
If intensity/abundance column is not from default PD output, it can be set here.

* **Returns**: Array of column names

## Preprocessing
