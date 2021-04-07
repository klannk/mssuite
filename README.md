# mssuite 

The mssuite python package provides a framework for streamlined data analysis of TMT based shotgun proteomics data. The package is of modular nature and besides already implemented pipelines you can easily build your one data analysis pipeline, including automated plots. The package defaults to PSM/Peptide output files from ProteomeDiscoverer Software. However, it is easily customizable to your input files, by changing a few paramters. 

## Install

### Using pip
mssuite has been uploaded to the PyPi repository and can be easily installed using the **pip** package manager:
```python
    pip install mssuite
```
Alternatively you can download the binaries and install locally using the following command:
```python
    !cd /PATH/TO/PACKAGE
    pip install .
```
### Compile from source
To compile from source, please download the package files and compile using python:
```python
    !cd /PATH/TO/PACKAGE
    python3 -m build
```
## Quickstart
Export your proteomics experiment on PSM or Peptide level from ProteomeDiscoverer Software as a tab-delimited text file. mssuite works with pandas dataframes, so you need to load your data as a pandas dataframe:
```python
import pandas as pd

psms = pd.read_csv("PATH/TO/FILE.txt",sep='\t',header=0)
```
To calculate the proper statistics you need to specify the experimental conditions in the order they appear in your input file:
```python
    conditions = ['Control','Control','Control','Treatment','Treatment','Treatment']
```
Lastly, you need to specify a working directory, where your output files are written:
```python
    wd = 'YOUR/PATH'
```
To start the analysis you initialize the Pipelines module from the mssuite package and run the analysis:
```python
pipe = Pipelines()
results = pipe.singlefile_lmm(psms, conditions, wd=wd)
```
The pipeline will now filter your input for contaminants and razor peptides, normalize the data, calculate differential expression analysis for all condition pairs by using a peptide-based linear mixed regression model and create automated plots (Clustermap, Abundances and Volcano plots). All output will be written to your specified working directory. The resulting dataframe will be written to a csv file and also returned to the results variable (see code above) for further use if intended.

## Credits and Citation


Link to Publication

## License

MIT License

Copyright (c) 2021 Kevin Klann

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.