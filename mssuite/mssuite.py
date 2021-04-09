'''
# mssuite - implementation of PBLMM algorithm and streamlined data analysis
# Copyright (C) 2021 Kevin Klann
#This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
'''
import os
import warnings
from collections import Counter

import DynaTMT.DynaTMT as mePROD
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.cluster.hierarchy import dendrogram, distance, fcluster, linkage
from scipy.spatial import distance_matrix
from scipy.stats import hypergeom, trim_mean, ttest_ind
from statsmodels.formula.api import ols
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")
mpl.style.use('tableau-colorblind10')
dirname = os.path.dirname(__file__)

class Defaults:
    '''
    This object contains PD specific default names for columns. Most functions will access this by default, but
    can be set manually. To use this package with input files other than Proteome Discoverer, initialise a instance of Defaults and set
    the column names according to your data layout.
    '''
    MasterProteinAccession = "Master Protein Accessions"
    labelsForLMM = [
        'Annotated Sequence',
        'Master Protein Accessions',
        'Abundance:',
    ]
    AbundanceColumn = "Abundance:"
    file_id = "File ID"

    def processor(self,list_of_df,function, *args,**kwargs):
        '''Processor function that applies a certain function to a list of dataframes, to allow rapid batch processing.
        Returns a list of processed dataframes
        '''
        results=[]
        for count, value in enumerate(list_of_df):
            results.append(function(value, *args,**kwargs))
        return results

    def get_channels(self, input_file, custom=None):
        '''Returns an array of all column names where Abundances are stored. Accesses Defaults object, but also custom string can be applied. However, its recommended
        to change defaults.AbundanceColumn for compatibility with all other functions. Its basically just a wrapper for python list comprehension.
        '''
        if custom is None:
            channels = [col for col in input_file.columns if self.AbundanceColumn in col]
        else:
            channels = [col for col in input_file.columns if custom in col]
        return channels

class Preprocessing:
    def __init__(self):
        pass

    def filter_peptides(self, input_file):
        '''
        Filters peptide files for non-unique peptides and contaminant proteins
        '''
        mpa1=Defaults.MasterProteinAccession
        mpa = [col for col in input_file.columns if mpa1 in col]
        mpa = mpa[0]
        input_file = input_file[~input_file[mpa].str.contains(';', na=False)]
        input_file = input_file[input_file['Contaminant'] == False]
        return input_file

    def psm_splitting(self, input_file):
        '''takes a file of combined PD analysis and splits it into separate dfs for normalization and processing
        '''
        try:
            input_file[Defaults.file_id] = input_file[Defaults.file_id].str.split('.')[0] #split filenames, since everything after the dot are fraction numbers
        except ValueError:
            pass
        grouped_input = input_file.groupby(by=Defaults.file_id)#groupby files
        
        arrayofdataframes = [grouped_input.get_group(x) for x in grouped_input.groups]#split dataframe into list of dataframes 
        return arrayofdataframes #Returns: Array of dataframes 
    
    def psm_joining(self, input_list):
        '''Joins all dataframes in a list for IRS normalisation. 
        '''
        for idx in range(0,len(input_list)): #reidexing on sequence and modifications
            print(idx)
            
            #print(input_list[idx].index)
            channels = [col for col in input_list[idx].columns if "Abundance:" in col]
            input_list[idx]['Identifier'] = input_list[idx]['Annotated Sequence'].map(str) + input_list[idx]['Modifications']

            input_list[idx][channels] = input_list[idx].groupby(['Identifier'])[channels].transform('sum')
            input_list[idx] = input_list[idx].drop_duplicates(subset=['Identifier'])
            input_list[idx].index = input_list[idx]['Identifier']
            input_list[idx]=input_list[idx].add_suffix(idx) #adds to all columns the number of dataset to avoid duplicates
            #print(input_list[idx].index)
        
        Input1=input_list[0].join(input_list[1:],how='inner')
        print("Done")
        return Input1

    def median_normalisation(self, input_file, channels):
        '''
        #Performs Median normalisation. Besides input file the function needs an array of all
        column names that contain the quantifications to be normalized (channels).
        '''
        input_file = input_file.dropna(subset=channels)
        print("Median Normalization")
        minimum = np.argmin(input_file[channels].median().values)
        summed = np.array(input_file[channels].median().values)
        minimum = summed[minimum]
        norm_factors = summed/minimum
        input_file[channels] = input_file[channels].divide(
            norm_factors, axis=1)
        print("Normalization done")
        return input_file

    def total_intensity(self, input_file, channels):
        '''
        #Performs total intensity normalisation. Besides input file the function needs an array of all
        column names that contain the quantifications to be normalized (channels).
        '''
        # remove missing value rows
        input_file = input_file.dropna(subset=channels)
        print("Normalization")
        # calculate summed intensity for each column and search minimum index
        minimum = np.argmin(input_file[channels].sum().values)
        summed = np.array(input_file[channels].sum().values)
        minimum = summed[minimum]
        # calculuate norm factors
        norm_factors = summed/minimum
        # normalize input with norm factors
        input_file[channels] = input_file[channels].divide(
            norm_factors, axis=1)
        print("Normalization done")
        return input_file

    def TMM(self, input_file, channels):
        '''This function implements TMM normalisation (Robinson & Oshlack, 2010, Genome Biology).
        '''
        input_file = input_file.dropna(subset=channels)
        # Trim the 5% Quantiles from dataset
        input_trim = input_file[input_file[channels]
                                < input_file[channels].quantile(.95)]
        print("TMM Normalization")
        # Calculate ratios to first channel
        input_trim[channels] = input_trim[channels].divide(
            input_trim[channels[0]], axis=0)
        tm = np.argmin(trim_mean(input_trim[channels], 0.25))
        summed = np.array(trim_mean(input_trim[channels], 0.25))
        minimum = summed[tm]
        norm_factors = summed/minimum
        input_file[channels] = input_file[channels].divide(
            norm_factors, axis=1)
        return input_file

    def chunks(self, l, n):
        """Yield successive n-sized chunks from l."""
        for i in range(0, len(l), n):
            yield l[i:i + n]

    def IRS_normalisation(self, input_file, bridge, plexes):
        '''
        This function performs IRS normalisation for a input pandas df. Bridge channels have to be same TMT channel and plexes must have same size
        bridge = String that defines bridge channel
        plexes = Number (INT) of plexes to normalize
        quant = String that is included in all Quantification columns
        '''
        print('Internal Reference scaling')
        abundance_column=Defaults.AbundanceColumn
        # search for quantification columns
        defaults=Defaults()
        channels = defaults.get_channels(input_file=input_file, custom=abundance_column)
        # remove missing values from input
        input_file = input_file.dropna(subset=channels)
        # search for bridge channels
        bridge_channels = [i for i in channels if str(bridge) in i]
        # calculate mean of all bridge channels
        bridge_mean = np.array(input_file[bridge_channels].mean(axis=1))
        # calculate norm factors for bridge
        cfs = input_file[bridge_channels].divide(bridge_mean, axis=0)
        cfs_cols = cfs.columns
        chunks_channels = list(self.chunks(
            channels, int(len(channels)/plexes)))
        for i in range(0, len(chunks_channels), 1):
            norms = np.array(cfs[cfs_cols[i]])
            input_file[chunks_channels[i]
                       ] = input_file[chunks_channels[i]].divide(norms, axis=0)
        input_file = input_file.dropna(subset=channels)
        print('Done')
        return input_file

class Annotation:
    # all functions related to annotation of protein/peptide files
    def __init__(self):
        pass

    def basic_annotation(self, input_file):
        '''
        Performs basic annotation and adds Gene symbols, protein names, taxonomy and MW to the df. The
        DFs index needs to be the accession. So far only for human proteins.
        '''

        # Reads annotation file

        database = pd.read_csv(os.path.join(dirname, '../data/Annotation_data_human.txt'), sep='\t', header=0)
        database.index = database['Entry']

        # Iterate over entries in input file and write annotations
        for index in input_file.index:
            try:
                input_file.loc[index,
                               'Gene_Symbol'] = str(database.loc[index, 'Gene names']).split(' ')[0]
                input_file.loc[index,
                               'Protein_Name'] = database.loc[index, 'Protein names']
                input_file.loc[index,
                               'Organism'] = database.loc[index, 'Organism']
                input_file.loc[index, 'Mass'] = database.loc[index, 'Mass']
            except KeyError:
                pass
        return input_file

class Rollup:
    def __init__(self):
        pass

    def protein_rollup_sum(self, input_file, channels):
        '''
        This function takes Peptide level (or PSM) dataframes and performs a sum based rollup to protein level.
        the channels variable takes an array of column names that contain the quantifictions. You can create such an
        array via this command:
        channels = [col for col in PSM.columns if 'Abundance:' in col]

        mpa1 variable contains a string that is included in the Accession column. The function will search for the column containing the string
        and use it for rollup.

        Returns Protein level DF.
        '''
        mpa1=Defaults().MasterProteinAccession
        print('Calculate Protein quantifications from PSM')
        mpa = [col for col in input_file.columns if mpa1 in col]
        mpa = mpa[0]

        PSM_grouped = input_file.groupby(by=[mpa])
        result = {}
        for group in PSM_grouped.groups:
            temp = PSM_grouped.get_group(group)
            sums = temp[channels].sum()
            result[group] = sums

        protein_df = pd.DataFrame.from_dict(
            result, orient='index', columns=channels)
        print("Combination done")

        return protein_df

    def protein_rollup_median(self, input_file, channels):
        '''
        This function takes Peptide level (or PSM) dataframes and performs a Median based rollup to protein level.
        the channels variable takes an array of column names that contain the quantifictions. You can create such an
        array via this command:
        channels = [col for col in PSM.columns if 'Abundance:' in col]

        mpa1 variable contains a string that is included in the Accession column. The function will search for the column containing the string
        and use it for rollup.

        Returns Protein level DF.
        '''
        mpa1=Defaults().MasterProteinAccession
        print('Calculate Protein quantifications from PSM')
        mpa = [col for col in input_file.columns if mpa1 in col]
        mpa = mpa[0]

        PSM_grouped = input_file.groupby(by=[mpa])
        result = {}
        for group in PSM_grouped.groups:
            temp = PSM_grouped.get_group(group)
            sums = temp[channels].median()
            result[group] = sums

        protein_df = pd.DataFrame.from_dict(
            result, orient='index', columns=channels)
        print("Combination done")

        return protein_df

    def protein_rollup_mean(self, input_file, channels):
        '''
        This function takes Peptide level (or PSM) dataframes and performs a Mean based rollup to protein level.
        the channels variable takes an array of column names that contain the quantifictions. You can create such an
        array via this command:
        channels = [col for col in PSM.columns if 'Abundance:' in col]

        mpa1 variable contains a string that is included in the Accession column. The function will search for the column containing the string
        and use it for rollup.

        Returns Protein level DF.
        '''
        mpa1=Defaults().MasterProteinAccession
        print('Calculate Protein quantifications from PSM')
        mpa = [col for col in input_file.columns if mpa1 in col]
        mpa = mpa[0]

        PSM_grouped = input_file.groupby(by=[mpa])
        result = {}
        for group in PSM_grouped.groups:
            temp = PSM_grouped.get_group(group)
            sums = temp[channels].mean()
            result[group] = sums

        protein_df = pd.DataFrame.from_dict(
            result, orient='index', columns=channels)
        print("Combination done")

        return protein_df

class HypothesisTesting:
    # Calculate two-sided t-test statistics for pairwise comparisons

    def __init__(self):
        self.pair_names = []
        self.comparison_data = {}

    def t_test(self, input_data, matrix1, matrix2, name=''):
        '''Calculates p values and corrected p values (q values, BH-FDR) according to a students t-test (two-sided) for each row/protein in a dataframe. Needs column names for the two matrices as arrays.
        Name defines the suffix that is added to the resulting columns e.g. comparison.
        '''
        self.pair_names.append(name)
        # Matrix 1 contains controls and Matrix 2 the treatments
        string = 'p_value'+str(name)
        string_fc = 'fold_change'+str(name)
        string_q = 'q_value' +str(name)
        for protein in input_data.index:

            m1 = input_data.loc[protein, matrix1].to_numpy()
            m2 = input_data.loc[protein, matrix2].to_numpy()
            matrix3 = np.log2(np.mean(m2)/np.mean(m1))
            t_stat, pvalue = ttest_ind(m1, m2)
            try:
                input_data.loc[protein, string] = pvalue
                input_data.loc[protein, string_fc] = matrix3
            except ValueError:
                print('Error with: '+protein)
        input_data[string] = input_data[string].fillna(
                value=1)
        pvals = input_data[string].to_numpy()

        reject, pvals_corrected, a, b = multipletests(
                pvals, method='fdr_bh')

        input_data[string_q] = pvals_corrected
        self.comparison_data = self.export_comparison_strings()
        return input_data

    def tessa(self, source):
        '''Returns all possible pairs from source as array
        '''
        result = []
        for p1 in range(len(source)):
            for p2 in range(p1+1, len(source)):
                result.append([source[p1], source[p2]])
        return result

    def peptide_based_lmm(self, input_file, conditions, norm=Preprocessing.total_intensity, pairs=None):
        
        columns=Defaults().labelsForLMM
        self.pair_names = []
        channels = [col for col in input_file.columns if columns[2] in col]
        print(channels)
        if norm is not None:
            input_file = norm(Preprocessing,input_file, channels)
        else:
            input_file = input_file.dropna(subset=channels)
            print('No Normalization applied')
        # Protein level quantifications
        roll = Rollup()
        protein_data = roll.protein_rollup_sum(
            input_file=input_file, channels=channels)
        # Prepare Peptide data for LMM
        Peptides_for_LM = input_file[channels]

        sequence = [col for col in input_file.columns if columns[0] in col]
        sequence = sequence[0]

        Peptides_for_LM['Sequence'] = input_file[sequence]

        Acc = [col for col in input_file.columns if columns[1] in col]
        Acc = Acc[0]
        Peptides_for_LM['Accession'] = input_file[Acc]

        melted_Peptides = Peptides_for_LM.melt(
            id_vars=['Accession', 'Sequence'], value_vars=channels)
        # Replace column names with conditions
        melted_Peptides.replace(to_replace=channels,
                                value=conditions, inplace=True)
        unique_conditions = list(set(conditions))
        if pairs == None:
            pairs = self.tessa(unique_conditions)
        else:
            pass
        print(pairs)
        for pair in pairs:
   
            
            pair.sort()
            print(pair)
            temp = melted_Peptides[(melted_Peptides['variable'].str.contains(pair[0], regex=False)) | (
                melted_Peptides['variable'].str.contains(pair[1], regex=False))]
            temp['value'] = np.log2(temp['value'])
            temp = temp.dropna()
           
            grouped = temp.groupby(by=['Accession'])
            result_dict = {}
            fold_changes = []
            counter = 0
            for i in grouped.groups:

                temp2 = grouped.get_group(i)
                vc = {'Sequence': '0+Sequence'}
                model = smf.mixedlm(
                    "value ~ variable + Sequence", temp2, groups='Accession', vc_formula=vc)
                try:
                    result = model.fit()
                    if counter == 0:
                        print(result.summary())
                        counter = counter + 1
                    else:
                        pass
                    fc = result.params[1]
                    pval = result.pvalues[1]
                    fold_changes.append(fc)
                    result_dict[i] = pval
                except:
                    pass

            result_df_peptides_LMM = pd.DataFrame.from_dict(
                result_dict, orient='index', columns=['p_value'])
            result_df_peptides_LMM['fold_change'] = np.array(fold_changes)
            # Multiple testing correction:
            result_df_peptides_LMM['p_value'] = result_df_peptides_LMM['p_value'].fillna(
                value=1)
            pvals = result_df_peptides_LMM['p_value'].to_numpy()

            reject, pvals_corrected, a, b = multipletests(
                pvals, method='fdr_bh')

            result_df_peptides_LMM['q_value'] = pvals_corrected

            comparison = str(pair[0]) + '_' + str(pair[1])
            self.pair_names.append(comparison)
            result_df_peptides_LMM = result_df_peptides_LMM.add_suffix(
                comparison)
            protein_data = protein_data.join(result_df_peptides_LMM)
        self.comparison_data = self.export_comparison_strings()
        return protein_data

    def get_comparisons(self):
        '''Returns all comparisons that have been performed on that dataframe for further use. 
        '''
        
        return set(self.pair_names)

    def export_comparison_strings(self):
        '''Returns a nested dictionary (json format) with all comparisons tested during hypothesis testing as keys and the column names for P values, q values and fold changes.
        '''
        data = {}
        for pair in list(set(self.pair_names)):
            data[pair] = {'name':pair,
                'pvalue':'p_value'+pair,
                'qvalue':'q_value'+pair,
                'fold_change':'fold_change'+pair
                }
        return data
    
    def get_columnnames_for_comparison(self, comparison):
        data = self.comparison_data[comparison]
        return data['fold_change'], data['pvalue'], data['qvalue']

    def get_significant_hits(self, input_data, comparison,fc_cutoff = 0.5, p_cutoff = 0.05, use_q = True):
        '''Returns all significantly regulated genes from hypothesis testing for further use in pathway enrichment analysis
        '''
        comparison_data = self.comparison_data
               
        comparison_dict = comparison_data[comparison]
        fold_change = comparison_dict['fold_change']
        if use_q == True:
            pval = comparison_dict['qvalue']
        else:
            pval = comparison_dict['pvalue']
        upregulated = input_data[(input_data[fold_change] > fc_cutoff)&(input_data[pval] < p_cutoff)]
        downregulated = input_data[(input_data[fold_change] < -fc_cutoff)&(input_data[pval] < p_cutoff)]
        genes_up = list(upregulated.index)
        genes_down = list(downregulated.index)
        data = {'up':genes_up,'down':genes_down}
        return data            

class PathwayEnrichment:
    def __init__(self):
        print("Pathway Enrichment Initialized")
        self.database = None
        self.counts = None
        self.total = 0

    def get_background_sizes(self, background=list):
        reactome_database = pd.read_csv(os.path.join(dirname, 
            "../data/UniProt2Reactome_PE_All_Levels.txt"), sep='\t', header=None)
        reactome_database.columns = ['Accession', 'Reactome_Protein_Name',
                                     'Compartment', 'Reactome_ID', 'URL', 'Description', 'Evidence_Code', 'Species']

        reactome_database['Checked'] = reactome_database['Accession'].isin(
            background)
        reactome_database = reactome_database[reactome_database['Checked'] == True]
        grouped_reactome = reactome_database.groupby(
            "Reactome_ID").agg('count')
        self.database = reactome_database
        self.counts = grouped_reactome.iloc[:, 0]
        self.total = len(background)

    def get_pathway_sizes(self, species="Homo sapiens"):
        # Read file
        reactome_database = pd.read_csv(os.path.join(dirname, 
            "../data/UniProt2Reactome_PE_All_Levels.txt"), sep='\t', header=None)
        reactome_database.columns = ['Accession', 'Reactome_Protein_Name',
                                     'Compartment', 'Reactome_ID', 'URL', 'Description', 'Evidence_Code', 'Species']
        # filter by species

        reactome_filtered = reactome_database[reactome_database['Species'].str.contains(
            species)]
        # Count proteins per pathway
        accessions = set(list(reactome_filtered['Accession']))
        grouped_reactome = reactome_filtered.groupby(
            "Reactome_ID").agg('count')
        self.database = reactome_filtered
        self.counts = grouped_reactome.iloc[:, 0]
        self.total = len(accessions)

    def get_enrichment(self, genes):
        pathways = []
        enrichmentResult = {}
        genes = list(set(genes))
        listLength = len(genes)
        database = self.database
        counts = self.counts
        for gene in genes:
            try:
                temp = database[database['Accession'].str.contains(gene)]
                listOfFoundPathways = list(temp['Reactome_ID'])
                pathways = pathways + listOfFoundPathways
            except KeyError:
                print("Gene not found in database")
        # Count pathway mentions
        pathwayMentions = Counter(pathways)
        # Calculate Hypergeometric test p value
        for pathway in list(set(pathways)):
            setSize = counts[pathway]
            foundSize = pathwayMentions[pathway]
            P_value = hypergeom.sf(
                (foundSize - 1), self.total, setSize, listLength)
            enrichmentResult[pathway] = P_value
        resultDf = pd.DataFrame.from_dict(
            enrichmentResult, orient='index', columns=['P value'])
        resultDf['Reactome_ID'] = resultDf.index
        # FDR correction
        pvals = resultDf['P value'].values.flatten()
        try:
            fdr = multipletests(pvals, method='fdr_bh')
            resultDf['FDR'] = fdr[1]
        except ZeroDivisionError:
            resultDf['FDR'] = pvals
       
        resultDf = resultDf.sort_values(by='FDR')
        resultDf = resultDf[resultDf['FDR']<0.1]
        temp_database = pd.read_csv(os.path.join(dirname, 
            "../data/ReactomePathways.txt"), sep='\t', header=None, index_col=0)
        temp_database.columns = ['Description', 'Species']
        for entry in resultDf.index:
            resultDf.loc[entry,
                         'Description'] = temp_database.loc[entry, 'Description']
        return resultDf

class Correlation:  # TODO finish coexpression clustering algo
    #PRELIMINARY DOES NOT WORK PROPERLY YET
    def __init__(self):
        pass

    def fancy_dendrogram(self, *args, **kwargs):
        max_d = kwargs.pop('max_d', None)
        if max_d and 'color_threshold' not in kwargs:
            kwargs['color_threshold'] = max_d
        annotate_above = kwargs.pop('annotate_above', 0)

        ddata = dendrogram(*args, **kwargs)

        if not kwargs.get('no_plot', False):
            plt.title('Hierarchical Clustering Dendrogram (truncated)')
            plt.xlabel('sample index or (cluster size)')
            plt.ylabel('distance')
            for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
                x = 0.5 * sum(i[1:3])
                y = d[1]
                if y > annotate_above:
                    plt.plot(x, y, 'o', c=c)
                    plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),
                                 textcoords='offset points',
                                 va='top', ha='center')
            if max_d:
                plt.axhline(y=max_d, c='k')
        return ddata

    def coexpression(self, data, columnsToUse, cutoff=0.99, conditions=[], anova=True,verbose_plotting = True):
        sns.boxplot(data=data, showfliers=False)
        # Sample array contains treatments for each column if applicable e.g. DMSO DMSO DMSO GTPP GTPP GTPP
        # Perform ANOVA between replicates to filter for high variance proteins
        data_copy = data[columnsToUse].copy()

        if anova == True:
            for protein in data.index:
                anova_df = data_copy[data_copy.index == protein]
                anova_df.columns = conditions
                anova_df_melted = pd.melt(anova_df, ignore_index=False)
                model = ols('value ~ C(variable)', data=anova_df_melted).fit()
                anova_table = sm.stats.anova_lm(model, typ=2)
                P_value = anova_table['PR(>F)'][0]
                data.loc[protein, 'P_value'] = P_value
            # Filter for changing proteins in any condition with P value 0.1 (relaxed, because only prefilter)
            data = data[data['P_value'] < 0.1]
            column_names = np.append(conditions, ['P_value'])
            data.columns = column_names
        else:
            data.columns = conditions
        data_grouped = data.groupby(by=data.columns, axis=1).mean()
        print(data_grouped.head(10))
        try:
            data_grouped = data_grouped.drop(labels=['P_value'], axis=1)
        except:
            pass
        log_data = data_grouped.transform(np.log2)
        log_data = log_data.replace([np.inf, -np.inf], np.nan).dropna(axis=0)
        trans_log_data = log_data.transpose()
        # Calcualte Correlation Matrix
        corr_matrix_pearson = trans_log_data.corr(method='pearson')
        flat_pearson = corr_matrix_pearson.values.flatten()
        if verbose_plotting == True:
            # Plot matrix
            plt.matshow(corr_matrix_pearson, cmap='RdBu')
            plt.title("Pearson")
            plt.show()
            # Plot correlation distribution
            plt.hist(flat_pearson, bins=100)
            plt.title('Pearson')
            plt.show()
        else:
            pass
        # Cluster Correlation MAtrix
        correlations_array = np.asarray(corr_matrix_pearson)
        row_linkage = linkage(correlations_array, method='average')
        plt.figure(figsize=(10, 7))
        plt.title("Dendogram of Correlation - Hierachichal clustering")

        self.fancy_dendrogram(row_linkage, truncate_mode='lastp',
                              p=200,
                              leaf_rotation=90.,
                              leaf_font_size=12.,
                              show_contracted=True,
                              annotate_above=10,
                              max_d=20)
        if verbose_plotting == True:
            plt.show()
        else:
            pass
        corr_matrix_pearson['Clusters'] = fcluster(
            row_linkage, 20, criterion='distance')
        corr_matrix_pearson['Accession'] = corr_matrix_pearson.index
        melted_corr_matrix = corr_matrix_pearson.melt(
            id_vars=['Accession', 'Clusters'], value_name='Correlation')
        filtered_melted_corr_Matrix = melted_corr_matrix[melted_corr_matrix['Correlation'] > cutoff]
        filtered_melted_corr_Matrix = filtered_melted_corr_Matrix[
            filtered_melted_corr_Matrix['Correlation'] != 1]
        Clusters = {}
        grouped_corr_matrix = corr_matrix_pearson.groupby(by=['Clusters'])
        for group in grouped_corr_matrix.groups:
            temp_df = grouped_corr_matrix.get_group(group)
            genes_in_cluster = list(temp_df.index)
            Clusters[group] = genes_in_cluster
        return filtered_melted_corr_Matrix, Clusters

class Visualization:
    #TODO add volcano plots and heatmaps for pipeline output
    def __init__(self):
        pass

    def volcano_plot(self, input_file, fold_change, pval,comparison,wd,mode='save'):
        '''Produces a volcano plot and saves it 
        '''
        temp = input_file.copy()
        
        temp[pval]=np.where(temp[pval] < 0.000000001, 0.000000001, temp[pval])#For visualization purpose
        temp['coloring']=1-np.log(temp[pval])+ abs(temp[fold_change]) * 5
        temp = temp.dropna()
        fig = sns.scatterplot(x=fold_change,y=pval, data=temp, hue='coloring',legend=False,alpha=0.6,s=12)
        fig.invert_yaxis()
        plt.axvline(x=0.5, linewidth=0.5,linestyle='dashed',color='black',alpha=0.5)
        plt.axvline(x=-0.5, linewidth=0.5,linestyle='dashed',color='black',alpha=0.5)
        plt.axhline(y=0.05,linewidth=0.5,linestyle='dashed',color='black',alpha=0.5)
        plt.yscale('log')
        plt.xlabel('Fold change (log2)')
        plt.ylabel('P value')
        plt.title(str(comparison))
        if mode == 'save':
            plt.savefig(wd+str(comparison)+'_Volcano.pdf',transparent=True)
        else:
            plt.show()
        plt.close()

    def boxplots(self,input_file, channels,wd,mode='save'):
        fig= sns.boxplot(data=input_file[channels],showfliers=False)
        plt.yscale('log')
        plt.xticks(rotation = 90)
        plt.xlabel('Sample')
        plt.ylabel('TMT intensity')
        plt.title('Sample abundances after processing')
        plt.subplots_adjust(bottom=0.35)
        if mode == 'save':
            plt.savefig(wd+str(comparison)+'_Volcano.pdf',transparent=True)
        else:
            plt.show()
        plt.close()
    
    def heatmap(self,input_file,channels,conditions,wd,mode='save'):
        temp = input_file[channels].dropna().copy()
        fig = sns.clustermap(data=temp[channels],z_score=0,xticklabels=conditions,yticklabels=False)
        if mode == 'save':
            plt.savefig(wd+str(comparison)+'_Volcano.pdf',transparent=True)
        else:
            plt.show()
        plt.close()

class Pipelines:
    def __init__(self):
        pass

    def singlefile_lmm(self, psms, conditions,pairs=None,wd=None,filter=True):
        defaults = Defaults()
        labels=defaults.labelsForLMM
        abundance_column=defaults.AbundanceColumn
        process = Preprocessing()
        hypo = HypothesisTesting()
        annot = Annotation()
        vis = Visualization()
        path = PathwayEnrichment()
        print("Initialized")
        
        channels=defaults.get_channels(psms,custom=abundance_column) #Get channel nammes
        
        if filter == True:
            print('Filtering')
            psms = process.filter_peptides(psms)
        else:
            pass
        print('Peptide based linear models for differential expression')
        
        result = hypo.peptide_based_lmm(psms,conditions=conditions,pairs=pairs)
        #Annotation
        print('Annotate')
        result = annot.basic_annotation(result)
        #vis
        print('Visualization')
        channels_02 = defaults.get_channels(result)
        vis.boxplots(result,channels_02,wd=wd)
        vis.heatmap(result,channels_02,conditions,wd=wd)
        comparisons = list(hypo.get_comparisons())
        for index in range(len(comparisons)):
            fc, p, q = hypo.get_columnnames_for_comparison(comparisons[index])
            vis.volcano_plot(result,fc,p,comparisons[index],wd=wd)
        #Pathway enrichment
        print('Pathway Enrichment')
        background = list(result.index)
        path.get_background_sizes(background)
        for index in range(len(comparisons)):
            hits = hypo.get_significant_hits(result, comparisons[index])
            up = hits['up']
            down = hits['down']
            up_pathways = path.get_enrichment(up)
            down_pathways = path.get_enrichment(down)
            up_pathways.to_csv(wd+str(comparisons[index])+'Pathways_UP.csv',line_terminator='\n')
            down_pathways.to_csv(wd+str(comparisons[index])+'Pathways_DOWN.csv',line_terminator='\n')
        print('Writing result file')
        result.to_csv(wd+"Result_groups_accession_0_Sequence.csv",line_terminator='\n')
        print('Done')
        return result

    def multifile_lmm(self, psms, conditions,bridge,pairs=None,wd=None,filter=True):
        
        defaults = Defaults()
        labels=defaults.labelsForLMM
        abundance_column=defaults.AbundanceColumn
        process = Preprocessing()
        hypo = HypothesisTesting()
        annot = Annotation()
        vis = Visualization()
        path = PathwayEnrichment()
        print("Initialized")
        channels=defaults.get_channels(psms,custom=abundance_column) #Get channel nammes
        if filter == True:
            print('Filtering')
            psms = process.filter_peptides(psms)
        else:
            pass
        print('Splitting PSMs')
        array_of_dfs=process.psm_splitting(psms)
        number_of_files = len(array_of_dfs)
        print('Number of Files:',number_of_files)
        #Normalization
        print('Normalize each file')
        array_of_dfs = defaults.processor(array_of_dfs,process.total_intensity, channels=channels)
        #join back for IRS
        print('Join for IRS')
        joined_df = process.psm_joining(array_of_dfs)
        #IRS
        print('Preparing for IRS')
        IRS_df = process.IRS_normalisation(joined_df,bridge,number_of_files)
        #LMM
        print('Peptide based linear models for differential expression')
        result = hypo.peptide_based_lmm(IRS_df,conditions=conditions,pairs=pairs)
        #Annotation
        print('Annotate')
        result = annot.basic_annotation(result)
        #vis
        print('Visualization')
        channels_02 = defaults.get_channels(result)
        vis.boxplots(result,channels_02,wd=wd)
        vis.heatmap(result,channels_02,conditions,wd=wd)
        comparisons = list(hypo.get_comparisons())
        for index in range(len(comparisons)):
            fc, p, q = hypo.get_columnnames_for_comparison(comparisons[index])
            vis.volcano_plot(result,fc,p,comparisons[index],wd=wd)
        #Pathway enrichment
        print('Pathway Enrichment')
        background = list(result.index)
        path.get_background_sizes(background)
        for index in range(len(comparisons)):
            hits = hypo.get_significant_hits(result, comparisons[index])
            up = hits['up']
            down = hits['down']
            up_pathways = path.get_enrichment(up)
            down_pathways = path.get_enrichment(down)
            up_pathways.to_csv(wd+str(comparisons[index])+'Pathways_UP.csv',line_terminator='\n')
            down_pathways.to_csv(wd+str(comparisons[index])+'Pathways_DOWN.csv',line_terminator='\n')
        print('Writing result file')
        result.to_csv(wd+"Result.csv",line_terminator='\n')
        print('Done')
        return result
    
    def wrapdynaTMT(self, psms,baseline_index=0):
        dyna = mePROD.PD_input(psms)
        dyna.IT_adjustment()
        dyna.total_intensity_normalisation()
        heavy = dyna.extract_heavy()
        peptides = dyna.baseline_correction_peptide_return(heavy,i_baseline=baseline_index)
        return peptides

    def multifile_meprod_lmm(self, psms, conditions,bridge,pairs=None,baseline_index=0,wd=None):
        defaults = Defaults()
        labels = defaults.labelsForLMM
        abundance_column=defaults.AbundanceColumn
        process = Preprocessing()
        hypo = HypothesisTesting()
        annot = Annotation()
        vis = Visualization()
        path = PathwayEnrichment()
        psms = process.filter_peptides(psms)
        print('Splitting PSMs')
        array_of_dfs=process.psm_splitting(psms)
        number_of_files = len(array_of_dfs)
        print('Number of Files:',number_of_files)
        #Normalization
        print('Normalize each file')
        array_of_dfs = defaults.processor(array_of_dfs,self.wrapdynaTMT,baseline_index=baseline_index)
        #join back for IRS
        print('Join for IRS')
        joined_df = process.psm_joining(array_of_dfs)
        #IRS
        print('Preparing for IRS')
        IRS_df = process.IRS_normalisation(joined_df,bridge,number_of_files,abundance_column=abundance_column)
        #LMM
        print('Peptide based linear models for differential expression')
        result = hypo.peptide_based_lmm(IRS_df,conditions=conditions,columns=labels,pairs=pairs,norm=None)
        #Annotation
        print('Annotate')
        result = annot.basic_annotation(result)
        #vis
        print('Visualization')
        channels_02 = defaults.get_channels(result)
        vis.boxplots(result,channels_02,wd=wd)
        vis.heatmap(result,channels_02,conditions,wd=wd)
        comparisons = list(hypo.get_comparisons())
        for index in range(len(comparisons)):
            fc, p, q = hypo.get_columnnames_for_comparison(comparisons[index])
            vis.volcano_plot(result,fc,p,comparisons[index],wd=wd)
        #Pathway enrichment
        print('Pathway Enrichment')
        background = list(result.index)
        path.get_background_sizes(background)
        for index in range(len(comparisons)):
            hits = hypo.get_significant_hits(result, comparisons[index])
            up = hits['up']
            down = hits['down']
            up_pathways = path.get_enrichment(up)
            down_pathways = path.get_enrichment(down)
            up_pathways.to_csv(wd+str(comparisons[index])+'Pathways_UP.csv',line_terminator='\n')
            down_pathways.to_csv(wd+str(comparisons[index])+'Pathways_DOWN.csv',line_terminator='\n')
        print('Writing result file')
        result.to_csv(wd+"Result.csv",line_terminator='\n')
        print('Done')
        return result

    def singlefile_meprod_lmm(self, psms, conditions,pairs=None,baseline_index=0,wd=None):
        defaults = Defaults()
        labels = defaults.labelsForLMM
        process = Preprocessing()
        hypo = HypothesisTesting()
        annot = Annotation()
        vis = Visualization()
        path = PathwayEnrichment()
        psms = process.filter_peptides(psms)
        dyna = mePROD.PD_input(psms)
        dyna.IT_adjustment()
        dyna.total_intensity_normalisation()
        heavy = dyna.extract_heavy()
        peptides = dyna.baseline_correction_peptide_return(heavy,i_baseline=baseline_index)

        result = hypo.peptide_based_lmm(peptides,conditions=conditions,columns=labels,pairs=pairs,norm=None)
        #Annotation
        print('Annotate')
        result = annot.basic_annotation(result)
        #vis
        print('Visualization')
        channels_02 = defaults.get_channels(result)
        vis.boxplots(result,channels_02,wd=wd)
        vis.heatmap(result,channels_02,conditions,wd=wd)
        comparisons = list(hypo.get_comparisons())
        for index in range(len(comparisons)):
            fc, p, q = hypo.get_columnnames_for_comparison(comparisons[index])
            vis.volcano_plot(result,fc,p,comparisons[index],wd=wd)
        #Pathway enrichment
        print('Pathway Enrichment')
        background = list(result.index)
        path.get_background_sizes(background)
        for index in range(len(comparisons)):
            hits = hypo.get_significant_hits(result, comparisons[index])
            up = hits['up']
            down = hits['down']
            up_pathways = path.get_enrichment(up)
            down_pathways = path.get_enrichment(down)
            up_pathways.to_csv(wd+str(comparisons[index])+'Pathways_UP.csv',line_terminator='\n')
            down_pathways.to_csv(wd+str(comparisons[index])+'Pathways_DOWN.csv',line_terminator='\n')
        print('Writing result file')
        result.to_csv(wd+"Result.csv",line_terminator='\n')
        print('Done')
        return result
       
        
def main():
    #Testing process for pipelines
   
    wd = 'C://Users/Kevin/Desktop/MassSpec/GroundTruth/Test/'
    psms = pd.read_csv(wd+"20210226_KKL_GroundT_F-(1)_PSMs.txt",sep='\t',header=0)
    conditions=['0C','0C','0C','1Mix1','1Mix1','1Mix','2Mix2','2Mix2','2Mix2','3Mix3','3Mix3','3Mix3']
    bridge = '129C'
    pairs = [['0C','2Mix2']]
    pipe = Pipelines()
    results = pipe.singlefile_lmm(psms,conditions,pairs=pairs,wd=wd)



if __name__ == '__main__':
   main()