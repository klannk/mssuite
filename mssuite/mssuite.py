# TODO Testing coexpression clustering
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
from scipy.stats import hypergeom, ttest_ind
from statsmodels.formula.api import ols
from statsmodels.stats.multitest import multipletests
from scipy.stats import trim_mean

warnings.filterwarnings("ignore")
import os
dirname = os.path.dirname(__file__)

class Defaults(object):
    '''
    This object contains PD specific default names for columns. Most functions will access this by default, but
    can be set manually.
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
        results=[]
        for count, value in enumerate(list_of_df):
            results.append(function(value, *args,**kwargs))
        return results

    def get_channels(self, input_file, custom=None):
        if custom is None:
            channels = [col for col in input_file.columns if self.AbundanceColumn in col]
        else:
            channels = [col for col in input_file.columns if custom in col]
        return channels

class Preprocessing:
    def __init__(self):
        pass

    def filter_peptides(self, input_file, mpa1=Defaults.MasterProteinAccession):
        '''
        Filters peptide files for non-unique peptides and contaminant proteins
        '''
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

    def IRS_normalisation(self, input_file, bridge, plexes, abundance_column=Defaults.AbundanceColumn):
        '''
        This function performs IRS normalisation for a input pandas df. Bridge channels have to be same TMT channel and plexes must have same size
        bridge = String that defines bridge channel
        plexes = Number (INT) of plexes to normalize
        quant = String that is included in all Quantification columns
        '''
        print('Internal Reference scaling')
        # search for quantification columns
        defaults=Defaults()
        channels = defaults  .get_channels(input_file=input_file, custom=abundance_column)
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

    def protein_rollup_sum(self, input_file, channels, mpa1=Defaults.MasterProteinAccession):
        '''
        This function takes Peptide level (or PSM) dataframes and performs a sum based rollup to protein level.
        the channels variable takes an array of column names that contain the quantifictions. You can create such an
        array via this command:
        channels = [col for col in PSM.columns if 'Abundance:' in col]

        mpa1 variable contains a string that is included in the Accession column. The function will search for the column containing the string
        and use it for rollup.

        Returns Protein level DF.
        '''
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

    def protein_rollup_median(self, input_file, channels,  mpa1=Defaults.MasterProteinAccession):
        '''
        This function takes Peptide level (or PSM) dataframes and performs a Median based rollup to protein level.
        the channels variable takes an array of column names that contain the quantifictions. You can create such an
        array via this command:
        channels = [col for col in PSM.columns if 'Abundance:' in col]

        mpa1 variable contains a string that is included in the Accession column. The function will search for the column containing the string
        and use it for rollup.

        Returns Protein level DF.
        '''
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

    def protein_rollup_mean(self, input_file, channels,  mpa1=Defaults.MasterProteinAccession):
        '''
        This function takes Peptide level (or PSM) dataframes and performs a Mean based rollup to protein level.
        the channels variable takes an array of column names that contain the quantifictions. You can create such an
        array via this command:
        channels = [col for col in PSM.columns if 'Abundance:' in col]

        mpa1 variable contains a string that is included in the Accession column. The function will search for the column containing the string
        and use it for rollup.

        Returns Protein level DF.
        '''
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

    def t_test(self, input_data, matrix1, matrix2, name=''):
        self.pair_names.append(name)
        # Matrix 1 contains controls and Matrix 2 the treatments
        string = 'P_value_t_test_'+str(name)
        string_fc = 'Fold_change'+str(name)
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
        return input_data

    def tessa(self, source):
        result = []
        for p1 in range(len(source)):
            for p2 in range(p1+1, len(source)):
                result.append([source[p1], source[p2]])
        return result

    def peptide_based_lmm(self, input_file, conditions, columns=Defaults.labelsForLMM, norm=Preprocessing.total_intensity, pairs=None):
        self.pair_names = []
        channels = [col for col in input_file.columns if columns[2] in col]
        if norm is not None:
            input_file = norm(Preprocessing,input_file, channels)
        else:
            input_file = input_file.dropna(subset=channels)
            print('No Normalization applied')
        # Protein level quantifications
        roll = Rollup()
        protein_data = roll.protein_rollup_sum(
            input_file=input_file, channels=channels, mpa1=columns[1])
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
                result_dict, orient='index', columns=['P value'])
            result_df_peptides_LMM['fold_change'] = np.array(fold_changes)
            # Multiple testing correction:
            result_df_peptides_LMM['P value'] = result_df_peptides_LMM['P value'].fillna(
                value=1)
            pvals = result_df_peptides_LMM['P value'].to_numpy()

            reject, pvals_corrected, a, b = multipletests(
                pvals, method='fdr_bh')

            result_df_peptides_LMM['q value'] = pvals_corrected

            comparison = str(pair[0]) + '_' + str(pair[1])
            self.pair_names.append(comparison)
            result_df_peptides_LMM = result_df_peptides_LMM.add_suffix(
                comparison)
            protein_data = protein_data.join(result_df_peptides_LMM)
        return protein_data

    def get_comparisons(self):
        '''Returns all comparisons that have been performed on that dataframe for further use. 
        '''
        return set(self.pair_names)


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
        fdr = multipletests(pvals, method='fdr_bh')
        resultDf['FDR'] = fdr[1]
        resultDf = resultDf.sort_values(by='FDR')
        temp_database = pd.read_csv(
            "../data/ReactomePathways.txt", sep='\t', header=None, index_col=0)
        temp_database.columns = ['Description', 'Species']
        for entry in resultDf.index:
            resultDf.loc[entry,
                         'Description'] = temp_database.loc[entry, 'Description']
        return resultDf


class Correlation:  # TODO finish coexpression clustering algo

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

    def coexpression(self, data, columnsToUse, cutoff=0.99, samples=[], anova=True):
        sns.boxplot(data=data, showfliers=False)
        # Sample array contains treatments for each column if applicable e.g. DMSO DMSO DMSO GTPP GTPP GTPP
        # Perform ANOVA between replicates to filter for high variance proteins
        data_copy = data.copy()

        if anova == True:
            for protein in data.index:
                anova_df = data_copy[data_copy.index == protein]
                anova_df.columns = samples
                anova_df_melted = pd.melt(anova_df, ignore_index=False)
                model = ols('value ~ C(variable)', data=anova_df_melted).fit()
                anova_table = sm.stats.anova_lm(model, typ=2)
                P_value = anova_table['PR(>F)'][0]
                data.loc[protein, 'P_value'] = P_value
            # Filter for changing proteins in any condition with P value 0.1 (relaxed, because only prefilter)
            data = data[data['P_value'] < 0.1]
            column_names = np.append(samples, ['P_value'])
            data.columns = column_names
        else:
            data.columns = samples
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

        # Plot matrix
        plt.matshow(corr_matrix_pearson, cmap='RdBu')
        plt.title("Pearson")
        plt.show()
        # Plot correlation distribution
        plt.hist(flat_pearson, bins=100)
        plt.title('Pearson')
        plt.show()
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
        plt.show()
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


class Pipelines:
    def __init__(self):
        pass

    def multifile_lmm(self, psms, conditions,bridge,pairs=None,labels=Defaults.labelsForLMM, mpa=Defaults.MasterProteinAccession,abundance_column=Defaults.AbundanceColumn):
        defaults = Defaults()
        process = Preprocessing()
        hypo = HypothesisTesting()

        channels=defaults.get_channels(psms,custom=abundance_column) #Get channel nammes
        array_of_dfs=process.psm_splitting(psms)
        number_of_files = len(array_of_dfs)
        #Normalization
        array_of_dfs = defaults.processor(array_of_dfs,process.total_intensity, channels=channels)
        #join back for IRS
        joined_df = process.psm_joining(array_of_dfs)
        #IRS
        IRS_df = process.IRS_normalisation(joined_df,bridge,number_of_files,abundance_column=abundance_column)
        #LMM
        result = hypo.peptide_based_lmm(IRS_df,conditions=conditions,columns=labels,pairs=pairs)
        return result
    
    
if __name__ == "__main__":



    '''
    
    wd="C://Users/Kevin/Desktop/MassSpec/IRS_TEST/"
    psms = pd.read_csv(wd+"20210317_KKL_Calu3_pool3_PSMs.txt",sep='\t',header=0)
    defaults = Defaults()
    process = Preprocessing()
    hypo = HypothesisTesting()
    annot= Annotation()
    psms=process.filter_peptides(psms)
    array = process.psm_splitting(psms)
    channels=defaults.get_channels(psms,custom=defaults.AbundanceColumn)
    number_of_files = len(array)
    array_of_dfs = defaults.processor(array,process.total_intensity, channels=channels)
    joined_df = process.psm_joining(array_of_dfs)
    joined_df.to_csv(wd+"Intermediate_Join.csv",line_terminator='\n')
    channels=defaults.get_channels(joined_df,custom=defaults.AbundanceColumn)

    IRS_df = process.IRS_normalisation(joined_df,'129C',number_of_files,abundance_column=defaults.AbundanceColumn)
    IRS_df.to_csv(wd+"Intermediate_IRS.csv",line_terminator='\n')

    conditions = ['0Mock','FFM1','FFM2','SA','Brasil','B117','Bridge','0Mock','FFM1','FFM2','SA','Brasil','B117','Bridge','0Mock','FFM1','FFM2','SA','Brasil','B117','Bridge']
    pairs = [['0Mock','FFM1'],['0Mock','FFM2'],['0Mock','SA'],['0Mock','Brasil'],['0Mock','B117']]
    result = hypo.peptide_based_lmm(IRS_df,conditions=conditions,pairs=pairs)
    result = annot.basic_annotation(result)
    result.to_csv(wd+"Result.csv",line_terminator='\n')
    '''
