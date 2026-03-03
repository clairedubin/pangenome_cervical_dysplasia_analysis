import lz4.frame
import pandas as pd
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sklearn.metrics
from scipy import stats
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy.optimize import curve_fit
import math
from Bio import SeqIO
import gzip
import json
from scipy.stats import mode
import networkx as nx



class GeneContent:
    def __init__(self, species_name, outdir, paths, sample_metadata, group_col='outcome', sample_list=[], 
                 species_id=None, min_frequency_threshold=0.05, max_frequency_threshold=0.95, orig_centroid_grouping=95,
                centroid_threshold=75, jaccard_threshold=0.9, min_mean_depth=5,
                min_median_marker_coverage=None,  min_unique_fraction_covered=1, min_copy_num=0.35,
                min_sample_read_count=100000,
                min_group_count=10, hide_plots=False, use_cache=True, save_cache=True, gene_profile_only=False):
        
        self.species_name = species_name
        self.outdir = outdir
        self.paths = paths  # Dictionary containing file paths
        self.species_id = species_id
        self.use_cache = use_cache
        self.save_cache = save_cache
        self.hide_plots = hide_plots
        self.min_copy_num = min_copy_num
        self.min_sample_read_count = min_sample_read_count

        if self.species_id is None:
            self.get_MIDAS_species_id()

        self.filter_suffix = f'min_copy_num{self.min_copy_num}'


        self.orig_centroid_grouping = orig_centroid_grouping
        self.centroid_threshold = centroid_threshold
        if self.centroid_threshold > self.orig_centroid_grouping:
            raise ValueError (f"Centroid grouping threshold {self.centroid_threshold} must be less than or "
                               f"equal to the original centroid grouping threshold {self.orig_centroid_grouping}")

        #load per-sample species profile and add total gene counts per sample
        self.sample_metadata = sample_metadata

        if self.min_sample_read_count:
            self.filter_suffix += f'min_read_count{self.min_sample_read_count}'
            old_sample_count = self.sample_metadata.shape[0]
            self.sample_metadata = self.sample_metadata[self.sample_metadata['num_reads'] >= self.min_sample_read_count]
            print(f'Analyzing {self.sample_metadata.shape[0]}/{old_sample_count} samples with at least {self.min_sample_read_count} reads')
        if len(sample_list) == 0:
            self.sample_list = self.sample_metadata.index.tolist()
        else:
            self.sample_list = [i for i in sample_list if i in self.sample_metadata.index.tolist()]

        self.gene_profile = self.load_gene_profile()


        self.min_mean_depth = min_mean_depth
        self.group_col = group_col
        self.min_group_count = min_group_count
        
        if not gene_profile_only:
        
            #filter by mean depth
            self.filtered_gene_profile = self.filter_by_depth()

            #calculate presence/absence based on copy number
            self.gene_presence_absence = self.calc_presence_absence()

            #group by centroid
            self.centroids = self.load_centroids()
            self.centroid_presence_absence = self.group_by_centroid()

            #filter centroid presence/absence matrix
            self.min_frequency_threshold = min_frequency_threshold
            self.max_frequency_threshold = max_frequency_threshold

            self.centroid_presence_absence_filt = self.apply_frequency_filter()

            # #calculate jaccard distances and collapse similar centroids
            self.jaccard_threshold = jaccard_threshold
            if jaccard_threshold < 1:
                self.centroid_jaccard_dists = self.calc_jaccard_dists()
                self.centroid_presence_absence_filt = self.apply_jaccard_dist_filter()

    def load_marker_gene_profile(self):

        #unique fraction covered
        file_path_ufc = f'{self.paths["MIDAS_results"]}/merge/species/species_unique_fraction_covered.tsv'
        unique_frac_cov = pd.read_csv(file_path_ufc, sep='\t', index_col=0)
        unique_frac_cov = unique_frac_cov.loc[self.species_id].to_frame(name = 'unique_fraction_covered')      

        #median marker coverage
        file_path_mmc = f'{self.paths["MIDAS_results"]}/merge/species/species_marker_median_coverage.tsv'
        marker_med_cov = pd.read_csv(file_path_mmc, sep='\t', index_col=0)
        marker_med_cov = marker_med_cov.loc[self.species_id].to_frame(name = 'median_marker_coverage')      

        df = unique_frac_cov.merge(marker_med_cov, left_index=True, right_index=True, how='outer').reset_index(names='sample')
        df = df[df['sample'].isin(self.sample_list)]
            
        samples_not_found = [i for i in  self.sample_list if i not in df['sample'].values]
        if len(samples_not_found) > 0:
            print(f'{len(samples_not_found)} samples in sample metadata file were not found in merged output files:'
                f'{samples_not_found}')

        print(f'Loaded marker gene profile for {df.shape[0]} samples')
        return df
    
    def load_genome_wide_gene_profile(self):

        to_concat = []

        for sample in self.marker_gene_profile['sample'].values:
            gs_path = f'{self.paths["MIDAS_results"]}/{sample}/genes/genes_summary.tsv'
            if not os.path.exists(gs_path):
                continue
            genes_summary = pd.read_csv(gs_path, sep='\t', index_col=0)
            if self.species_id not in genes_summary.index.tolist():
                continue

            t = genes_summary.loc[self.species_id][['covered_genes', 'fraction_covered','marker_depth', 'mean_depth']]
            to_concat += [t.rename(sample)]
        
        print(f'Loaded genome-wide depth info for {self.marker_gene_profile.shape[0]} samples')
        return pd.concat(to_concat, axis=1).T

    def calc_gene_counts_per_sample(self):
        

        num_genes_per_sample = {}
        for sample in self.marker_gene_profile['sample'].values:
            file_path = f'{self.paths["MIDAS_results"]}/{sample}/genes/{self.species_id}.genes.tsv.lz4'
            if not os.path.exists(file_path):
                continue
            try:
                temp = self.decompress_tsv_lz4_to_dataframe(file_path)
                num_genes_per_sample[sample] = temp[temp['copy_number'] >= self.min_copy_num].shape[0]
            except RuntimeError:
                continue
        
        df = pd.DataFrame.from_dict(num_genes_per_sample, orient='index')
        df.columns = ['num_genes_detected']
        print(f'Loaded gene count TSVs for {self.marker_gene_profile.shape[0]} samples')

        return df
    
    def load_gene_profile(self):


        out_path = f'{self.outdir}/{self.species_id}/gene_profile.{self.filter_suffix}.csv'

        if self.use_cache:
            try:
                df = pd.read_csv(out_path, index_col=0)
                print(f'Gene counts loaded from {out_path}')
                return df
            except FileNotFoundError:
                pass
                
        self.marker_gene_profile = self.load_marker_gene_profile()
        self.genome_wide_gene_profile = self.load_genome_wide_gene_profile()
        self.gene_counts = self.calc_gene_counts_per_sample()

        merged = self.marker_gene_profile.merge(self.genome_wide_gene_profile, how='outer', left_on='sample', right_index=True)
        merged = merged.merge(self.gene_counts, how='outer', left_on='sample', right_index=True)

        if self.save_cache:
            os.makedirs(os.path.dirname(out_path), exist_ok=True)
            merged.to_csv(out_path)

        return merged

    def filter_by_depth(self):

        df = self.gene_profile[self.gene_profile['mean_depth'] >= self.min_mean_depth]

        print(f'{df.shape[0]} samples pass filters: '
        f'mean_depth >= {self.min_mean_depth}')
        
        self.filter_suffix += f'.mean_depth{self.min_mean_depth}'
        remaining_samples = df['sample'].tolist()
        group_counts = self.sample_metadata.loc[remaining_samples][self.group_col].value_counts()
        if min(group_counts) < self.min_group_count:
            raise ValueError(f'Too few samples remaining:\n {group_counts}')
        return df
    

    @staticmethod
    def decompress_tsv_lz4_to_dataframe(compressed_file):
        decompressed_file = compressed_file.replace('.lz4', '')
        try:
            with lz4.frame.open(compressed_file, 'rb') as compressed, open(decompressed_file, 'wb') as decompressed:
                decompressed.write(compressed.read())
            df = pd.read_csv(decompressed_file, sep='\t')
            os.remove(decompressed_file)
            return df
        except Exception as e:
            print(f"Error loading {compressed_file}: {e}")
            if os.path.exists(decompressed_file):
                os.remove(decompressed_file)
            


    def plot_marker_coverage_vs_num_genes(self, vert_lines=[5, 10]):

        plt.scatter(self.gene_profile['median_marker_coverage'], 
                    self.gene_profile['num_genes_detected'], s=3)
        plt.xlabel('Median Marker Coverage')
        plt.ylabel('Num Genes Detected')
        for v in vert_lines:
            plt.axvline(math.ceil(v), c='red', label=f'x={math.ceil(v)}', linestyle='--')
        
        plt.title(self.species_name)
        plt.legend()
        plt.show()

    def calc_marker_threshold(self):
        x = self.gene_profile['median_marker_coverage']
        y = self.gene_profile['num_genes_detected']

        # Find the piecewise regression-based inflection point
        x_piecewise_inflection = find_piecewise_inflection(x, y)

        # Plot results
        if ~self.hide_plots:
            self.plot_marker_coverage_vs_num_genes(vert_lines=[x_piecewise_inflection])

        return math.ceil(x_piecewise_inflection)

    def filter_by_marker_coverage(self):

        high_cov = self.gene_profile[(self.gene_profile['unique_fraction_covered'] >= self.min_unique_fraction_covered)\
                                        & (self.gene_profile['median_marker_coverage'] >= self.min_median_marker_coverage)]
        
        print(f'{high_cov.shape[0]} samples pass marker coverage filters: '
                f'unique_fraction_covered >= {self.min_unique_fraction_covered} and '
                f'median_marker_coverage >= {self.min_median_marker_coverage}')
        
        self.filter_suffix += f'.unique_fraction_covered{self.min_unique_fraction_covered}.median_marker_coverage{self.min_median_marker_coverage}'
        high_cov_samples = high_cov['sample'].tolist()
        group_counts = self.sample_metadata.loc[high_cov_samples][self.group_col].value_counts()
        if min(group_counts) < self.min_group_count:
            raise ValueError(f'Too few samples remaining:\n {group_counts}')
        return high_cov
    
    def calc_presence_absence(self):

        file_path = f'{self.paths["MIDAS_results"]}/merge/genes/{self.species_id}/{self.species_id}.genes_copynum.tsv.lz4'
        copy_num = self.decompress_tsv_lz4_to_dataframe(file_path)        
        copy_num = copy_num.set_index(f'cluster_{self.orig_centroid_grouping}_id')
        pres_abs = copy_num.map(lambda x: 0 if x < self.min_copy_num else 1)
        pres_abs.columns = [i for i in pres_abs.columns]
        pres_abs = pres_abs[[s for s in self.filtered_gene_profile['sample'].unique() if s in pres_abs.columns]].T
        return pres_abs

    def load_centroids(self):

        centroid_path = f'{self.paths["pangenomes"]}/{self.species_id}/clusters_99_info.tsv'

        if os.path.exists(centroid_path):
            
            centroids = pd.read_csv(centroid_path, sep='\t')
            return centroids



    def group_by_centroid(self):

        if self.centroid_threshold == self.orig_centroid_grouping:
            gene_detected_c = (self.gene_presence_absence > 0).astype(int)

        else:
        # if True:
            grouped_centroids = self.centroids.groupby(f'centroid_{self.centroid_threshold}')[f'centroid_{self.orig_centroid_grouping}'].apply(list).to_dict()
            to_concat  = []
            for c, c_orig_list in grouped_centroids.items():

                cols = [i for i in c_orig_list if i in self.gene_presence_absence.columns]

                new_col = self.gene_presence_absence[cols].sum(axis=1)
                new_col.name = c
                to_concat += [new_col.to_frame().T]
                
            gene_detected_c = pd.concat(to_concat)
            gene_detected_c = (gene_detected_c > 0).astype(int).T

            print(f'{self.gene_presence_absence.shape[1]} C{self.orig_centroid_grouping} centroids have been grouped into ' \
                f'{gene_detected_c.shape[1]} C{self.centroid_threshold} centroids')

        gene_detected_c = gene_detected_c[gene_detected_c.sum(axis=1) > 0]
        print(f'{gene_detected_c.shape[1]} C{self.centroid_threshold} centroids' \
         f' are detected in at least one sample')
        self.filter_suffix += f'.C{self.centroid_threshold}'

        return gene_detected_c


    def apply_frequency_filter(self):
        freqs = self.centroid_presence_absence.sum(axis=0) / self.centroid_presence_absence.shape[0]
        freq_filt = self.centroid_presence_absence.T.loc[freqs[(freqs >= self.min_frequency_threshold) & (freqs < self.max_frequency_threshold)].index.tolist()]

        print(f'{freq_filt.shape[0]} of {self.centroid_presence_absence.shape[1]} centroids remain '
              f'after applying {self.min_frequency_threshold*100}/{(self.max_frequency_threshold)*100}% min/max frequency filter')
        self.filter_suffix += f'.frequency{self.min_frequency_threshold}-{self.max_frequency_threshold}'
        return freq_filt.T

    def calc_jaccard_dists(self):

        outfile = f'jaccard_dists.{self.filter_suffix}.csv.gz'
        out_path = f'{self.outdir}/{self.species_id}/{outfile}'

        if self.use_cache:
            try:
                df = pd.read_csv(out_path, index_col=0, compression='gzip')
                print(f'Jaccard distances loaded from {out_path}')
                return df
            except FileNotFoundError:
                pass
        
        print(f'Calculating Jaccard distances for {self.centroid_presence_absence_filt.shape[1]} centroids')
        #transpose to iterate through rows
        centroid_presence_absence_filt_T = self.centroid_presence_absence_filt.T
        jaccards = np.zeros((centroid_presence_absence_filt_T.shape[0], centroid_presence_absence_filt_T.shape[0]))
        for i, (_, row1) in enumerate(centroid_presence_absence_filt_T.iterrows()):
            for j, (_, row2) in enumerate(centroid_presence_absence_filt_T.iloc[i:].iterrows()):
                j_idx = i + j
                jaccard = sklearn.metrics.jaccard_score(row1, row2)
                jaccards[i, j_idx] = jaccard
                jaccards[j_idx, i] = jaccard
        
        jaccard_df = pd.DataFrame(jaccards, columns=centroid_presence_absence_filt_T.index, index=centroid_presence_absence_filt_T.index)
        
        if self.save_cache:
            os.makedirs(os.path.dirname(out_path), exist_ok=True)
            jaccard_df.to_csv(out_path, compression='gzip')
        
        return jaccard_df

    def apply_jaccard_dist_filter(self):

        jaccard_dists = self.centroid_jaccard_dists
        
        # Build an undirected graph
        G = nx.Graph()

        # Add all genes as nodes
        G.add_nodes_from(jaccard_dists.index)

        for i, gene1 in enumerate(jaccard_dists.index):
            for j in range(i + 1, len(jaccard_dists.columns)):
                gene2 = jaccard_dists.columns[j]
                if jaccard_dists.iat[i, j] > self.jaccard_threshold:
                    G.add_edge(gene1, gene2)

        # Find connected components (each one is a group)
        groups = list(nx.connected_components(G))

        # Map each gene to its group number
        gene_to_group = {}
        for group_id, group in enumerate(groups):
            for gene in group:
                if self.jaccard_threshold < 1:
                    gene_to_group[gene] = f'CAG_{group_id}'
                else:
                    gene_to_group[gene] = gene

        group_df = pd.DataFrame.from_dict(gene_to_group, orient='index', columns=['group']).reset_index()
        group_df = group_df.rename(columns={'index': 'gene'})

        self.filter_suffix += f'.jaccardmax{self.jaccard_threshold}'
        if self.save_cache:
            outfile = f'{self.outdir}/{self.species_id}/centroid_groupings.{self.filter_suffix}.csv'
            group_df.to_csv(outfile, index=False)

        
        centroid_to_group = group_df.set_index('gene')['group'].to_dict()
        df = self.centroid_presence_absence_filt.copy()
        df_grouped = df.rename(columns=centroid_to_group)


        collapsed_df = (
            self.centroid_presence_absence_filt
            .rename(columns=centroid_to_group)
            .T
            .groupby(level=0)
            .agg(lambda x: x.mode().iloc[0] if not x.mode().empty else pd.NA)
            .T
        )


        print(f'Collapsed {jaccard_dists.shape[0]} centroids into '
              f'{collapsed_df.shape[1]} co-abundant groups with Jaccard similarity >= {self.jaccard_threshold}')



        return collapsed_df



    def get_gene_groups(self, genes):
        
        gene_group_dict = {}

        if isinstance(genes, str):
            genes = [genes]

        for g in genes:
            gene_group_dict[g] = self.jaccard_collapsed_genes[g]

        return gene_group_dict
    
