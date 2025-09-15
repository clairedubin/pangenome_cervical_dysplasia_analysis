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
        # gene_info_path = f'{self.paths["pangenomes"]}/{self.species_id}/genes_info.tsv'

        if os.path.exists(centroid_path):
            
            centroids = pd.read_csv(centroid_path, sep='\t')
            return centroids

        # elif os.path.exists(gene_info_path):
        #     centroids = pd.read_csv(gene_info_path, sep='\t')
        #     return centroids

        # return centroids

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
                gene_to_group[gene] = f'CAG_{group_id}'

        group_df = pd.DataFrame.from_dict(gene_to_group, orient='index', columns=['group']).reset_index()
        group_df = group_df.rename(columns={'index': 'gene'})

        self.filter_suffix += f'.jaccardmax{self.jaccard_threshold}'
        if self.save_cache:
            outfile = f'{self.outdir}/{self.species_id}/centroid_groupings.{self.filter_suffix}.csv'
            group_df.to_csv(outfile, index=False)

        
        centroid_to_group = group_df.set_index('gene')['group'].to_dict()
        df = self.centroid_presence_absence_filt.copy()
        df_grouped = df.rename(columns=centroid_to_group)

        # Step 3: Group columns by group number and take the median across each group
        # Transpose, group by column group (originally centroids), take median, then transpose back
        # collapsed_df = (
        #     self.centroid_presence_absence_filt
        #     .rename(columns=centroid_to_group)
        #     .T
        #     .groupby(level=0)
        #     .median()
        #     .T
        # )

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

        # all_genes = self.centroid_jaccard_dists.index.tolist()

        # similar_genes = {}
        # dropped_genes = []

        # for i, gene1 in enumerate(all_genes):
        #     if gene1 in dropped_genes:
        #         continue
            
        #     for j, gene2 in enumerate(all_genes[i:]):
        #         if gene1 == gene2:
        #             continue
                
        #         if self.centroid_jaccard_dists.loc[gene1, gene2] > self.jaccard_threshold:
        #             if gene1 in similar_genes:
        #                 similar_genes[gene1] += [gene2]
        #             else:
        #                 similar_genes[gene1] = [gene2]
                        
        #             dropped_genes += [gene2]

        # genes_to_keep = [i for i in all_genes if i not in dropped_genes]
        # self.jaccard_collapsed_genes = similar_genes
        
        # print(f'Collapsing {self.centroid_jaccard_dists.shape[0]} centroids into '
        #       f'{len(genes_to_keep)} co-abundant groups with Jaccard distance < {self.jaccard_threshold}')
        
        # self.filter_suffix += f'.jaccardmax{self.jaccard_threshold}'

        # if self.save_cache:
        #     outfile = f'{self.outdir}/{self.species_id}/coabundant_centroids.{self.filter_suffix}.json'
        #     with open(outfile, 'w') as f:
        #         json.dump(similar_genes, f)

        # return self.centroid_presence_absence_filt[genes_to_keep]
    
    def run_GEE(self, formula, sample_metadata_df, group_by='subject_id', \
            cov_structure='independence',  \
            family=sm.genmod.families.Binomial(link=sm.genmod.families.links.Logit())):
        

        if cov_structure.lower() == 'independence':
            cov_struct = sm.cov_struct.Independence()

        elif cov_structure.lower() == 'exchangeable':
            cov_struct = sm.cov_struct.Exchangeable()


        outfile = f'GEE_results.{self.filter_suffix}.{cov_structure.lower()}.csv'
        out_path = f'{self.outdir}/{self.species_id}/{outfile}' 
        
        if self.save_cache:
            try:
                res = pd.read_csv(out_path, index_col=0)
                print(f'GEE results loaded from {out_path}')
                sig = res[res['BH_corr_p']<0.05]
                print(f'{sig.shape[0]} gene groups with corrected p<0.05')
                return res
            except FileNotFoundError:
                pass
        
        results = {}
        
        data = self.centroid_presence_absence_filt.merge(sample_metadata_df, left_index=True, right_index=True, how='left')
        groups = data[group_by].to_numpy()
        metadata_cols = sample_metadata_df.columns.tolist()
        print(f'Running GEE on {self.centroid_presence_absence_filt.shape[1]} gene groups')
        for gene in self.centroid_presence_absence_filt.columns:

            #temporarily replace column name for compatibility with GEE
            temp = data[[gene]+metadata_cols]
            temp_gene_name = gene.replace('.','_')
            if temp_gene_name != gene:
                temp = temp.rename(columns={gene:temp_gene_name})

            gene_formula = formula.replace('gene', temp_gene_name)
            model = sm.GEE.from_formula(gene_formula, 
                                        groups = groups,
                                        data = temp,
                                        cov_struct = cov_struct,
                                        family = family)
            result = model.fit()

            #counts of samples with genes pres/abs by outcome
            t = data[['outcome', gene]].value_counts()
            sample_counts = []
            for outcome  in ['full term', 'preterm']:
                for pres in [0,1]:

                    try: count = t.loc[outcome][pres]
                    except KeyError: count = 0

                    sample_counts += [count]

            results[gene] = [result.params[temp_gene_name], result.pvalues[temp_gene_name]]  + sample_counts

        res = pd.DataFrame.from_dict(results, orient='index', columns=['coef', 'p', 'ftb_absent', 'ftb_present', 'ptb_absent', 'ptb_present'])
        res['exp(coef)'] = np.exp(res['coef'])
        # res['log2(exp(coef))'] = np.log2(res['exp(coef)'])

        # print(res.head())
        res = res[~res['p'].isna()]
        if res.shape[0] == 0:
            print('results are NA for all genes')
            print(res.head())
            return pd.DataFrame()
        res['BH_corr_p'] = multipletests(res['p'], method='fdr_bh')[1]
        res['-log10(BH_corr_p)'] = -1*np.log10(res['BH_corr_p'])
        res['-log10(p)'] = -1*np.log10(res['p'])
        res = res.sort_values('p')

        #add info about other centroids in group
        res['num_centroids_in_group'] = res.index.map({k:1+len(v) for k,v in self.jaccard_collapsed_genes.items()})
        res['centroids_in_group'] = res.index.map({k:';'.join([k]+v) for k,v in self.jaccard_collapsed_genes.items()})

        #fill values for groups with only one member
        res['num_centroids_in_group'] = res['num_centroids_in_group'].fillna(1)
        res['centroids_in_group'] = res['centroids_in_group'].fillna(res.index.to_series())

        if self.save_cache:
            os.makedirs(os.path.dirname(out_path), exist_ok=True)
            res.to_csv(out_path)
            print(f'GEE results saved to {out_path}')

        self.GEE_result = res
        sig = res[res['BH_corr_p']<0.05]
        
        print(f'{sig.shape[0]} gene groups with corrected p<0.05')
        return res
    
    def run_LME(self, formula, sample_metadata_df, group_by='subject_id', \
            cov_structure='independence',  \
            family=sm.genmod.families.Binomial(link=sm.genmod.families.links.Logit())):
        
        #i dont know why but it errors if this isn't imported within the function
        import statsmodels.formula.api as smf



        outfile = f'LME_results.{self.filter_suffix}.csv'
        out_path = f'{self.outdir}/{self.species_id}/{outfile}' 
        
        # if self.save_cache:
        #     try:
        #         res = pd.read_csv(out_path, index_col=0)
        #         print(f'LME results loaded from {out_path}')
        #         sig = res[res['BH_corr_p']<0.05]
        #         print(f'{sig.shape[0]} gene groups with corrected p<0.05')
        #         return res
        #     except FileNotFoundError:
        #         pass
        
        results = {}
        
        data = self.centroid_presence_absence_filt.merge(sample_metadata_df, left_index=True, right_index=True, how='left')
        groups = data[group_by].to_numpy()
        metadata_cols = sample_metadata_df.columns.tolist()
        print(f'Running LME on {self.centroid_presence_absence_filt.shape[1]} gene groups')
        for gene in self.centroid_presence_absence_filt.columns:

            #temporarily replace column name for compatibility with GEE
            temp = data[[gene]+metadata_cols]
            temp_gene_name = gene.replace('.','_')
            if temp_gene_name != gene:
                temp = temp.rename(columns={gene:temp_gene_name})

            gene_formula = formula.replace('gene', temp_gene_name)
            model = smf.mixedlm(gene_formula, 
                                        groups = groups,
                                        data = temp,
                                        )
            result = model.fit()

            #counts of samples with genes pres/abs by outcome
            t = data[['outcome', gene]].value_counts()
            sample_counts = []
            for outcome  in ['full term', 'preterm']:
                for pres in [0,1]:

                    try: count = t.loc[outcome][pres]
                    except KeyError: count = 0

                    sample_counts += [count]

            results[gene] = [result.params[temp_gene_name], result.pvalues[temp_gene_name]]  + sample_counts

        res = pd.DataFrame.from_dict(results, orient='index', columns=['coef', 'p', 'ftb_absent', 'ftb_present', 'ptb_absent', 'ptb_present'])
        res['exp(coef)'] = np.exp(res['coef'])

        res = res[~res['p'].isna()]
        if res.shape[0] == 0:
            print('results are NA for all genes')
            print(res.head())
            return pd.DataFrame()
        res['BH_corr_p'] = multipletests(res['p'], method='fdr_bh')[1]
        res['-log10(BH_corr_p)'] = -1*np.log10(res['BH_corr_p'])
        res['-log10(p)'] = -1*np.log10(res['p'])
        res = res.sort_values('p')

        #add info about other centroids in group
        res['num_centroids_in_group'] = res.index.map({k:1+len(v) for k,v in self.jaccard_collapsed_genes.items()})
        res['centroids_in_group'] = res.index.map({k:';'.join([k]+v) for k,v in self.jaccard_collapsed_genes.items()})

        #fill values for groups with only one member
        res['num_centroids_in_group'] = res['num_centroids_in_group'].fillna(1)
        res['centroids_in_group'] = res['centroids_in_group'].fillna(res.index.to_series())

        if self.save_cache:
            os.makedirs(os.path.dirname(out_path), exist_ok=True)
            res.to_csv(out_path)
            print(f'GEE results saved to {out_path}')

        self.GEE_result = res
        sig = res[res['BH_corr_p']<0.05]
        
        print(f'{sig.shape[0]} gene groups with corrected p<0.05')
        return res
    
    def get_gene_sequences(self, gene_list, print_seqs=True, as_AA=False, save=None):

        to_write = ''
        if isinstance(gene_list, str):
            gene_list = [gene_list]

        centroids_path = f'{self.paths["pangenomes"]}/{self.species_id}/centroids.ffn'
        if os.path.exists(centroids_path):
            record_dict = SeqIO.index(centroids_path, "fasta")

        elif os.path.exists(centroids_path+'.gz'):
            record_dict = {}
            with gzip.open(centroids_path+'.gz', "rt") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    record_dict[record.id] = record
        else:
            raise ValueError('Centroid sequence file not found')

        seq_dict = {}
        for g in gene_list:

            seq = record_dict[g].seq
            if as_AA:
                seq = str(seq.translate())
            else:
                seq = str(seq)
                seq_dict[g] = str(seq)
            seq_dict[g] = seq
            if print_seqs:
                print(f'>{g}')
                print(seq)
                print()

            if save:
                to_write += f'>{g}\n{seq}\n'

        if save:
            with open(save, 'w') as f:
                f.write(to_write)

        return seq_dict

    def get_gene_groups(self, genes):
        
        gene_group_dict = {}

        if isinstance(genes, str):
            genes = [genes]

        for g in genes:
            gene_group_dict[g] = self.jaccard_collapsed_genes[g]

        return gene_group_dict
    
def plot_GEE(GEE_result, sig_threshold=0.05, coef_max=None, title=None):

    if coef_max:
        GEE_result = GEE_result[abs(GEE_result['coef']) < coef_max]

    not_sig = GEE_result[GEE_result['BH_corr_p']>=sig_threshold]
    sig = GEE_result[GEE_result['BH_corr_p']<sig_threshold]
    
    fig, ax = plt.subplots(1,1, figsize=(7,5))

    sns.scatterplot(x=sig['exp(coef)'], y=sig['-log10(p)'], s=80, color='darkslateblue',
                    # hue=sig['coef']>0, 
            # palette=['lightskyblue','darkslateblue'], 
            legend=False,  linewidth=0, marker='*', ax=ax)
    sns.scatterplot(x=not_sig['exp(coef)'], y=not_sig['-log10(p)'], s=20, color='darkgrey', 
            legend=False, linewidth=0, alpha=0.7, ax=ax)
    
    if title:
            ax.set_title(title)

    ax.set_xlabel('Odds Ratio')
    ax.set_ylabel(r'-log$_{10}$(p)')
    plt.show()

def plot_frequency(GEE_result, gene_list, annotate=True):

    if isinstance(gene_list, str):
        gene_list = [gene_list]

    num_rows, num_cols = 1+(len(gene_list)//5), min(len(gene_list),5)
    # Setting up subplots for each gene
    fig, axes = plt.subplots(num_rows, num_cols,
                            figsize=(3*num_cols, 4*num_rows), sharey=True)
    if len(gene_list) == 1:
        axes = np.array(axes)
    for ax, gene in zip(axes.reshape(-1), gene_list):
        
        orig = GEE_result.loc[gene]
        gene_res = orig.copy()

        ftb_total = gene_res[['ftb_present','ftb_absent']].sum()
        ptb_total = gene_res[['ptb_present','ptb_absent']].sum()

        gene_res[['ftb_present','ftb_absent']] = gene_res[['ftb_present','ftb_absent']]/ftb_total
        gene_res[['ptb_present','ptb_absent']] = gene_res[['ptb_present','ptb_absent']]/ptb_total

        ax.bar([0,1], gene_res[['ftb_present', 'ptb_present']], color='purple', label='gene present')
        ax.bar([0,1], gene_res[['ftb_absent', 'ptb_absent']], bottom=gene_res[['ftb_present', 'ptb_present']],
            color='thistle', label='gene absent')
        

        # # Create the stacked bars
        # bar1 = ax.bar([0, 1], gene_res[['ftb_present', 'ptb_present']], color='purple', label='gene present')
        # bar2 = ax.bar([0, 1], gene_res[['ftb_absent', 'ptb_absent']], bottom=gene_res[['ftb_present', 'ptb_present']],
        #               color='thistle', label='gene absent')

        # Adding text on top of bars

        if annotate:
            for i in range(2):  # Loop through Full Term (0) and Preterm (1)

                present_height = gene_res['ftb_present' if i == 0 else 'ptb_present']
                absent_height = gene_res['ftb_absent' if i == 0 else 'ptb_absent']

                if present_height == 0:
                    present_height += 0.05
                if absent_height == 0:
                    absent_height -= 0.05

                present_text = orig['ftb_present' if i == 0 else 'ptb_present']
                absent_text = orig['ftb_absent' if i == 0 else 'ptb_absent']

                # Adding text for 'gene present'
                ax.text(i, (present_height / 2), int(present_text), 
                        ha='center', va='center', color='white', fontsize=12,)

                # Adding text for 'gene absent'
                ax.text(i, (present_height + absent_height / 2), int(absent_text),
                        ha='center', va='center', color='black', fontsize=12, )

        
        ax.set_title(f'{gene}\np={round(gene_res["BH_corr_p"],4)}, coef={round(gene_res["coef"],4)}')
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['Full Term', 'Preterm'])
        # ax.set_xlabel('Outcome')
        ax.spines[['right', 'top']].set_visible(False)
        
    axes.reshape(-1)[len(gene_list)-1].legend()
    #delete empty plots
    for i in range(len(gene_list),len(axes.reshape(-1))):
        fig.delaxes(axes.reshape(-1)[i])
    plt.tight_layout()
    plt.show()

# def piecewise_linear(x, x0, y0, k1, k2):
#         """
#         Defines a piecewise linear function with two segments.
        
#         Parameters:
#             x  : Independent variable
#             x0 : Breakpoint (inflection point)
#             y0 : Value at the breakpoint
#             k1 : Slope before the breakpoint
#             k2 : Slope after the breakpoint

#         Returns:
#             Piecewise linear function values.
#         """
#         return np.piecewise(x, [x < x0], 
#                             [lambda x: k1 * (x - x0) + y0, 
#                             lambda x: k2 * (x - x0) + y0])

# def find_piecewise_inflection(x, y):
#         """
#         Finds the inflection point using piecewise linear regression.
        
#         Parameters:
#             x (array-like): X-axis data (independent variable)
#             y (array-like): Y-axis data (dependent variable)

#         Returns:
#             float: X value where the inflection point occurs
#         """
#         # Sort x and y to ensure they are in ascending order
#         sorted_indices = np.argsort(x)
#         x = np.array(x)[sorted_indices]
#         y = np.array(y)[sorted_indices]

#         # Initial guesses for fitting (median x as breakpoint, median y, slopes)
#         p0 = [np.median(x), np.median(y), -1, 0]

#         # Fit the piecewise linear model
#         params, _ = curve_fit(piecewise_linear, x, y, p0=p0)

#         # Extract the inflection point (breakpoint)
#         x_inflection = params[0]
        
#         return x_inflection


def piecewise_linear(x, x0, y0, k1, k2):
    return np.piecewise(
        x, [x < x0],
        [lambda x: k1 * (x - x0) + y0,
         lambda x: k2 * (x - x0) + y0]
    )

def find_piecewise_inflection(x, y, x_min=5):
    """
    More robust inflection finder that ignores low-x noise and encourages breakpoint in the flatter region.
    """
    # Filter out very low x values (noise region)
    x = np.array(x)
    y = np.array(y)
    mask = x >= x_min
    x = x[mask]
    y = y[mask]

    # Sort by x
    sorted_idx = np.argsort(x)
    x = x[sorted_idx]
    y = y[sorted_idx]

    # Set initial guess
    x0_initial = np.percentile(x, 50)
    y0_initial = np.percentile(y, 50)
    k1_initial = 10
    k2_initial = 0.5

    p0 = [x0_initial, y0_initial, k1_initial, k2_initial]

    # Set bounds to avoid low-x breakpoints
    lower_bounds = [np.percentile(x, 25), min(y), -np.inf, -np.inf]
    upper_bounds = [np.percentile(x, 90), max(y), np.inf, np.inf]

    params, _ = curve_fit(
        piecewise_linear, x, y, p0=p0,
        bounds=(lower_bounds, upper_bounds)
    )

    return params[0]  # x0 = inflection point
