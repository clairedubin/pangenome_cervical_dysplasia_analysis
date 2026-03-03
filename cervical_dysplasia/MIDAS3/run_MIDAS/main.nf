#!/usr/bin/env nextflow

nextflow.enable.dsl=2

//// PARAMS
params.merge_genes_cluster_level = 95
params.results_dir = workflow.launchDir + "/results_C${params.merge_genes_cluster_level}"

params.bowtie_db_type = 'per_sample'

// Ensure directories end with trailing "/" characters
params.findAll { key, _ -> key.endsWith("_dir")}
    .each { key, path ->
    if (!path.endsWith("/")) {
        params[key.replace("_dir", "_path")] = "${path}/"
    } else {
        params[key.replace("_dir", "_path")] = "${path}"
    }
}

workflow {

    samples = Channel
        .fromPath(params.fastq_csv)
        .splitCsv(header: false)

    RunSpecies(samples)
    MergeSpecies(RunSpecies.out.sample_id.collect())

    if( params.db_name == 'localdb' ) {
        species_to_prune = MergeSpecies.out.species_list.splitCsv(header: true, sep: '\t')
                                                    .map { it.values() }
                                                    .flatten()
    }
    else {
        DownloadDB(MergeSpecies.out.species_list)
        species_to_prune = DownloadDB.out.splitCsv(header: true, sep: '\t')
                                                    .map { it.values() }
                                                    .flatten()
    }
    
    // pruned_species = PruneCentroids(species_to_prune).collect()
    // BuildBowtieDB(pruned_species)

    // pruned_species.view()
    // PruneCentroids(species_to_prune).collect()
    // RunGenes(PruneCentroids.out,
    //         samples.join(RunSpecies.out.sample_profile_tuple)
    //         // .combine(BuildBowtieDB.out.pangenomes_species)
    //     )

    pruned_species = PruneCentroids(species_to_prune).species.collect()
    pruned_done = PruneCentroids.out.pruned_done.collect()

    if( params.bowtie_db_type == 'per_species' ) {
        
        BuildBowtieDB(pruned_species)
        RunGenes_PreBuiltBowtieDB(
            samples.join(RunSpecies.out.sample_profile_tuple)
            .combine(BuildBowtieDB.out.pangenomes_species)
            )
        run_genes_out = RunGenes_PreBuiltBowtieDB.out.sample_id.collect()
        
    }
    else {
        RunGenes(
            pruned_done,
            samples.join(RunSpecies.out.sample_profile_tuple)
        )
        run_genes_out = RunGenes.out.sample_id.collect()

    }

     MergeGenes(run_genes_out)


}

process RunSpecies {

    label 'mem_medium'
    publishDir "${params.results_path}/", mode: "copy"

    input:
    tuple val(sample_id), path(R1), path(R2)

    output:
    val(sample_id), emit: sample_id
    tuple val(sample_id), path("${sample_id}/species/species_profile.tsv"), emit: sample_profile_tuple
    path("${sample_id}/species/log.txt")

    script:
    """
    #! /usr/bin/env bash
    set -e
    set -x
    midas run_species \
        --sample_name ${sample_id} \
        -1 ${R1} -2 ${R2} \
        --midasdb_name ${params.db_name} \
        --midasdb_dir ${params.db_path} \
        --num_cores ${task.cpus} \
        .
    """
}

process MergeSpecies {

    label 'single_cpu'
    publishDir "${params.results_path}/merge/", mode: "copy"

    input:
    val(sample_id_list)

    output:
    path('species/species_prevalence.tsv'), emit: species_prevalence
    path('species_list.tsv'), emit: species_list
    path('species/*.tsv')
    path('species/species_log.txt')

    shell:
    """
    #! /usr/bin/env bash

    echo -e "sample_name\tmidas_outdir" > sample_list.tsv
    echo !{sample_id_list} | sed 's/[][]//g' | tr ',' '\n' | tr -d ' ' | \
    while read sample; do echo -e "\$sample\t!{params.results_path}"; done >> sample_list.tsv

    midas merge_species \
        --samples_list sample_list.tsv \
        --min_cov 0.75 \
        .

    echo "species" >> species_list.tsv
    awk '\$6 >= 1 {print \$1}' species/species_prevalence.tsv | awk '(NR>1)' >> species_list.tsv
    """
}

process DownloadDB {
    
    label 'local'
    
    input:
    path(species_list)

    output:
    path("species_list.tsv")

    script:
    """
    #! /usr/bin/env bash

    tail -n +2 ${species_list} > species_list_no_header.tsv

    
    echo "downloading"

    midas database \
        --download \
        --midasdb_name ${params.db_name} \
        --midasdb_dir ${params.db_path} \
        --species_list species_list_no_header.tsv \
        --num_cores ${task.cpus}
    """

}

process PruneCentroids {

    label 'single_cpu'
    publishDir "${params.db_path}/pangenomes/${species}/pruned", mode: "copy"
    
    input:
    val(species)

    output:
    val(species), emit: 'species'
    val(true), emit: 'pruned_done'
    
    """
    midas prune_centroids \
        --midasdb_name ${params.db_name} \
        --midasdb_dir ${params.db_path} \
        --species ${species} \
        --remove_singleton \
        -t ${task.cpus} \
        --force \
        --debug
    """
}

process BuildBowtieDB {

    label 'mem_medium'
    publishDir "${params.db_path}/bt2_indexes", mode: "copy"

    input:
    val(species_list)

    output:
    path('pangenomes.species'), emit: pangenomes_species
    path('pangenomes*.bt2*') //will end in .bt2 for small index and .bt2l for large index
    path('pangenomes.fa')
    path('bt2-db-build-pangenomes.log')
    
    """
    midas build_bowtie2db \
        --midasdb_name ${params.db_name} \
        --midasdb_dir ${params.db_path} \
        --species_list ${species_list.join(',' )} \
        --bt2_indexes_name pangenomes \
        --bt2_indexes_dir . \
        --num_cores ${task.cpus}
    """
}

process RunGenes {
    
    label 'mem_high'
    publishDir "${params.results_path}/", mode: "copy"
    maxRetries 2
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }

    input:
    val(pruned_done) //placeholder to make sure pruning is done before this RunGenes starts
    tuple val(sample_id), \
    path(R1), path(R2), \
    path(species_profile)
    // \
    // path(pangenomes_species)

    output:
    val("${sample_id}"), emit: sample_id
    path("${sample_id}/genes/log.txt")
    path("${sample_id}/genes/*.genes.tsv.lz4"), optional: true
    path("${sample_id}/genes/genes_summary.tsv"), optional: true
    
    """
    #! /usr/bin/env bash

    mkdir -p ${sample_id}/species/
    cp species_profile.tsv ${sample_id}/species/

    module load CBI samtools

    {
        midas run_genes \
            --sample_name ${sample_id} \
            -1 ${R1} -2 ${R2} \
            --midasdb_name ${params.db_name} \
            --midasdb_dir ${params.db_path} \
            --num_cores ${task.cpus} \
            --select_by median_marker_coverage,unique_fraction_covered \
            --select_threshold 2,0.5 \
            --debug \
            --prune_centroids \
            --remove_singleton \
            .
    } || {
        # Check if the specific error occurs, and handle it
        if grep -q "AssertionError: No (specified) species pass the marker_depth filter" .command.log; then
            echo "Warning: No species passed the marker_depth filter. Skipping sample ${sample_id}."
            exit 0  # Continue with the next process without failure
        else
            echo "An unexpected error occurred." >&2
            exit 1  # Fail if some other error occurred
        fi
    }

    """

}



process RunGenes_PreBuiltBowtieDB {
    
    label 'mem_high'
    publishDir "${params.results_path}/", mode: "copy"
    maxRetries 2
    errorStrategy  { task.attempt <= maxRetries  ? 'retry' : 'ignore' }

    input:
    tuple val(sample_id), \
    path(R1), path(R2), \
    path(species_profile), \
    path(pangenomes_species)

    output:
    val("${sample_id}"), emit: sample_id
    path("${sample_id}/genes/log.txt")
    path("${sample_id}/genes/*.genes.tsv.lz4"), optional: true
    path("${sample_id}/genes/genes_summary.tsv"), optional: true
    
    """
    #! /usr/bin/env bash

    mkdir -p ${sample_id}/species/
    cp species_profile.tsv ${sample_id}/species/

    module load CBI samtools

    {
        midas run_genes \
            --sample_name ${sample_id} \
            -1 ${R1} -2 ${R2} \
            --midasdb_name ${params.db_name} \
            --midasdb_dir ${params.db_path} \
            --num_cores ${task.cpus} \
            --prebuilt_bowtie2_indexes "${params.db_path}/bt2_indexes/pangenomes" \
            --prebuilt_bowtie2_species ${pangenomes_species} \
            --select_by median_marker_coverage,unique_fraction_covered \
            --select_threshold 2,0.5 \
            --debug \
            --prune_centroids \
            --remove_singleton \
            .
    } || {
        # Check if the specific error occurs, and handle it
        if grep -q "AssertionError: No (specified) species pass the marker_depth filter" .command.log; then
            echo "Warning: No species passed the marker_depth filter. Skipping sample ${sample_id}."
            exit 0  # Continue with the next process without failure
        else
            echo "An unexpected error occurred." >&2
            exit 1  # Fail if some other error occurred
        fi
    }

    """

}

process MergeGenes {

    label 'mem_veryhigh'
    publishDir "${params.results_path}/merge/", mode: "copy"

    input:
    val(sample_id_list)

    output:
    path('genes/*/*.tsv.lz4')
    path('genes/genes_summary.tsv')
    path('genes/genes_log.txt')

    shell:
    """
    #! /usr/bin/env bash

    echo -e "sample_name\tmidas_outdir" > sample_list_for_gene_merge.tsv
    echo !{sample_id_list} | sed 's/[][]//g' | tr ',' '\n' | tr -d ' ' | \
    while read sample; do 
        if [ -f "!{params.results_path}/\$sample/genes/genes_summary.tsv" ]; then
            echo -e "\$sample\t!{params.results_path}"
        fi
    done >> sample_list_for_gene_merge.tsv



    midas merge_genes \
        --samples_list sample_list_for_gene_merge.tsv \
        --midasdb_name ${params.db_name} \
        --midasdb_dir ${params.db_path} \
        --num_cores ${task.cpus} \
	    --cluster_level_out ${params.merge_genes_cluster_level} \
        .
    """
}

