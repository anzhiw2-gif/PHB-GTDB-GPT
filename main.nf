#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

workflow {
    log.info "PHB-GTDB-GPT | selected stage: ${params.entry}"
    log.info "Project directory: ${params.project_dir}"
    log.info "Output directory: ${params.outdir}"
    log.info "Workflow modules will be enabled stage by stage after P01 input auditing."
}
