nextflow.enable.dsl=2

include { gwas_vcf_regenie_1 } from './modules/gwas_vcf_regenie_1/module.nf'
include { twas_vcf_1 } from './modules/twas_vcf_1/module.nf'

workflow {
input1 = Channel.fromPath(params.gwas_vcf_regenie_1.ch_user_input_vcf)
input2 = Channel.fromPath(params.gwas_vcf_regenie_1.ch_king_reference_data).splitCsv()
input3 = Channel.fromPath(params.gwas_vcf_regenie_1.ch_input_pheno_transform)
input4 = Channel.fromPath(params.gwas_vcf_regenie_1.ch_high_ld_regions)
input5 = Channel.fromPath(params.gwas_vcf_regenie_1.ch_gwas_cat)
input6 = Channel.fromPath(params.gwas_vcf_regenie_1.ch_ld_scores)
input7 = Channel.fromPath(params.gwas_vcf_regenie_1.ch_pheno)
input8 = Channel.fromPath(params.twas_vcf_1.ch_gene_annotations)
input9 = Channel.fromPath(params.twas_vcf_1.ch_codon_file)
input10 = Channel.fromPath(params.twas_vcf_1.ch_priority_file)
input11 = Channel.fromPath(params.twas_vcf_1.ch_ref_fasta)
input12 = Channel.fromPath(params.twas_vcf_1.ch_ref_fasta_index)
input13 = Channel.fromPath(params.twas_vcf_1.ch_ld_reference)
input14 = Channel.fromPath(params.twas_vcf_1.ch_eqtl_weights)
input15 = Channel.fromPath(params.twas_vcf_1.ch_gwas_sumstats)
gwas_vcf_regenie_1(input1, input2, input3, input4, input5, input6, input7)
twas_vcf_1(input8, input9, input10, input11, input12, input13, input14, input15)
}
