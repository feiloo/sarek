
// Create samplesheets to restart from different steps
include { CHANNEL_ALIGN_CREATE_CSV                    } from '../subworkflows/local/channel_align_create_csv/main'
include { CHANNEL_MARKDUPLICATES_CREATE_CSV           } from '../subworkflows/local/channel_markduplicates_create_csv/main'
include { CHANNEL_BASERECALIBRATOR_CREATE_CSV         } from '../subworkflows/local/channel_baserecalibrator_create_csv/main'
include { CHANNEL_APPLYBQSR_CREATE_CSV                } from '../subworkflows/local/channel_applybqsr_create_csv/main'
include { CHANNEL_VARIANT_CALLING_CREATE_CSV          } from '../subworkflows/local/channel_variant_calling_create_csv/main'

// Download annotation cache if needed
include { PREPARE_CACHE                               } from '../subworkflows/local/prepare_cache/main'

// Build indices if needed
include { PREPARE_GENOME                              } from '../subworkflows/local/prepare_genome/main'

// Build intervals if needed
include { PREPARE_INTERVALS                           } from '../subworkflows/local/prepare_intervals/main'

// Build CNVkit reference if needed
include { PREPARE_REFERENCE_CNVKIT                    } from '../subworkflows/local/prepare_reference_cnvkit/main'

// Convert BAM files to FASTQ files
include { BAM_CONVERT_SAMTOOLS as CONVERT_FASTQ_INPUT } from '../subworkflows/local/bam_convert_samtools/main'
include { BAM_CONVERT_SAMTOOLS as CONVERT_FASTQ_UMI   } from '../subworkflows/local/bam_convert_samtools/main'

// Run FASTQC
include { FASTQC                                      } from '../modules/nf-core/fastqc/main'

// TRIM/SPLIT FASTQ Files
include { FASTP                                       } from '../modules/nf-core/fastp/main'

// Create umi consensus bams from fastq
include { FASTQ_CREATE_UMI_CONSENSUS_FGBIO            } from '../subworkflows/local/fastq_create_umi_consensus_fgbio/main'

// Map input reads to reference genome
include { FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON    } from '../subworkflows/local/fastq_align_bwamem_mem2_dragmap_sentieon/main'

// Merge and index BAM files (optional)
include { BAM_MERGE_INDEX_SAMTOOLS                    } from '../subworkflows/local/bam_merge_index_samtools/main'

// Convert BAM files
include { SAMTOOLS_CONVERT as BAM_TO_CRAM             } from '../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_CONVERT as BAM_TO_CRAM_MAPPING     } from '../modules/nf-core/samtools/convert/main'

// Convert CRAM files (optional)
include { SAMTOOLS_CONVERT as CRAM_TO_BAM             } from '../modules/nf-core/samtools/convert/main'
include { SAMTOOLS_CONVERT as CRAM_TO_BAM_RECAL       } from '../modules/nf-core/samtools/convert/main'

// Mark Duplicates (+QC)
include { BAM_MARKDUPLICATES                          } from '../subworkflows/local/bam_markduplicates/main'
include { BAM_MARKDUPLICATES_SPARK                    } from '../subworkflows/local/bam_markduplicates_spark/main'
include { BAM_SENTIEON_DEDUP                          } from '../subworkflows/local/bam_sentieon_dedup/main'

// QC on CRAM
include { CRAM_QC_MOSDEPTH_SAMTOOLS as CRAM_QC_NO_MD  } from '../subworkflows/local/cram_qc_mosdepth_samtools/main'
include { CRAM_QC_MOSDEPTH_SAMTOOLS as CRAM_QC_RECAL  } from '../subworkflows/local/cram_qc_mosdepth_samtools/main'

// Create recalibration tables
include { BAM_BASERECALIBRATOR                        } from '../subworkflows/local/bam_baserecalibrator/main'
include { BAM_BASERECALIBRATOR_SPARK                  } from '../subworkflows/local/bam_baserecalibrator_spark/main'

// Create recalibrated cram files to use for variant calling (+QC)
include { BAM_APPLYBQSR                               } from '../subworkflows/local/bam_applybqsr/main'
include { BAM_APPLYBQSR_SPARK                         } from '../subworkflows/local/bam_applybqsr_spark/main'

// Variant calling on a single normal sample
include { BAM_VARIANT_CALLING_GERMLINE_ALL            } from '../subworkflows/local/bam_variant_calling_germline_all/main'

// Variant calling on a single tumor sample
include { BAM_VARIANT_CALLING_TUMOR_ONLY_ALL          } from '../subworkflows/local/bam_variant_calling_tumor_only_all/main'

// Variant calling on tumor/normal pair
include { BAM_VARIANT_CALLING_SOMATIC_ALL             } from '../subworkflows/local/bam_variant_calling_somatic_all/main'

// POST VARIANTCALLING: e.g. merging
include { POST_VARIANTCALLING                         } from '../subworkflows/local/post_variantcalling/main'

// QC on VCF files
include { VCF_QC_BCFTOOLS_VCFTOOLS                    } from '../subworkflows/local/vcf_qc_bcftools_vcftools/main'

// Annotation
include { VCF_ANNOTATE_ALL                            } from '../subworkflows/local/vcf_annotate_all/main'

// REPORTING VERSIONS OF SOFTWARE USED
include { CUSTOM_DUMPSOFTWAREVERSIONS                 } from '../modules/nf-core/custom/dumpsoftwareversions/main'

// MULTIQC
include { MULTIQC                                     } from '../modules/nf-core/multiqc/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
nextflow.enable.dsl = 2

workflow SAREK {
take:
  args

main:
  println args.dbsnp
  println params.dbsnp


// Check input path parameters to see if they exist
def checkPathParamList = [
    args.ascat_alleles,
    args.ascat_loci,
    args.ascat_loci_gc,
    args.ascat_loci_rt,
    args.bwa,
    args.bwamem2,
    args.cf_chrom_len,
    args.chr_dir,
    args.cnvkit_reference,
    args.dbnsfp,
    args.dbnsfp_tbi,
    args.dbsnp,
    args.dbsnp_tbi,
    args.dict,
    args.dragmap,
    args.fasta,
    args.fasta_fai,
    args.germline_resource,
    args.germline_resource_tbi,
    args.input,
    args.intervals,
    args.known_indels,
    args.known_indels_tbi,
    args.known_snps,
    args.known_snps_tbi,
    args.mappability,
    args.multiqc_config,
    args.pon,
    args.pon_tbi,
    args.spliceai_indel,
    args.spliceai_indel_tbi,
    args.spliceai_snv,
    args.spliceai_snv_tbi
]

// only check if we are using the tools
if (args.tools && (args.tools.split(',').contains('snpeff') || args.tools.split(',').contains('merge'))) checkPathParamList.add(args.snpeff_cache)
if (args.tools && (args.tools.split(',').contains('vep')    || args.tools.split(',').contains('merge'))) checkPathParamList.add(args.vep_cache)


for (param in checkPathParamList) if (param) file(param, checkIfExists: true)

// Set input, can either be from --input or from automatic retrieval in WorkflowSarek.groovy

ch_from_samplesheet = args.build_only_index ? Channel.empty() : Channel.fromPath(args.input).splitCsv(skip: 1)

input_sample = ch_from_samplesheet
        .map{ meta, fastq_1, fastq_2, table, cram, crai, bam, bai, vcf, variantcaller ->
            // generate patient_sample key to group lanes together
            [ meta.patient + meta.sample, [meta, fastq_1, fastq_2, table, cram, crai, bam, bai, vcf, variantcaller] ]
        }
        .tap{ ch_with_patient_sample } // save the channel
        .groupTuple() //group by patient_sample to get all lanes
        .map { patient_sample, ch_items ->
            // get number of lanes per sample
            [ patient_sample, ch_items.size() ]
        }
        .combine(ch_with_patient_sample, by: 0) // for each entry add numLanes
        .map { patient_sample, num_lanes, ch_items ->

            (meta, fastq_1, fastq_2, table, cram, crai, bam, bai, vcf, variantcaller) = ch_items
            if (meta.lane && fastq_2) {
                meta           = meta + [id: "${meta.sample}-${meta.lane}".toString()]
                def CN         = args.seq_center ? "CN:${args.seq_center}\\t" : ''

                def flowcell   = flowcellLaneFromFastq(fastq_1)
                // Don't use a random element for ID, it breaks resuming
                def read_group = "\"@RG\\tID:${flowcell}.${meta.sample}.${meta.lane}\\t${CN}PU:${meta.lane}\\tSM:${meta.patient}_${meta.sample}\\tLB:${meta.sample}\\tDS:${args.fasta}\\tPL:${args.seq_platform}\""

                meta           = meta - meta.subMap('lane') + [num_lanes: num_lanes.toInteger(), read_group: read_group.toString(), data_type: 'fastq', size: 1]

                if (args.step == 'mapping') return [ meta, [ fastq_1, fastq_2 ] ]
                else {
                    error("Samplesheet contains fastq files but step is `$args.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations")
                }

            // start from BAM
            } else if (meta.lane && bam) {
                if (args.step != 'mapping' && !bai) {
                    error("BAM index (bai) should be provided.")
                }
                meta            = meta + [id: "${meta.sample}-${meta.lane}".toString()]
                def CN          = args.seq_center ? "CN:${args.seq_center}\\t" : ''
                def read_group  = "\"@RG\\tID:${meta.sample}_${meta.lane}\\t${CN}PU:${meta.lane}\\tSM:${meta.patient}_${meta.sample}\\tLB:${meta.sample}\\tDS:${args.fasta}\\tPL:${args.seq_platform}\""

                meta            = meta - meta.subMap('lane') + [num_lanes: num_lanes.toInteger(), read_group: read_group.toString(), data_type: 'bam', size: 1]

                if (args.step != 'annotate') return [ meta - meta.subMap('lane'), bam, bai ]
                else {
                    error("Samplesheet contains bam files but step is `annotate`. The pipeline is expecting vcf files for the annotation. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations")
                }

            // recalibration
            } else if (table && cram) {
                meta = meta + [id: meta.sample, data_type: 'cram']

                if (!(args.step == 'mapping' || args.step == 'annotate')) return [ meta - meta.subMap('lane'), cram, crai, table ]
                else {
                    error("Samplesheet contains cram files but step is `$args.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations")
                }

            // recalibration when skipping MarkDuplicates
            } else if (table && bam) {
                meta = meta + [id: meta.sample, data_type: 'bam']

                if (!(args.step == 'mapping' || args.step == 'annotate')) return [ meta - meta.subMap('lane'), bam, bai, table ]
                else {
                    error("Samplesheet contains bam files but step is `$args.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations")
                }

            // prepare_recalibration or variant_calling
            } else if (cram) {
                meta = meta + [id: meta.sample, data_type: 'cram']

                if (!(args.step == 'mapping' || args.step == 'annotate')) return [ meta - meta.subMap('lane'), cram, crai ]
                else {
                    error("Samplesheet contains cram files but step is `$args.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations")
                }

            // prepare_recalibration when skipping MarkDuplicates or `--step markduplicates`
            } else if (bam) {
                meta = meta + [id: meta.sample, data_type: 'bam']

                if (!(args.step == 'mapping' || args.step == 'annotate')) return [ meta - meta.subMap('lane'), bam, bai ]
                else {
                    error("Samplesheet contains bam files but step is `$args.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations")
                }

            // annotation
            } else if (vcf) {
                meta = meta + [id: meta.sample, data_type: 'vcf', variantcaller: variantcaller ?: '']

                if (args.step == 'annotate') return [ meta - meta.subMap('lane'), vcf ]
                else {
                    error("Samplesheet contains vcf files but step is `$args.step`. Please check your samplesheet or adjust the step parameter.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations")
                }
            } else {
                error("Missing or unknown field in csv file header. Please check your samplesheet")
            }
        }

if (args.step != 'annotate' && args.tools && !args.build_only_index) {
    // Two checks for ensuring that the pipeline stops with a meaningful error message if
    // 1. the sample-sheet only contains normal-samples, but some of the requested tools require tumor-samples, and
    // 2. the sample-sheet only contains tumor-samples, but some of the requested tools require normal-samples.
    input_sample.filter{ it[0].status == 1 }.ifEmpty{ // In this case, the sample-sheet contains no tumor-samples
        if (!args.build_only_index) {
            def tools_tumor = ['ascat', 'controlfreec', 'mutect2', 'msisensorpro']
            def tools_tumor_asked = []
            tools_tumor.each{ tool ->
                if (args.tools.split(',').contains(tool)) tools_tumor_asked.add(tool)
            }
            if (!tools_tumor_asked.isEmpty()) {
                error('The sample-sheet only contains normal-samples, but the following tools, which were requested with "--tools", expect at least one tumor-sample : ' + tools_tumor_asked.join(", "))
            }
        }
    }
    input_sample.filter{ it[0].status == 0 }.ifEmpty{ // In this case, the sample-sheet contains no normal/germline-samples
        def tools_requiring_normal_samples = ['ascat', 'deepvariant', 'haplotypecaller', 'msisensorpro']
        def requested_tools_requiring_normal_samples = []
        tools_requiring_normal_samples.each{ tool_requiring_normal_samples ->
            if (args.tools.split(',').contains(tool_requiring_normal_samples)) requested_tools_requiring_normal_samples.add(tool_requiring_normal_samples)
        }
        if (!requested_tools_requiring_normal_samples.isEmpty()) {
            error('The sample-sheet only contains tumor-samples, but the following tools, which were requested by the option "tools", expect at least one normal-sample : ' + requested_tools_requiring_normal_samples.join(", "))
        }
    }
}

// Fails when wrongfull extension for intervals file
if (args.wes && !args.step == 'annotate') {
    if (args.intervals && !args.intervals.endsWith("bed"))  error("Target file specified with `--intervals` must be in BED format for targeted data")
    else log.warn("Intervals file was provided without parameter `--wes`: Pipeline will assume this is Whole-Genome-Sequencing data.")
} else if (args.intervals && !args.intervals.endsWith("bed") && !args.intervals.endsWith("list")) error("Intervals file must end with .bed, .list, or .interval_list")

if (args.step == 'mapping' && args.aligner.contains("dragmap") && !(args.skip_tools && args.skip_tools.split(',').contains("baserecalibrator"))) {
    log.warn("DragMap was specified as aligner. Base recalibration is not contained in --skip_tools. It is recommended to skip baserecalibration when using DragMap\nhttps://gatk.broadinstitute.org/hc/en-us/articles/4407897446939--How-to-Run-germline-single-sample-short-variant-discovery-in-DRAGEN-mode")
}

if (args.step == 'mapping' && args.aligner.contains("sentieon-bwamem") && args.umi_read_structure) {
    error("Sentieon BWA is currently not compatible with FGBio UMI handeling. Please choose a different aligner.")
}

if (args.tools && args.tools.split(',').contains("sentieon_haplotyper") && args.joint_germline && (!args.sentieon_haplotyper_emit_mode || !(args.sentieon_haplotyper_emit_mode.contains('gvcf')))) {
    error("When setting the option `--joint_germline` and including `sentieon_haplotyper` among the requested tools, please set `--sentieon_haplotyper_emit_mode` to include `gvcf`.")
}

// Fails or warns when missing files or args.for ascat
if (args.tools && args.tools.split(',').contains('ascat')) {
    if (!args.ascat_alleles) {
        error("No allele files were provided for running ASCAT. Please provide a zip folder with allele files.")
    }
    if (!args.ascat_loci) {
        error("No loci files were provided for running ASCAT. Please provide a zip folder with loci files.")
    }
    if (!args.ascat_loci_gc && !args.ascat_loci_rt) {
        log.warn("No LogRCorrection performed in ASCAT. For LogRCorrection to run, please provide either loci gc files or both loci gc files and loci rt files.")
    }
    if (args.wes) {
        log.warn("Default reference files not suited for running ASCAT on WES data. It's recommended to use the reference files provided here: https://github.com/Wedge-lab/battenberg#required-reference-files")
    }
}

// Warns when missing files or args.for mutect2
if (args.tools && args.tools.split(',').contains('mutect2')) {
    if (!args.pon) {
        log.warn("No Panel-of-normal was specified for Mutect2.\nIt is highly recommended to use one: https://gatk.broadinstitute.org/hc/en-us/articles/5358911630107-Mutect2\nFor more information on how to create one: https://gatk.broadinstitute.org/hc/en-us/articles/5358921041947-CreateSomaticPanelOfNormals-BETA-")
    }
    if (!args.germline_resource) {
        log.warn("If Mutect2 is specified without a germline resource, no filtering will be done.\nIt is recommended to use one: https://gatk.broadinstitute.org/hc/en-us/articles/5358911630107-Mutect2")
    }
    if (args.pon && args.pon.contains("/Homo_sapiens/GATK/GRCh38/Annotation/GATKBundle/1000g_pon.hg38.vcf.gz")) {
        log.warn("The default Panel-of-Normals provided by GATK is used for Mutect2.\nIt is highly recommended to generate one from normal samples that are technical similar to the tumor ones.\nFor more information: https://gatk.broadinstitute.org/hc/en-us/articles/360035890631-Panel-of-Normals-PON-")
    }
}

// Fails when missing resources for baserecalibrator
// Warns when missing resources for haplotypecaller
if (!args.dbsnp && !args.known_indels) {
    if (args.step in ['mapping', 'markduplicates', 'prepare_recalibration', 'recalibrate'] && (!args.skip_tools || (args.skip_tools && !args.skip_tools.split(',').contains('baserecalibrator')))) {
        error("Base quality score recalibration requires at least one resource file. Please provide at least one of `--dbsnp` or `--known_indels`\nYou can skip this step in the workflow by adding `--skip_tools baserecalibrator` to the command.")
    }
    if (args.tools && (args.tools.split(',').contains('haplotypecaller') || args.tools.split(',').contains('sentieon_haplotyper'))) {
        log.warn "If GATK's Haplotypecaller or Sentieon's Haplotyper is specified, without `--dbsnp` or `--known_indels no filtering will be done. For filtering, please provide at least one of `--dbsnp` or `--known_indels`.\nFor more information see FilterVariantTranches (single-sample, default): https://gatk.broadinstitute.org/hc/en-us/articles/5358928898971-FilterVariantTranches\nFor more information see VariantRecalibration (--joint_germline): https://gatk.broadinstitute.org/hc/en-us/articles/5358906115227-VariantRecalibrator\nFor more information on GATK Best practice germline variant calling: https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-"
    }
}
if (args.joint_germline && (!args.tools || !(args.tools.split(',').contains('haplotypecaller') || args.tools.split(',').contains('sentieon_haplotyper')))) {
    error("The GATK's Haplotypecaller or Sentieon's Haplotyper should be specified as one of the tools when doing joint germline variant calling.) ")
}

if (args.joint_germline && (!args.dbsnp || !args.known_indels || !args.known_snps || args.no_intervals)) {
    log.warn "If GATK's Haplotypecaller or Sentieon's Haplotyper is specified, without `--dbsnp`, `--known_snps`, `--known_indels` or the associated resource labels (ie `known_snps_vqsr`), no variant recalibration will be done. For recalibration you must provide all of these resources.\nFor more information see VariantRecalibration: https://gatk.broadinstitute.org/hc/en-us/articles/5358906115227-VariantRecalibrator \nJoint germline variant calling also requires intervals in order to genotype the samples. As a result, if `--no_intervals` is set to `true` the joint germline variant calling will not be performed."
}

// Fails when --joint_mutect2 is used without enabling mutect2
if (args.joint_mutect2 && (!args.tools || !args.tools.split(',').contains('mutect2'))) {
    error("The mutect2 should be specified as one of the tools when doing joint somatic variant calling with Mutect2. (The mutect2 could be specified by adding `--tools mutect2` to the nextflow command.)")
}

// Fails when missing tools for variant_calling or annotate
if ((args.step == 'variant_calling' || args.step == 'annotate') && !args.tools) {
    error("Please specify at least one tool when using `--step ${args.step}`.\nhttps://nf-co.re/sarek/parameters#tools")
}

// Fails when missing sex information for CNV tools
if (args.tools && (args.tools.split(',').contains('ascat') || args.tools.split(',').contains('controlfreec'))) {
    input_sample.map{
        if (it[0].sex == 'NA' ) {
            error("Please specify sex information for each sample in your samplesheet when using '--tools' with 'ascat' or 'controlfreec'.\nhttps://nf-co.re/sarek/usage#input-samplesheet-configurations")
        }
    }
}

if ((args.download_cache) && (args.snpeff_cache || args.vep_cache)) {
    error("Please specify either `--download_cache` or `--snpeff_cache`, `--vep_cache`.\nhttps://nf-co.re/sarek/usage#how-to-customise-snpeff-and-vep-annotation")
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Initialize file channels based on args. defined in the args.genomes[args.genome] scope
ascat_alleles      = args.ascat_alleles      ? Channel.fromPath(args.ascat_alleles).collect()     : Channel.empty()
ascat_loci         = args.ascat_loci         ? Channel.fromPath(args.ascat_loci).collect()        : Channel.empty()
ascat_loci_gc      = args.ascat_loci_gc      ? Channel.fromPath(args.ascat_loci_gc).collect()     : Channel.value([])
ascat_loci_rt      = args.ascat_loci_rt      ? Channel.fromPath(args.ascat_loci_rt).collect()     : Channel.value([])
cf_chrom_len       = args.cf_chrom_len       ? Channel.fromPath(args.cf_chrom_len).collect()      : []
chr_dir            = args.chr_dir            ? Channel.fromPath(args.chr_dir).collect()           : Channel.value([])
dbsnp              = args.dbsnp              ? Channel.fromPath(args.dbsnp).collect()             : Channel.value([])
fasta              = args.fasta              ? Channel.fromPath(args.fasta).first()               : Channel.empty()
fasta_fai          = args.fasta_fai          ? Channel.fromPath(args.fasta_fai).collect()         : Channel.empty()
germline_resource  = args.germline_resource  ? Channel.fromPath(args.germline_resource).collect() : Channel.value([]) // Mutect2 does not require a germline resource, so set to optional input
known_indels       = args.known_indels       ? Channel.fromPath(args.known_indels).collect()      : Channel.value([])
known_snps         = args.known_snps         ? Channel.fromPath(args.known_snps).collect()        : Channel.value([])
mappability        = args.mappability        ? Channel.fromPath(args.mappability).collect()       : Channel.value([])
pon                = args.pon                ? Channel.fromPath(args.pon).collect()               : Channel.value([]) // PON is optional for Mutect2 (but highly recommended)

// Initialize value channels based on args. defined in the args.genomes[args.genome] scope
ascat_genome       = args.ascat_genome       ?: Channel.empty()
dbsnp_vqsr         = args.dbsnp_vqsr         ? Channel.value(args.dbsnp_vqsr) : Channel.empty()
known_indels_vqsr  = args.known_indels_vqsr  ? Channel.value(args.known_indels_vqsr) : Channel.empty()
known_snps_vqsr    = args.known_snps_vqsr    ? Channel.value(args.known_snps_vqsr) : Channel.empty()
snpeff_db          = args.snpeff_db          ?: Channel.empty()
vep_cache_version  = args.vep_cache_version  ?: Channel.empty()
vep_genome         = args.vep_genome         ?: Channel.empty()
vep_species        = args.vep_species        ?: Channel.empty()

// Initialize files channels based on args. not defined within the args.genomes[args.genome] scope
if (args.snpeff_cache && args.tools && (args.tools.split(',').contains("snpeff") || args.tools.split(',').contains('merge'))) {
    def snpeff_annotation_cache_key = ''
    if (args.snpeff_cache == "s3://annotation-cache/snpeff_cache") {
        snpeff_annotation_cache_key = "${args.snpeff_genome}.${args.snpeff_db}/"
    } else {
        snpeff_annotation_cache_key = args.use_annotation_cache_keys ? "${args.snpeff_genome}.${args.snpeff_db}/" : ""
    }
    def snpeff_cache_dir =  "${snpeff_annotation_cache_key}${args.snpeff_genome}.${args.snpeff_db}"
    def snpeff_cache_path_full = file("$args.snpeff_cache/$snpeff_cache_dir", type: 'dir')
    if ( !snpeff_cache_path_full.exists() || !snpeff_cache_path_full.isDirectory() ) {
        if (args.snpeff_cache == "s3://annotation-cache/snpeff_cache") {
            error("This path is not available within annotation-cache. Please check https://annotation-cache.github.io/ to create a request for it.")
        } else {
            error("Files within --snpeff_cache invalid. Make sure there is a directory named ${snpeff_cache_dir} in ${args.snpeff_cache}.\nhttps://nf-co.re/sarek/usage#how-to-customise-snpeff-and-vep-annotation")
        }
    }
    snpeff_cache = Channel.fromPath(file("${args.snpeff_cache}/${snpeff_annotation_cache_key}"), checkIfExists: true).collect()
        .map{ cache -> [ [ id:"${args.snpeff_genome}.${args.snpeff_db}" ], cache ] }
    } else if (args.tools && (args.tools.split(',').contains("snpeff") || args.tools.split(',').contains('merge')) && !args.download_cache) {
        error("No cache for SnpEff or automatic download of said cache has been detected.\nPlease refer to https://nf-co.re/sarek/docs/usage/#how-to-customise-snpeff-and-vep-annotation for more information.")
    } else snpeff_cache = []

if (args.vep_cache && args.tools && (args.tools.split(',').contains("vep") || args.tools.split(',').contains('merge'))) {
    def vep_annotation_cache_key = ''
    if (args.vep_cache == "s3://annotation-cache/vep_cache") {
        vep_annotation_cache_key = "${args.vep_cache_version}_${args.vep_genome}/"
    } else {
        vep_annotation_cache_key = args.use_annotation_cache_keys ? "${args.vep_cache_version}_${args.vep_genome}/" : ""
    }
    def vep_cache_dir = "${vep_annotation_cache_key}${args.vep_species}/${args.vep_cache_version}_${args.vep_genome}"
    def vep_cache_path_full = file("$args.vep_cache/$vep_cache_dir", type: 'dir')
    if ( !vep_cache_path_full.exists() || !vep_cache_path_full.isDirectory() ) {
        if (args.vep_cache == "s3://annotation-cache/vep_cache") {
            error("This path is not available within annotation-cache. Please check https://annotation-cache.github.io/ to create a request for it.")
        } else {
            error("Files within --vep_cache invalid. Make sure there is a directory named ${vep_cache_dir} in ${args.vep_cache}.\nhttps://nf-co.re/sarek/usage#how-to-customise-snpeff-and-vep-annotation")
        }
    }
    vep_cache = Channel.fromPath(file("${args.vep_cache}/${vep_annotation_cache_key}"), checkIfExists: true).collect()
    } else if (args.tools && (args.tools.split(',').contains("vep") || args.tools.split(',').contains('merge')) && !args.download_cache) {
        error("No cache for VEP or automatic download of said cache has been detected.\nPlease refer to https://nf-co.re/sarek/docs/usage/#how-to-customise-snpeff-and-vep-annotation for more information.")
    } else vep_cache = []

vep_extra_files = []

if (args.dbnsfp && args.dbnsfp_tbi) {
    vep_extra_files.add(file(args.dbnsfp, checkIfExists: true))
    vep_extra_files.add(file(args.dbnsfp_tbi, checkIfExists: true))
}

if (args.spliceai_snv && args.spliceai_snv_tbi && args.spliceai_indel && args.spliceai_indel_tbi) {
    vep_extra_files.add(file(args.spliceai_indel, checkIfExists: true))
    vep_extra_files.add(file(args.spliceai_indel_tbi, checkIfExists: true))
    vep_extra_files.add(file(args.spliceai_snv, checkIfExists: true))
    vep_extra_files.add(file(args.spliceai_snv_tbi, checkIfExists: true))
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL/NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/



    // MULTIQC
    ch_multiqc_config                     = Channel.fromPath("$moduleDir/../assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config              = args.multiqc_config ? Channel.fromPath( args.multiqc_config, checkIfExists: true ) : Channel.empty()
    ch_multiqc_logo                       = args.multiqc_logo   ? Channel.fromPath( args.multiqc_logo, checkIfExists: true ) : Channel.empty()
    ch_multiqc_custom_methods_description = args.multiqc_methods_description ? file(args.multiqc_methods_description, checkIfExists: true) : file("$moduleDir/../assets/methods_description_template.yml", checkIfExists: true)

    // To gather all QC reports for MultiQC
    reports  = Channel.empty()
    // To gather used softwares versions for MultiQC
    versions = Channel.empty()

    // Download cache if needed
    // Assuming that if the cache is provided, the user has already downloaded it
    ensemblvep_info = args.vep_cache    ? [] : Channel.of([ [ id:"${args.vep_cache_version}_${args.vep_genome}" ], args.vep_genome, args.vep_species, args.vep_cache_version ])
    snpeff_info     = args.snpeff_cache ? [] : Channel.of([ [ id:"${args.snpeff_genome}.${args.snpeff_db}" ], args.snpeff_genome, args.snpeff_db ])

    if (args.download_cache) {
        PREPARE_CACHE(ensemblvep_info, snpeff_info)
        snpeff_cache = PREPARE_CACHE.out.snpeff_cache
        vep_cache    = PREPARE_CACHE.out.ensemblvep_cache.map{ meta, cache -> [ cache ] }

        versions = versions.mix(PREPARE_CACHE.out.versions)
    }

    // Build indices if needed
    PREPARE_GENOME(
        ascat_alleles,
        ascat_loci,
        ascat_loci_gc,
        ascat_loci_rt,
        chr_dir,
        dbsnp,
        fasta,
        fasta_fai,
        germline_resource,
        known_indels,
        known_snps,
        pon)

    // Gather built indices or get them from the params
    // Built from the fasta file:
    dict       = args.dict        ? Channel.fromPath(args.dict).map{ it -> [ [id:'dict'], it ] }.collect()
                                    : PREPARE_GENOME.out.dict
    fasta_fai  = args.fasta_fai   ? Channel.fromPath(args.fasta_fai).collect()
                                    : PREPARE_GENOME.out.fasta_fai
    bwa        = args.bwa         ? Channel.fromPath(args.bwa).collect()
                                    : PREPARE_GENOME.out.bwa
    bwamem2    = args.bwamem2     ? Channel.fromPath(args.bwamem2).collect()
                                    : PREPARE_GENOME.out.bwamem2
    dragmap    = args.dragmap     ? Channel.fromPath(args.dragmap).collect()
                                    : PREPARE_GENOME.out.hashtable

    // Gather index for mapping given the chosen aligner
    index_alignement = (args.aligner == "bwa-mem" || args.aligner == "sentieon-bwamem") ? bwa :
        args.aligner == "bwa-mem2" ? bwamem2 :
        dragmap

    // TODO: add a args.for msisensorpro_scan
    msisensorpro_scan      = PREPARE_GENOME.out.msisensorpro_scan

    // For ASCAT, extracted from zip or tar.gz files:
    allele_files           = PREPARE_GENOME.out.allele_files
    chr_files              = PREPARE_GENOME.out.chr_files
    gc_file                = PREPARE_GENOME.out.gc_file
    loci_files             = PREPARE_GENOME.out.loci_files
    rt_file                = PREPARE_GENOME.out.rt_file

    // Tabix indexed vcf files:
    dbsnp_tbi              = args.dbsnp                   ? args.dbsnp_tbi             ? Channel.fromPath(args.dbsnp_tbi).collect()             : PREPARE_GENOME.out.dbsnp_tbi             : Channel.value([])
    germline_resource_tbi  = args.germline_resource       ? args.germline_resource_tbi ? Channel.fromPath(args.germline_resource_tbi).collect() : PREPARE_GENOME.out.germline_resource_tbi : [] //do not change to Channel.value([]), the check for its existence then fails for Getpileupsumamries
    known_indels_tbi       = args.known_indels            ? args.known_indels_tbi      ? Channel.fromPath(args.known_indels_tbi).collect()      : PREPARE_GENOME.out.known_indels_tbi      : Channel.value([])
    known_snps_tbi         = args.known_snps              ? args.known_snps_tbi        ? Channel.fromPath(args.known_snps_tbi).collect()        : PREPARE_GENOME.out.known_snps_tbi        : Channel.value([])
    pon_tbi                = args.pon                     ? args.pon_tbi               ? Channel.fromPath(args.pon_tbi).collect()               : PREPARE_GENOME.out.pon_tbi               : Channel.value([])

    // known_sites is made by grouping both the dbsnp and the known snps/indels resources
    // Which can either or both be optional
    known_sites_indels     = dbsnp.concat(known_indels).collect()
    known_sites_indels_tbi = dbsnp_tbi.concat(known_indels_tbi).collect()

    known_sites_snps       = dbsnp.concat(known_snps).collect()
    known_sites_snps_tbi   = dbsnp_tbi.concat(known_snps_tbi).collect()

    // Build intervals if needed
    PREPARE_INTERVALS(fasta_fai, args.intervals, args.no_intervals)

    // Intervals for speed up preprocessing/variant calling by spread/gather
    // [interval.bed] all intervals in one file
    intervals_bed_combined         = args.no_intervals ? Channel.value([])      : PREPARE_INTERVALS.out.intervals_bed_combined
    intervals_bed_gz_tbi_combined  = args.no_intervals ? Channel.value([])      : PREPARE_INTERVALS.out.intervals_bed_gz_tbi_combined

    // For QC during preprocessing, we don't need any intervals (MOSDEPTH doesn't take them for WGS)
    intervals_for_preprocessing = args.wes ?
        intervals_bed_combined.map{it -> [ [ id:it.baseName ], it ]}.collect() :
        Channel.value([ [ id:'null' ], [] ])

    intervals            = PREPARE_INTERVALS.out.intervals_bed        // [ interval, num_intervals ] multiple interval.bed files, divided by useful intervals for scatter/gather
    intervals_bed_gz_tbi = PREPARE_INTERVALS.out.intervals_bed_gz_tbi // [ interval_bed, tbi, num_intervals ] multiple interval.bed.gz/.tbi files, divided by useful intervals for scatter/gather

    intervals_and_num_intervals = intervals.map{ interval, num_intervals ->
        if ( num_intervals < 1 ) [ [], num_intervals ]
        else [ interval, num_intervals ]
    }

    intervals_bed_gz_tbi_and_num_intervals = intervals_bed_gz_tbi.map{ intervals, num_intervals ->
        if ( num_intervals < 1 ) [ [], [], num_intervals ]
        else [ intervals[0], intervals[1], num_intervals ]
    }

    if (args.tools && args.tools.split(',').contains('cnvkit')) {
        if (args.cnvkit_reference) {
            cnvkit_reference = Channel.fromPath(args.cnvkit_reference).collect()
        } else {
            PREPARE_REFERENCE_CNVKIT(fasta, intervals_bed_combined)
            cnvkit_reference = PREPARE_REFERENCE_CNVKIT.out.cnvkit_reference

            versions = versions.mix(PREPARE_REFERENCE_CNVKIT.out.versions)
        }
    } else {
        cnvkit_reference = Channel.value([])
    }

    // Gather used softwares versions
    versions = versions.mix(PREPARE_GENOME.out.versions)
    versions = versions.mix(PREPARE_INTERVALS.out.versions)

    // PREPROCESSING

    if (args.step == 'mapping') {

        // Figure out if input is bam or fastq
        input_sample_type = input_sample.branch{
            bam:   it[0].data_type == "bam"
            fastq: it[0].data_type == "fastq"
        }

        // Convert any bam input to fastq
        // fasta are not needed when converting bam to fastq -> [ id:"fasta" ], []
        // No need for fasta.fai -> []
        interleave_input = false // Currently don't allow interleaved input
        CONVERT_FASTQ_INPUT(
            input_sample_type.bam,
            [ [ id:"fasta" ], [] ], // fasta
            [ [ id:'null' ], [] ],  // fasta_fai
            interleave_input)

        // Gather fastq (inputed or converted)
        // Theorically this could work on mixed input (fastq for one sample and bam for another)
        // But not sure how to handle that with the samplesheet
        // Or if we really want users to be able to do that
        input_fastq = input_sample_type.fastq.mix(CONVERT_FASTQ_INPUT.out.reads)

        // STEP 0: QC & TRIM
        // `--skip_tools fastqc` to skip fastqc
        // Trim only with `--trim_fastq`
        // Additional options to be set up

        // QC
        if (!(args.skip_tools && args.skip_tools.split(',').contains('fastqc'))) {
            FASTQC(input_fastq)

            reports = reports.mix(FASTQC.out.zip.collect{ meta, logs -> logs })
            versions = versions.mix(FASTQC.out.versions.first())
        }

        // UMI consensus calling
        if (args.umi_read_structure) {
            FASTQ_CREATE_UMI_CONSENSUS_FGBIO(
                input_fastq,
                fasta,
                fasta_fai,
                index_alignement,
                args.group_by_umi_strategy)

            bam_converted_from_fastq = FASTQ_CREATE_UMI_CONSENSUS_FGBIO.out.consensusbam.map{ meta, bam -> [ meta, bam, [] ] }

            // Convert back to fastq for further preprocessing
            // fasta are not needed when converting bam to fastq -> [ id:"fasta" ], []
            // No need for fasta.fai -> []
            interleave_input = false // Currently don't allow interleaved input
            CONVERT_FASTQ_UMI(
                bam_converted_from_fastq,
                [ [ id:"fasta" ], [] ], // fasta
                [ [ id:'null' ], [] ],  // fasta_fai
                interleave_input)

            reads_for_fastp = CONVERT_FASTQ_UMI.out.reads

            // Gather used softwares versions
            versions = versions.mix(CONVERT_FASTQ_UMI.out.versions)
            versions = versions.mix(FASTQ_CREATE_UMI_CONSENSUS_FGBIO.out.versions)
        } else {
            reads_for_fastp = input_fastq
        }

        // Trimming and/or splitting
        if (args.trim_fastq || args.split_fastq > 0) {

            save_trimmed_fail = false
            save_merged = false
            FASTP(
                reads_for_fastp,
                [], // we are not using any adapter fastas at the moment
                save_trimmed_fail,
                save_merged
            )

            reports = reports.mix(FASTP.out.json.collect{ meta, json -> json })
            reports = reports.mix(FASTP.out.html.collect{ meta, html -> html })

            if (args.split_fastq) {
                reads_for_alignment = FASTP.out.reads.map{ meta, reads ->
                    read_files = reads.sort(false) { a,b -> a.getName().tokenize('.')[0] <=> b.getName().tokenize('.')[0] }.collate(2)
                    [ meta + [ size:read_files.size() ], read_files ]
                }.transpose()
            } else reads_for_alignment = FASTP.out.reads

            versions = versions.mix(FASTP.out.versions)

        } else {
            reads_for_alignment = reads_for_fastp
        }

        // STEP 1: MAPPING READS TO REFERENCE GENOME
        // reads will be sorted
        reads_for_alignment = reads_for_alignment.map{ meta, reads ->
            // Update meta.id to meta.sample no multiple lanes or splitted fastqs
            if (meta.size * meta.num_lanes == 1) [ meta + [ id:meta.sample ], reads ]
            else [ meta, reads ]
        }

        sort_bam = true
        FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON(reads_for_alignment, index_alignement, sort_bam, fasta, fasta_fai)

        // Grouping the bams from the same samples not to stall the workflow
        bam_mapped = FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON.out.bam.map{ meta, bam ->

            // Update meta.id to be meta.sample, ditching sample-lane that is not needed anymore
            // Update meta.data_type
            // Remove no longer necessary fields:
            //   read_group: Now in the BAM header
            //    num_lanes: only needed for mapping
            //         size: only needed for mapping

            // Use groupKey to make sure that the correct group can advance as soon as it is complete
            // and not stall the workflow until all reads from all channels are mapped
            [ groupKey( meta - meta.subMap('num_lanes', 'read_group', 'size') + [ data_type:'bam', id:meta.sample ], (meta.num_lanes ?: 1) * (meta.size ?: 1)), bam ]
        }.groupTuple()

        bai_mapped = FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON.out.bai.map{ meta, bai ->
            [ groupKey( meta - meta.subMap('num_lanes', 'read_group', 'size') + [ data_type:'bai', id:meta.sample ], (meta.num_lanes ?: 1) * (meta.size ?: 1)), bai ]
        }.groupTuple()

        // gatk4 markduplicates can handle multiple bams as input, so no need to merge/index here
        // Except if and only if save_mapped or (skipping markduplicates and sentieon-dedup)
        if (
            args.save_mapped ||
            (
                (args.skip_tools && args.skip_tools.split(',').contains('markduplicates')) &&
                !(args.tools && args.tools.split(',').contains('sentieon_dedup'))
            )
        ) {
            // bams are merged (when multiple lanes from the same sample), indexed and then converted to cram
            BAM_MERGE_INDEX_SAMTOOLS(bam_mapped)

            BAM_TO_CRAM_MAPPING(BAM_MERGE_INDEX_SAMTOOLS.out.bam_bai, fasta, fasta_fai)
            // Create CSV to restart from this step
            args.save_output_as_bam ? CHANNEL_ALIGN_CREATE_CSV(BAM_MERGE_INDEX_SAMTOOLS.out.bam_bai) : CHANNEL_ALIGN_CREATE_CSV(BAM_TO_CRAM_MAPPING.out.alignment_index)

            // Gather used softwares versions
            versions = versions.mix(BAM_MERGE_INDEX_SAMTOOLS.out.versions)
            versions = versions.mix(BAM_TO_CRAM_MAPPING.out.versions)
        }

        // Gather used softwares versions
        versions = versions.mix(CONVERT_FASTQ_INPUT.out.versions)
        versions = versions.mix(FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON.out.versions)
    }

    if (args.step in ['mapping', 'markduplicates']) {

        // ch_cram_no_markduplicates_restart = Channel.empty()
        cram_markduplicates_no_spark = Channel.empty()
        cram_sentieon_dedup          = Channel.empty()
        cram_markduplicates_spark    = Channel.empty()

        // STEP 2: markduplicates (+QC) + convert to CRAM

        // ch_bam_for_markduplicates will contain bam mapped with FASTQ_ALIGN_BWAMEM_MEM2_DRAGMAP_SENTIEON when step is mapping
        // Or bams that are specified in the samplesheet.csv when step is prepare_recalibration
        cram_for_markduplicates = args.step == 'mapping' ? bam_mapped : input_sample.map{ meta, input, index -> [ meta, input ] }
        // if no MD is done, then run QC on mapped & converted CRAM files
        // or the input BAM (+converted) or CRAM files
        cram_skip_markduplicates = Channel.empty()

        // Should it be possible to restart from converted crams?
        // For now, conversion from bam to cram is only done when skipping markduplicates

        if (
            args.skip_tools &&
            args.skip_tools.split(',').contains('markduplicates') &&
            !(args.tools && args.tools.split(',').contains('sentieon_dedup'))
        ) {
            if (args.step == 'mapping') {
                cram_skip_markduplicates = BAM_TO_CRAM_MAPPING.out.alignment_index
            } else {
                input_markduplicates_convert = input_sample.branch{
                    bam:  it[0].data_type == "bam"
                    cram: it[0].data_type == "cram"
                }

                // Convert any input BAMs to CRAM
                BAM_TO_CRAM(input_markduplicates_convert.bam, fasta, fasta_fai)
                versions = versions.mix(BAM_TO_CRAM.out.versions)

                cram_skip_markduplicates = Channel.empty().mix(input_markduplicates_convert.cram, BAM_TO_CRAM.out.alignment_index)
            }

            CRAM_QC_NO_MD(cram_skip_markduplicates, fasta, intervals_for_preprocessing)

            // Gather QC reports
            reports = reports.mix(CRAM_QC_NO_MD.out.reports.collect{ meta, report -> report })

            // Gather used softwares versions
            versions = versions.mix(CRAM_QC_NO_MD.out.versions)
        } else if (args.use_gatk_spark && args.use_gatk_spark.contains('markduplicates')) {
            BAM_MARKDUPLICATES_SPARK(
                cram_for_markduplicates,
                dict.map{ meta, dict -> [ dict ] },
                fasta,
                fasta_fai,
                intervals_for_preprocessing)
            cram_markduplicates_spark = BAM_MARKDUPLICATES_SPARK.out.cram

            // Gather QC reports
            reports = reports.mix(BAM_MARKDUPLICATES_SPARK.out.reports.collect{ meta, report -> report })

            // Gather used softwares versions
            versions = versions.mix(BAM_MARKDUPLICATES_SPARK.out.versions)
        } else if (args.tools && args.tools.split(',').contains('sentieon_dedup')) {
            crai_for_markduplicates = args.step == 'mapping' ? bai_mapped : input_sample.map{ meta, input, index -> [ meta, index ] }
            BAM_SENTIEON_DEDUP(
                cram_for_markduplicates,
                crai_for_markduplicates,
                fasta,
                fasta_fai,
                intervals_for_preprocessing)

            cram_sentieon_dedup = BAM_SENTIEON_DEDUP.out.cram

            // Gather QC reports
            reports = reports.mix(BAM_SENTIEON_DEDUP.out.reports.collect{ meta, report -> report })

            // Gather used softwares versions
            versions = versions.mix(BAM_SENTIEON_DEDUP.out.versions)
        } else {
            BAM_MARKDUPLICATES(
                cram_for_markduplicates,
                fasta,
                fasta_fai,
                intervals_for_preprocessing)

            cram_markduplicates_no_spark = BAM_MARKDUPLICATES.out.cram

            // Gather QC reports
            reports = reports.mix(BAM_MARKDUPLICATES.out.reports.collect{ meta, report -> report })

            // Gather used softwares versions
            versions = versions.mix(BAM_MARKDUPLICATES.out.versions)
        }

        // ch_md_cram_for_restart contains either:
        // - crams from markduplicates
        // - crams from sentieon_dedup
        // - crams from markduplicates_spark
        // - crams from input step markduplicates --> from the converted ones only?
        ch_md_cram_for_restart = Channel.empty().mix(cram_markduplicates_no_spark, cram_markduplicates_spark, cram_sentieon_dedup)
            // Make sure correct data types are carried through
            .map{ meta, cram, crai -> [ meta + [data_type: "cram"], cram, crai ] }

        // If args.save_output_as_bam, then convert CRAM files to BAM
        CRAM_TO_BAM(ch_md_cram_for_restart, fasta, fasta_fai)
        versions = versions.mix(CRAM_TO_BAM.out.versions)

        // CSV should be written for the file actually out, either CRAM or BAM
        // Create CSV to restart from this step
        csv_subfolder = (args.tools && args.tools.split(',').contains('sentieon_dedup')) ? 'sentieon_dedup' : 'markduplicates'

        args.save_output_as_bam ? CHANNEL_MARKDUPLICATES_CREATE_CSV(CRAM_TO_BAM.out.alignment_index, csv_subfolder, args.outdir, args.save_output_as_bam) : CHANNEL_MARKDUPLICATES_CREATE_CSV(ch_md_cram_for_restart, csv_subfolder, args.outdir, args.save_output_as_bam)
    }

    if (args.step in ['mapping', 'markduplicates', 'prepare_recalibration']) {

        // Run if starting from step "prepare_recalibration"
        if (args.step == 'prepare_recalibration') {

            // Support if starting from BAM or CRAM files
            input_prepare_recal_convert = input_sample.branch{
                bam: it[0].data_type == "bam"
                cram: it[0].data_type == "cram"
            }

            // BAM files first must be converted to CRAM files since from this step on we base everything on CRAM format
            BAM_TO_CRAM(input_prepare_recal_convert.bam, fasta, fasta_fai)
            versions = versions.mix(BAM_TO_CRAM.out.versions)

            ch_cram_from_bam = BAM_TO_CRAM.out.alignment_index
                // Make sure correct data types are carried through
                .map{ meta, cram, crai -> [ meta + [data_type: "cram"], cram, crai ] }

            ch_cram_for_bam_baserecalibrator = Channel.empty().mix(ch_cram_from_bam, input_prepare_recal_convert.cram)
            ch_md_cram_for_restart = ch_cram_from_bam

        } else {

            // ch_cram_for_bam_baserecalibrator contains either:
            // - crams from markduplicates
            // - crams from markduplicates_spark
            // - crams converted from bam mapped when skipping markduplicates
            // - input cram files, when start from step markduplicates
            ch_cram_for_bam_baserecalibrator = Channel.empty().mix(ch_md_cram_for_restart, cram_skip_markduplicates )
                // Make sure correct data types are carried through
                .map{ meta, cram, crai -> [ meta + [data_type: "cram"], cram, crai ] }

        }

        // STEP 3: Create recalibration tables
        if (!(args.skip_tools && args.skip_tools.split(',').contains('baserecalibrator'))) {

            ch_table_bqsr_no_spark = Channel.empty()
            ch_table_bqsr_spark    = Channel.empty()

            if (args.use_gatk_spark && args.use_gatk_spark.contains('baserecalibrator')) {
            BAM_BASERECALIBRATOR_SPARK(
                ch_cram_for_bam_baserecalibrator,
                dict,
                fasta,
                fasta_fai,
                intervals_and_num_intervals,
                known_sites_indels,
                known_sites_indels_tbi)

                ch_table_bqsr_spark = BAM_BASERECALIBRATOR_SPARK.out.table_bqsr

                // Gather used softwares versions
                versions = versions.mix(BAM_BASERECALIBRATOR_SPARK.out.versions)
            } else {

            BAM_BASERECALIBRATOR(
                ch_cram_for_bam_baserecalibrator,
                dict,
                fasta,
                fasta_fai,
                intervals_and_num_intervals,
                known_sites_indels,
                known_sites_indels_tbi)

                ch_table_bqsr_no_spark = BAM_BASERECALIBRATOR.out.table_bqsr

                // Gather used softwares versions
                versions = versions.mix(BAM_BASERECALIBRATOR.out.versions)
            }

            // ch_table_bqsr contains either:
            // - bqsr table from baserecalibrator
            // - bqsr table from baserecalibrator_spark
            ch_table_bqsr = Channel.empty().mix(
                ch_table_bqsr_no_spark,
                ch_table_bqsr_spark)

            reports = reports.mix(ch_table_bqsr.collect{ meta, table -> table })

            cram_applybqsr = ch_cram_for_bam_baserecalibrator.join(ch_table_bqsr, failOnDuplicate: true, failOnMismatch: true)

            // Create CSV to restart from this step
            CHANNEL_BASERECALIBRATOR_CREATE_CSV(ch_md_cram_for_restart.join(ch_table_bqsr, failOnDuplicate: true), args.tools, args.skip_tools, args.save_output_as_bam, args.outdir)
        }
    }

    // STEP 4: RECALIBRATING
    if (args.step in ['mapping', 'markduplicates', 'prepare_recalibration', 'recalibrate']) {

        // Run if starting from step "prepare_recalibration"
        if (args.step == 'recalibrate') {

            // Support if starting from BAM or CRAM files
            input_recal_convert = input_sample.branch{
                bam:  it[0].data_type == "bam"
                cram: it[0].data_type == "cram"
            }

            // If BAM file, split up table and mapped file to convert BAM to CRAM
            input_only_table = input_recal_convert.bam.map{ meta, bam, bai, table -> [ meta, table ] }
            input_only_bam   = input_recal_convert.bam.map{ meta, bam, bai, table -> [ meta, bam, bai ] }

            // BAM files first must be converted to CRAM files since from this step on we base everything on CRAM format
            BAM_TO_CRAM(input_only_bam, fasta, fasta_fai)
            versions = versions.mix(BAM_TO_CRAM.out.versions)

            cram_applybqsr = Channel.empty().mix(
                BAM_TO_CRAM.out.alignment_index.join(input_only_table, failOnDuplicate: true, failOnMismatch: true),
                input_recal_convert.cram)
                // Join together converted cram with input tables
                .map{ meta, cram, crai, table -> [ meta + [data_type: "cram"], cram, crai, table ]}
        }

        if (!(args.skip_tools && args.skip_tools.split(',').contains('baserecalibrator'))) {
            cram_variant_calling_no_spark = Channel.empty()
            cram_variant_calling_spark    = Channel.empty()

            if (args.use_gatk_spark && args.use_gatk_spark.contains('baserecalibrator')) {

                BAM_APPLYBQSR_SPARK(
                    cram_applybqsr,
                    dict,
                    fasta,
                    fasta_fai,
                    intervals_and_num_intervals)

                cram_variant_calling_spark = BAM_APPLYBQSR_SPARK.out.cram

                // Gather used softwares versions
                versions = versions.mix(BAM_APPLYBQSR_SPARK.out.versions)

            } else {

                BAM_APPLYBQSR(
                    cram_applybqsr,
                    dict,
                    fasta,
                    fasta_fai,
                    intervals_and_num_intervals)

                cram_variant_calling_no_spark = BAM_APPLYBQSR.out.cram

                // Gather used softwares versions
                versions = versions.mix(BAM_APPLYBQSR.out.versions)
            }

            cram_variant_calling = Channel.empty().mix(
                cram_variant_calling_no_spark,
                cram_variant_calling_spark)

            CRAM_QC_RECAL(
                cram_variant_calling,
                fasta,
                intervals_for_preprocessing)

            // Gather QC reports
            reports = reports.mix(CRAM_QC_RECAL.out.reports.collect{ meta, report -> report })

            // Gather used softwares versions
            versions = versions.mix(CRAM_QC_RECAL.out.versions)

            // If args.save_output_as_bam, then convert CRAM files to BAM
            CRAM_TO_BAM_RECAL(cram_variant_calling, fasta, fasta_fai)
            versions = versions.mix(CRAM_TO_BAM_RECAL.out.versions)

            // CSV should be written for the file actually out out, either CRAM or BAM
            csv_recalibration = Channel.empty()
            csv_recalibration = args.save_output_as_bam ?  CRAM_TO_BAM_RECAL.out.alignment_index : cram_variant_calling

            // Create CSV to restart from this step
            CHANNEL_APPLYBQSR_CREATE_CSV(csv_recalibration)

        } else if (args.step == 'recalibrate') {
            // cram_variant_calling contains either:
            // - input bams converted to crams, if started from step recal + skip BQSR
            // - input crams if started from step recal + skip BQSR
            cram_variant_calling = Channel.empty().mix(
                BAM_TO_CRAM.out.alignment_index,
                input_recal_convert.cram.map{ meta, cram, crai, table -> [ meta, cram, crai ] })
        } else {
            // cram_variant_calling contains either:
            // - crams from markduplicates = ch_cram_for_bam_baserecalibrator if skip BQSR but not started from step recalibration
            cram_variant_calling = Channel.empty().mix(ch_cram_for_bam_baserecalibrator)
        }
    }

    if (args.step == 'variant_calling') {

        input_variant_calling_convert = input_sample.branch{
            bam:  it[0].data_type == "bam"
            cram: it[0].data_type == "cram"
        }

        // BAM files first must be converted to CRAM files since from this step on we base everything on CRAM format
        BAM_TO_CRAM(input_variant_calling_convert.bam, fasta, fasta_fai)
        versions = versions.mix(BAM_TO_CRAM.out.versions)

        cram_variant_calling = Channel.empty().mix(BAM_TO_CRAM.out.alignment_index, input_variant_calling_convert.cram)

    }

    if (args.tools) {

        if (args.step == 'annotate') cram_variant_calling = Channel.empty()

        //
        // Logic to separate germline samples, tumor samples with no matched normal, and combine tumor-normal pairs
        //
        cram_variant_calling_status = cram_variant_calling.branch{
            normal: it[0].status == 0
            tumor:  it[0].status == 1
        }

        // All Germline samples
        cram_variant_calling_normal_to_cross = cram_variant_calling_status.normal.map{ meta, cram, crai -> [ meta.patient, meta, cram, crai ] }

        // All tumor samples
        cram_variant_calling_pair_to_cross = cram_variant_calling_status.tumor.map{ meta, cram, crai -> [ meta.patient, meta, cram, crai ] }

        // Tumor only samples
        // 1. Group together all tumor samples by patient ID [ patient1, [ meta1, meta2 ], [ cram1, crai1, cram2, crai2 ] ]

        // Downside: this only works by waiting for all tumor samples to finish preprocessing, since no group size is provided
        cram_variant_calling_tumor_grouped = cram_variant_calling_pair_to_cross.groupTuple()

        // 2. Join with normal samples, in each channel there is one key per patient now. Patients without matched normal end up with: [ patient1, [ meta1, meta2 ], [ cram1, crai1, cram2, crai2 ], null ]
        cram_variant_calling_tumor_joined = cram_variant_calling_tumor_grouped.join(cram_variant_calling_normal_to_cross, failOnDuplicate: true, remainder: true)

        // 3. Filter out entries with last entry null
        cram_variant_calling_tumor_filtered = cram_variant_calling_tumor_joined.filter{ it ->  !(it.last()) }

        // 4. Transpose [ patient1, [ meta1, meta2 ], [ cram1, crai1, cram2, crai2 ] ] back to [ patient1, meta1, [ cram1, crai1 ], null ] [ patient1, meta2, [ cram2, crai2 ], null ]
        // and remove patient ID field & null value for further processing [ meta1, [ cram1, crai1 ] ] [ meta2, [ cram2, crai2 ] ]
        cram_variant_calling_tumor_only = cram_variant_calling_tumor_filtered.transpose().map{ it -> [it[1], it[2], it[3]] }

        if (args.only_paired_variant_calling) {
            // Normal only samples

            // 1. Join with tumor samples, in each channel there is one key per patient now. Patients without matched tumor end up with: [ patient1, [ meta1 ], [ cram1, crai1 ], null ] as there is only one matched normal possible
            cram_variant_calling_normal_joined = cram_variant_calling_normal_to_cross.join(cram_variant_calling_tumor_grouped, failOnDuplicate: true, remainder: true)

            // 2. Filter out entries with last entry null
            cram_variant_calling_normal_filtered = cram_variant_calling_normal_joined.filter{ it ->  !(it.last()) }

            // 3. Remove patient ID field & null value for further processing [ meta1, [ cram1, crai1 ] ] [ meta2, [ cram2, crai2 ] ] (no transposing needed since only one normal per patient ID)
            cram_variant_calling_status_normal = cram_variant_calling_normal_filtered.map{ it -> [it[1], it[2], it[3]] }

        } else {
            cram_variant_calling_status_normal = cram_variant_calling_status.normal
        }

        // Tumor - normal pairs
        // Use cross to combine normal with all tumor samples, i.e. multi tumor samples from recurrences
        cram_variant_calling_pair = cram_variant_calling_normal_to_cross.cross(cram_variant_calling_pair_to_cross)
            .map { normal, tumor ->
                def meta = [:]

                meta.id         = "${tumor[1].sample}_vs_${normal[1].sample}".toString()
                meta.normal_id  = normal[1].sample
                meta.patient    = normal[0]
                meta.sex        = normal[1].sex
                meta.tumor_id   = tumor[1].sample

                [ meta, normal[2], normal[3], tumor[2], tumor[3] ]
            }

        // GERMLINE VARIANT CALLING
        BAM_VARIANT_CALLING_GERMLINE_ALL(
            args.tools,
            args.skip_tools,
            cram_variant_calling_status_normal,
            [ [ id:'bwa' ], [] ], // bwa_index for tiddit; not used here
            dbsnp,
            dbsnp_tbi,
            dbsnp_vqsr,
            dict,
            fasta,
            fasta_fai,
            intervals_and_num_intervals,
            intervals_bed_combined, // [] if no_intervals, else interval_bed_combined.bed,
            intervals_bed_gz_tbi_combined, // [] if no_intervals, else interval_bed_combined_gz, interval_bed_combined_gz_tbi
            PREPARE_INTERVALS.out.intervals_bed_combined, // no_intervals.bed if no intervals, else interval_bed_combined.bed; Channel operations possible
            intervals_bed_gz_tbi_and_num_intervals,
            known_indels_vqsr,
            known_sites_indels,
            known_sites_indels_tbi,
            known_sites_snps,
            known_sites_snps_tbi,
            known_snps_vqsr,
            args.joint_germline,
            args.skip_tools && args.skip_tools.split(',').contains('haplotypecaller_filter'), // true if filtering should be skipped
            args.sentieon_haplotyper_emit_mode)

        // TUMOR ONLY VARIANT CALLING
        BAM_VARIANT_CALLING_TUMOR_ONLY_ALL(
            args.tools,
            cram_variant_calling_tumor_only,
            [ [ id:'bwa' ], [] ], // bwa_index for tiddit; not used here
            cf_chrom_len,
            chr_files,
            cnvkit_reference,
            dbsnp,
            dbsnp_tbi,
            dict,
            fasta,
            fasta_fai,
            germline_resource,
            germline_resource_tbi,
            intervals_and_num_intervals,
            intervals_bed_gz_tbi_and_num_intervals,
            intervals_bed_combined,
            intervals_bed_gz_tbi_combined, // [] if no_intervals, else interval_bed_combined_gz, interval_bed_combined_gz_tbi
            mappability,
            pon,
            pon_tbi,
            args.joint_mutect2,
            args.wes
        )

        // PAIR VARIANT CALLING
        BAM_VARIANT_CALLING_SOMATIC_ALL(
            args.tools,
            cram_variant_calling_pair,
            [ [ id:'bwa' ], [] ], // bwa_index for tiddit; not used here
            cf_chrom_len,
            chr_files,
            dbsnp,
            dbsnp_tbi,
            dict,
            fasta,
            fasta_fai,
            germline_resource,
            germline_resource_tbi,
            intervals_and_num_intervals,
            intervals_bed_gz_tbi_and_num_intervals,
            intervals_bed_combined,
            intervals_bed_gz_tbi_combined, // [] if no_intervals, else interval_bed_combined_gz, interval_bed_combined_gz_tbi
            mappability,
            msisensorpro_scan,
            pon,
            pon_tbi,
            allele_files,
            loci_files,
            gc_file,
            rt_file,
            args.joint_mutect2,
            args.wes
        )

        // POST VARIANTCALLING
        POST_VARIANTCALLING(BAM_VARIANT_CALLING_GERMLINE_ALL.out.vcf_all,
                            args.concatenate_vcfs)

        // Gather vcf files for annotation and QC
        vcf_to_annotate = Channel.empty()
        vcf_to_annotate = vcf_to_annotate.mix(BAM_VARIANT_CALLING_GERMLINE_ALL.out.vcf_deepvariant)
        vcf_to_annotate = vcf_to_annotate.mix(BAM_VARIANT_CALLING_GERMLINE_ALL.out.vcf_freebayes)
        vcf_to_annotate = vcf_to_annotate.mix(BAM_VARIANT_CALLING_GERMLINE_ALL.out.vcf_haplotypecaller)
        vcf_to_annotate = vcf_to_annotate.mix(BAM_VARIANT_CALLING_GERMLINE_ALL.out.vcf_manta)
        vcf_to_annotate = vcf_to_annotate.mix(BAM_VARIANT_CALLING_GERMLINE_ALL.out.vcf_sentieon_haplotyper)
        vcf_to_annotate = vcf_to_annotate.mix(BAM_VARIANT_CALLING_GERMLINE_ALL.out.vcf_strelka)
        vcf_to_annotate = vcf_to_annotate.mix(BAM_VARIANT_CALLING_GERMLINE_ALL.out.vcf_tiddit)
        vcf_to_annotate = vcf_to_annotate.mix(BAM_VARIANT_CALLING_GERMLINE_ALL.out.vcf_mpileup)
        vcf_to_annotate = vcf_to_annotate.mix(BAM_VARIANT_CALLING_TUMOR_ONLY_ALL.out.vcf_all)
        vcf_to_annotate = vcf_to_annotate.mix(BAM_VARIANT_CALLING_SOMATIC_ALL.out.vcf_all)

        // QC
        VCF_QC_BCFTOOLS_VCFTOOLS(vcf_to_annotate, intervals_bed_combined)

        reports = reports.mix(VCF_QC_BCFTOOLS_VCFTOOLS.out.bcftools_stats.collect{ meta, stats -> stats })
        reports = reports.mix(VCF_QC_BCFTOOLS_VCFTOOLS.out.vcftools_tstv_counts.collect{ meta, counts -> counts })
        reports = reports.mix(VCF_QC_BCFTOOLS_VCFTOOLS.out.vcftools_tstv_qual.collect{ meta, qual -> qual })
        reports = reports.mix(VCF_QC_BCFTOOLS_VCFTOOLS.out.vcftools_filter_summary.collect{ meta, summary -> summary })

        CHANNEL_VARIANT_CALLING_CREATE_CSV(vcf_to_annotate)

        // Gather used variant calling softwares versions
        versions = versions.mix(BAM_VARIANT_CALLING_GERMLINE_ALL.out.versions)
        versions = versions.mix(BAM_VARIANT_CALLING_SOMATIC_ALL.out.versions)
        versions = versions.mix(BAM_VARIANT_CALLING_TUMOR_ONLY_ALL.out.versions)
        versions = versions.mix(POST_VARIANTCALLING.out.versions)
        versions = versions.mix(VCF_QC_BCFTOOLS_VCFTOOLS.out.versions)

        // ANNOTATE
        if (args.step == 'annotate') vcf_to_annotate = input_sample

        if (args.tools.split(',').contains('merge') || args.tools.split(',').contains('snpeff') || args.tools.split(',').contains('vep')) {

            vep_fasta = (args.vep_include_fasta) ? fasta.map{ fasta -> [ [ id:fasta.baseName ], fasta ] } : [[id: 'null'], []]

            VCF_ANNOTATE_ALL(
                vcf_to_annotate.map{meta, vcf -> [ meta + [ file_name: vcf.baseName ], vcf ] },
                vep_fasta,
                args.tools,
                args.snpeff_genome ? "${args.snpeff_genome}.${args.snpeff_db}" : "${args.genome}.${args.snpeff_db}",
                snpeff_cache,
                vep_genome,
                vep_species,
                vep_cache_version,
                vep_cache,
                vep_extra_files)

            // Gather used softwares versions
            versions = versions.mix(VCF_ANNOTATE_ALL.out.versions)
            reports = reports.mix(VCF_ANNOTATE_ALL.out.reports)
        }
    }

    version_yaml = Channel.empty()
    if (!(args.skip_tools && args.skip_tools.split(',').contains('versions'))) {
        CUSTOM_DUMPSOFTWAREVERSIONS(versions.unique().collectFile(name: 'collated_versions.yml'))
        version_yaml = CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect()
    }

    if (!(args.skip_tools && args.skip_tools.split(',').contains('multiqc'))) {
        multiqc_files = Channel.empty()
        multiqc_files = multiqc_files.mix(version_yaml)
        multiqc_files = multiqc_files.mix(reports.collect().ifEmpty([]))

        MULTIQC(multiqc_files.collect(), ch_multiqc_config.collect().ifEmpty([]), ch_multiqc_custom_config.collect().ifEmpty([]), ch_multiqc_logo.collect().ifEmpty([]))

        multiqc_report = MULTIQC.out.report.toList()
        versions = versions.mix(MULTIQC.out.versions)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Parse first line of a FASTQ file, return the flowcell id and lane number.
def flowcellLaneFromFastq(path) {
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven fields)
    // or
    // FLOWCELLID:LANE:xx:... (five fields)
    def line
    path.withInputStream {
        InputStream gzipStream = new java.util.zip.GZIPInputStream(it)
        Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
        BufferedReader buffered = new BufferedReader(decoder)
        line = buffered.readLine()
    }
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(':')
    String fcid

    if (fields.size() >= 7) {
        // CASAVA 1.8+ format, from  https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/FileFormat_FASTQ-files_swBS.htm
        // "@<instrument>:<run number>:<flowcell ID>:<lane>:<tile>:<x-pos>:<y-pos>:<UMI> <read>:<is filtered>:<control number>:<index>"
        fcid = fields[2]
    } else if (fields.size() == 5) {
        fcid = fields[0]
    }
    return fcid
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
