package org.molgenis.asterix.config;

/**
 * @author OelenR
 *
 * class to store the constants in that are used for configuration
 */
public class ConfigConstants {

    //output dir staralleles of every individual
    public static final String STAR_ALLELE_OUTPUT_DIR = "STAR_ALLELE_OUTPUT_DIR";
    protected static final String OPTION_STAR_ALLELE_OUTPUT_DIR = "star_out";
    protected static final String DEFAULT_STAR_ALLELE_OUTPUT_DIR = "star_alleles/";

    //input dir conversion/mapping table snps to haplotypes
    public static final String SNP_HAPLO_TABLE_DIR = "SNP_HAPLO_TABLE_DIR";
    protected static final String OPTION_SNP_HAPLO_TABLE_DIR = "snp_haplo_table_dir";
    protected static final String DEFAULT_SNP_HAPLO_TABLE_DIR = "snp_to_haplo/";

    //genotype data, .haps file that needs to be converted
    public static final String HAPLOTYPE_DIR = "HAPLOTYPE_DIR";
    protected static final String OPTION_HAPLOTYPE_DIR = "haplo_type_dir";
    protected static final String DEFAULT_HAPLOTYPE_DIR = "ll_phased_active/";

    //input dir for converting haplotype to phenotype (star alleles to functions)
    public static final String HAPLO_PHENO_TABLE_DIR = "HAPLO_PHENO_TABLE_DIR";
    protected static final String OPTION_HAPLO_PHENO_TABLE_DIR = "haplo_pheno_table_dir";
    protected static final String DEFAULT_HAPLO_PHENO_TABLE_DIR = "haplo_to_pheno/";

    //output dir for converted haplotypes to predicted functions (per gene as file, with rows as persons)
    public static final String PREDICTED_PHENOTYPES_OUTPUT_DIR = "PREDICTED_PHENOTYPES_OUTPUT_DIR";
    protected static final String OPTION_PREDICTED_PHENOTYPES_OUTPUT_DIR = "pheno_out_dir";
    protected static final String DEFAULT_PREDICTED_PHENOTYPES_OUTPUT_DIR = "predicted_phenotypes/";

    //output file for sample matrix
    public static final String SAMPLE_MATRIX_OUT = "SAMPLE_MATRIX_OUT";
    protected static final String OPTION_SAMPLE_MATRIX_OUT = "sample_matrix_out";
    protected static final String DEFAULT_SAMPLE_MATRIX_OUT = "sample_matrix.csv";

    //whether to split the out file per person
    public static final String SPLIT_SAMPLES_PP = "SPLIT_SAMPLES_PP";
    protected static final String OPTION_SPLIT_SAMPLES_PP = "split_samples";
    protected static final String DEFAULT_SPLIT_SAMPLES_PP = "false";

    //properties loading via external .properties file
    protected static final String OPTION_PROPERTIES_FILE = "properties";

    //help string cli
    protected static final String OPTION_HELP = "help";

    //private constructor so Object cannot be created
    private ConfigConstants() {

    }
}
