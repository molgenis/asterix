package org.molgenis.asterix;

import org.molgenis.asterix.config.ConfigConstants;
import org.molgenis.asterix.config.ConfigProvider;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;

public class HaplotypeToStarAllele {

//    private static final String STAR_ALLELE_OUTPUT_DIR = "/groups/ll-geno/tmp04/ll-ov18-0434/pgx_profiles/star_alleles/";
//    private static final String SNP_HAPLO_TABLE_DIR = "/groups/ll-geno/tmp04/ll-ov18-0434/data/translation_tables/snp_to_haplo/";
//    private static final String HAPLOTYPE_DIR = "/groups/ll-geno/tmp04/ll-ov18-0434/data/phased_genotypes/";

    //output dir staralleles of every individual
    //private static final String STAR_ALLELE_OUTPUT_DIR = "/Users/harmbrugge/Documents/PGx/data/richtlijn/star_alleles/";
    //private static final String STAR_ALLELE_OUTPUT_DIR = "C:\\molgenis\\asterix_data\\star_alleles\\";
    private static String STAR_ALLELE_OUTPUT_DIR;
    //input dir conversion/mapping table snps to haplotypes
    //private static final String SNP_HAPLO_TABLE_DIR = "/Users/harmbrugge/Documents/PGx/data/richtlijn/snp_to_haplo/";
    //private static final String SNP_HAPLO_TABLE_DIR = "C:\\molgenis\\asterix_data\\snp_to_haplo\\";
    private static String SNP_HAPLO_TABLE_DIR;

    //genotype data, .haps file that needs to be converted
    //private static final String HAPLOTYPE_DIR = "/Users/harmbrugge/Documents/PGx/data/richtlijn/ll_phased_active/";
    //private static final String HAPLOTYPE_DIR = "C:\\molgenis\\asterix_data\\ll_phased_active\\";
    private static String HAPLOTYPE_DIR;

    private ConfigProvider configProvider = null;

    //gene as values, gene identifiers as keys
    private SortedMap<String, PgxGene> genes;

    //samples as values, sample identifiers as keys
    private Map<String, Sample> samples;



    public HaplotypeToStarAllele() throws IOException {
        //set config
        this.loadConfig();

        HaploTableReader haploTableReader = new HaploTableReader(SNP_HAPLO_TABLE_DIR);
        haploTableReader.readSnpHaploTables();
        this.genes = haploTableReader.getGenes();

        HaplotypeReader haplotypeReader = new HaplotypeReader(HAPLOTYPE_DIR);
        haplotypeReader.readHaplotypes();
        this.samples = haplotypeReader.getSamples();
    }

    public SortedMap<String, PgxGene> getGenes() {
        return genes;
    }

    public void determineStarAlleles() {

        for (String key : genes.keySet()) {
            PgxGene pgxGene = genes.get(key);

            for (Sample sample : samples.values()) {

                sample.getGenes().put(pgxGene.getName(), new PgxGene(pgxGene.getName()));

                String starAllele0 = null;
                String starAllele1 = null;

                // TODO: Is wild type always the first key?
                PgxHaplotype wildType = pgxGene.getPgxHaplotypes().get(pgxGene.getPgxHaplotypes().firstKey());

                for (PgxHaplotype pgxHaplotype : pgxGene.getPgxHaplotypes().values()) {

                    Set<Snp> snpsOnPgxHaploType = new HashSet<>(pgxHaplotype.getSnps().values());

                    Set<Snp> snpsOnHaplotype0 = new HashSet<>();
                    Set<Snp> snpsOnHaplotype1 = new HashSet<>();

                    for (Snp snp : snpsOnPgxHaploType) {
                        snpsOnHaplotype0.add(sample.getHaplotype0().get(snp.getId()));
                        snpsOnHaplotype1.add(sample.getHaplotype1().get(snp.getId()));
                    }

                    if (snpsOnHaplotype0.equals(snpsOnPgxHaploType)) {
                        starAllele0 = pgxHaplotype.getName();
                        sample.getGenes().get(pgxGene.getName()).setAllele0(starAllele0);
                    }

                    if (snpsOnHaplotype1.equals(snpsOnPgxHaploType)) {
                        starAllele1 = pgxHaplotype.getName();
                        sample.getGenes().get(pgxGene.getName()).setAllele1(starAllele1);
                    }

                    if (starAllele0 != null & starAllele1 != null) {
                        break;
                    }
                }

                //TODO: Here we allow a mismatch and fallback to the wildtype
                if (starAllele0 == null) {
                    sample.getGenes().get(pgxGene.getName()).setAllele0(wildType.getName());
                }

                if (starAllele1 == null) {
                    sample.getGenes().get(pgxGene.getName()).setAllele1(wildType.getName());
                }

            }

        }

    }

    public void writeStarAlleles() throws IOException {

        new File(STAR_ALLELE_OUTPUT_DIR).mkdirs();

        for (PgxGene pgxGene : genes.values()) {

            File startAlleleOutputFile = new File(STAR_ALLELE_OUTPUT_DIR + pgxGene.getName() + "_star_alleles.txt");
            if (!startAlleleOutputFile.exists()) startAlleleOutputFile.createNewFile();

            FileWriter fileWriter = new FileWriter(startAlleleOutputFile.getAbsoluteFile());
            BufferedWriter bw = new BufferedWriter(fileWriter);

            for (Sample sample : samples.values()) {
                String line = sample.getId() + "\t" + sample.getGenes().get(pgxGene.getName()).getAllele0() + "\t" + sample.getGenes().get(pgxGene.getName()).getAllele1() + "\n";
                bw.write(line);
            }

            bw.close();
        }

    }

    public Map<String, Sample> getSamples() {
        return samples;
    }

    /**
     * load the configuration
     */
    private void loadConfig(){
        //set config provider
        this.configProvider = ConfigProvider.getInstance();
        //load dirs
        STAR_ALLELE_OUTPUT_DIR = this.configProvider.getConfigParam(ConfigConstants.STAR_ALLELE_OUTPUT_DIR);
        SNP_HAPLO_TABLE_DIR = this.configProvider.getConfigParam(ConfigConstants.SNP_HAPLO_TABLE_DIR);
        HAPLOTYPE_DIR = this.configProvider.getConfigParam(ConfigConstants.HAPLOTYPE_DIR);
    }

}
