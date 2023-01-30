package org.molgenis.asterix.pipeline;

import org.molgenis.asterix.config.ConfigConstants;
import org.molgenis.asterix.config.ConfigProvider;
import org.molgenis.asterix.io.HaplotypeReader;
import org.molgenis.asterix.io.SnpToHaploTableReader;
import org.molgenis.asterix.io.VcfHaplotypeReader;
import org.molgenis.asterix.model.PgxGene;
import org.molgenis.asterix.model.PgxHaplotype;
import org.molgenis.asterix.model.PgxSample;
import org.molgenis.asterix.model.Snp;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;

public class HaplotypeToStarAllele {

    private static String STAR_ALLELE_OUTPUT_DIR;
    private static String SNP_HAPLO_TABLE_DIR;
    private static HaplotypeReader HAPLOTYPE_READER;

    // Gene identifiers as keys
    private SortedMap<String, PgxGene> genes;
    // Sample identifiers as keys
    private Map<String, PgxSample> samples;

    public HaplotypeToStarAllele() throws IOException {
        this.loadConfig();

        SnpToHaploTableReader snpToHaploTableReader = new SnpToHaploTableReader(SNP_HAPLO_TABLE_DIR);
        snpToHaploTableReader.readSnpHaploTables();
        this.genes = snpToHaploTableReader.getGenes();

        HaplotypeReader haplotypeReader = HAPLOTYPE_READER;
        haplotypeReader.readHaplotypes(genes.values());
        this.samples = haplotypeReader.getSamples();
    }

    public void determineStarAlleles() {
        for (String key : genes.keySet()) {
            PgxGene pgxGene = genes.get(key);

            for (PgxSample sample : samples.values()) {

                sample.getGenes().put(pgxGene.getName(), new PgxGene(pgxGene.getName()));

                String starAllele0 = null;
                String starAllele1 = null;

                for (PgxHaplotype pgxHaplotype : pgxGene.getPgxHaplotypes().values()) {

                    Set<Snp> snpsOnPgxHaploType = new HashSet<>(pgxHaplotype.getSnps().values());

                    Set<Snp> snpsOnHaplotype0 = new HashSet<>();
                    Set<Snp> snpsOnHaplotype1 = new HashSet<>();
                    for (Snp snp : snpsOnPgxHaploType) {
                        snpsOnHaplotype0.add(sample.getHaplotype0().get(snp.getId()));
                        snpsOnHaplotype1.add(sample.getHaplotype1().get(snp.getId()));
                    }

                    if (snpsOnHaplotype0.equals(snpsOnPgxHaploType)) {
                        if (starAllele0 != null) {
                            System.out.println("Duplicate star allele: " + starAllele0 + " " + pgxHaplotype.getName());
                        }
                        starAllele0 = pgxHaplotype.getName();
                        sample.getGenes().get(pgxGene.getName()).setAllele0(starAllele0);
                    }

                    if (snpsOnHaplotype1.equals(snpsOnPgxHaploType)) {
                        if (starAllele1 != null) {
                            System.out.println("Duplicate star allele: " + starAllele1 + " " + pgxHaplotype.getName());
                        }
                        starAllele1 = pgxHaplotype.getName();
                        sample.getGenes().get(pgxGene.getName()).setAllele1(starAllele1);
                    }

                    if (starAllele0 != null & starAllele1 != null) {
                        break;
                    }
                }

                //TODO: Here we allow a mismatch and fallback to the wild type
//                if (starAllele0 == null) sample.getGenes().get(pgxGene.getName()).setAllele0(pgxGene.getWildType().getName());
//                if (starAllele1 == null) sample.getGenes().get(pgxGene.getName()).setAllele1(pgxGene.getWildType().getName());

                // For now don't allow fallback and set allele to NA
                if (starAllele0 == null) sample.getGenes().get(pgxGene.getName()).setAllele0("NA");
                if (starAllele1 == null) sample.getGenes().get(pgxGene.getName()).setAllele1("NA");

            }

        }

    }

    public void writeStarAlleles() throws IOException {
        System.out.println("Writing output files");

        new File(STAR_ALLELE_OUTPUT_DIR).mkdirs();

        for (PgxGene pgxGene : genes.values()) {

            File startAlleleOutputFile = new File(STAR_ALLELE_OUTPUT_DIR + pgxGene.getName() + "_star_alleles.txt");
            if (!startAlleleOutputFile.exists()) startAlleleOutputFile.createNewFile();

            FileWriter fileWriter = new FileWriter(startAlleleOutputFile.getAbsoluteFile());
            BufferedWriter bw = new BufferedWriter(fileWriter);

            for (PgxSample sample : samples.values()) {
                String line = sample.getId() + "\t" + sample.getGenes().get(pgxGene.getName()).getAllele0() + "\t" +
                        sample.getGenes().get(pgxGene.getName()).getAllele1() + "\n";
                bw.write(line);
            }

            bw.close();
        }
    }

    public SortedMap<String, PgxGene> getGenes() {
        return genes;
    }

    public Map<String, PgxSample> getSamples() {
        return samples;
    }

    /**
     * Loads the configuration
     */
    private void loadConfig() {
        //set config provider
        ConfigProvider configProvider = ConfigProvider.getInstance();
        //load dirs
        STAR_ALLELE_OUTPUT_DIR = configProvider.getConfigParam(ConfigConstants.STAR_ALLELE_OUTPUT_DIR);
        SNP_HAPLO_TABLE_DIR = configProvider.getConfigParam(ConfigConstants.SNP_HAPLO_TABLE_DIR);
        HAPLOTYPE_READER = new VcfHaplotypeReader(configProvider.getConfigParam(ConfigConstants.HAPLOTYPE_DIR));
    }

}
