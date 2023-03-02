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
import java.util.*;

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

    public void determineStarAlleles() throws IOException {

        for (String key : genes.keySet()) {
            File logFile = new File(STAR_ALLELE_OUTPUT_DIR + key + "_star_alleles_log.txt");
            if (!logFile.exists()) logFile.createNewFile();
            FileWriter fileWriter = new FileWriter(logFile.getAbsoluteFile());
            BufferedWriter logWriter = new BufferedWriter(fileWriter);

            // Keys are the non wild-type snps, values are the samples+hap that matches the set of snps
            Map<Set<Snp>, List<String>> nanUniqueCombinations = new HashMap<>();

            PgxGene pgxGene = genes.get(key);

            int duplicateCount = 0;
            int naCount = 0;

            for (PgxSample sample : samples.values()) {
                sample.getGenes().put(pgxGene.getName(), new PgxGene(pgxGene.getName()));

                String starAllele0 = null;
                String starAllele1 = null;

                boolean hasDuplicate = false;

                Set<Snp> snpsOnHaplotype0 = sample.getHaplotype0(pgxGene.getVariants());
                Set<Snp> snpsOnHaplotype1 = sample.getHaplotype1(pgxGene.getVariants());

                for (PgxHaplotype pgxHaplotype : pgxGene.getPgxHaplotypes().values()) {
                    Set<Snp> snpsOnPgxHaploType = new HashSet<>(pgxHaplotype.getSnpsWithWildType().values());

                    if (snpsOnHaplotype0.equals(snpsOnPgxHaploType)) {
                        if (starAllele0 != null) {
                            logWriter.write(sample.getId() + ": duplicate haplo 0 " + starAllele0 + " " + pgxHaplotype.getName() + "\n");
                            hasDuplicate = true;
                        }
                        starAllele0 = pgxHaplotype.getName();
                        sample.getGenes().get(pgxGene.getName()).setAllele0(starAllele0);
                    }

                    if (snpsOnHaplotype1.equals(snpsOnPgxHaploType)) {
                        if (starAllele1 != null) {
                            logWriter.write(sample.getId() + ": duplicate haplo 1 " + starAllele1 + " " + pgxHaplotype.getName() + "\n");
                            hasDuplicate = true;
                        }
                        starAllele1 = pgxHaplotype.getName();
                        sample.getGenes().get(pgxGene.getName()).setAllele1(starAllele1);
                    }

                    if (starAllele0 != null & starAllele1 != null) {
                        //break; // this shouldn't happen when there are no similar entries in the translation tables
                    }
                }

                if (starAllele0 == null) {
                    sample.getGenes().get(pgxGene.getName()).setAllele0("NA");
                    //sample.getGenes().get(pgxGene.getName()).setAllele0(pgxGene.getWildType().getName());
                    determineNaHaplotypes(nanUniqueCombinations, snpsOnHaplotype0, sample, 0);
                    naCount += 1;
                }
                if (starAllele1 == null) {
                    sample.getGenes().get(pgxGene.getName()).setAllele1("NA");
                    //sample.getGenes().get(pgxGene.getName()).setAllele1(pgxGene.getWildType().getName());
                    determineNaHaplotypes(nanUniqueCombinations, snpsOnHaplotype1, sample, 1);
                    naCount += 1;
                }
                if (hasDuplicate) duplicateCount += 1;
            }
            System.out.println(pgxGene.getName() + ": " + duplicateCount + " duplicates");
            System.out.println(pgxGene.getName() + ": " + naCount + " nans");

            logNaHaplotype(logWriter, nanUniqueCombinations);

            logWriter.close();
        }
    }

    private void determineNaHaplotypes(Map<Set<Snp>, List<String>> nanUniqueCombinations, Set<Snp> snps, PgxSample sample, int haplotype) {
        Set<Snp> nonWildTypeSnps = new HashSet<>();
        for (Snp snp : snps) {
            if (!snp.getVariantAllele().equals(snp.getReferenceAllele())) {
                nonWildTypeSnps.add(snp);
            }
        }
        String samplesHap = sample.getId() + "_hap_" + haplotype;
        if (nanUniqueCombinations.containsKey(nonWildTypeSnps)) {
            List<String> samples = nanUniqueCombinations.get(nonWildTypeSnps);
            samples.add(samplesHap);
        } else {
            List<String> samples = new ArrayList<>();
            samples.add(samplesHap);
            nanUniqueCombinations.put(nonWildTypeSnps, samples);
        }
    }

    private void logNaHaplotype(BufferedWriter logWriter, Map<Set<Snp>, List<String>> nanUniqueCombinations) throws IOException {
        TreeMap<Set<Snp>, List<String>> sorter = new TreeMap(new ListSizeComparator(nanUniqueCombinations));
        sorter.putAll(nanUniqueCombinations);

        for (Set<Snp> snps : sorter.keySet()) {
            logWriter.write("\nNo translation for SNPs:\n");
            for (Snp snp : snps) {
                logWriter.write(snp.getId() + ", reference: " + snp.getReferenceAllele() + ", variant: " + snp.getVariantAllele() + "\n");
            }
            List<String> sampleHaplotypes = nanUniqueCombinations.get(snps);
            logWriter.write(sampleHaplotypes.size() + " affected samples:\n");
            for (String sampleHaplotype : sampleHaplotypes) {
                logWriter.write(sampleHaplotype + ", ");
            }
            logWriter.write("\n");

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
