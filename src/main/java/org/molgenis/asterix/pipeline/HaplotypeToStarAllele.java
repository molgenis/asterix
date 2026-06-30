package org.molgenis.asterix.pipeline;

import org.molgenis.asterix.config.ConfigConstants;
import org.molgenis.asterix.config.ConfigProvider;
import org.molgenis.asterix.io.HaplotypeReader;
import org.molgenis.asterix.io.SnpToHaploTableReader;
import org.molgenis.asterix.io.VcfHaplotypeReader;
import org.molgenis.asterix.model.*;
import org.molgenis.asterix.utils.ListSizeComparator;
import org.molgenis.asterix.utils.LogSizeComparator;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

public class HaplotypeToStarAllele {

    private static final Double PROBABILITY_CUT_OFF = 0.97;
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

            // Keep track of the samples that do not match a translation.
            // Keys are the non wild-type snps, values are the samples that matches the set of snps without a translation
            Map<Set<Snp>, List<String>> nanUniqueCombinations = new HashMap<>();
            Map<Set<Snp>, List<SnpLog>> lowProbSnpsCombinations = new HashMap<>();
            // Keys are the reliably detected variant snps, values are the affected sample haplotypes
            Map<Set<Snp>, List<SnpLog>> detectedMutationsCombinations = new HashMap<>();

            PgxGene pgxGene = genes.get(key);

            int duplicateCount = 0;
            int naCount = 0;

            for (PgxSample sample : samples.values()) {

                sample.getGenes().put(pgxGene.getName(), new PgxGene(pgxGene));

                String starAllele0 = null;
                String starAllele1 = null;

                boolean hasDuplicate = false;

                Set<Snp> snpsOnHaplotype0 = sample.getHaplotype0(pgxGene.getVariants());
                Set<Snp> snpsOnHaplotype1 = sample.getHaplotype1(pgxGene.getVariants());

                Set<Snp> lowProbSnpsOn0 = new HashSet<>();
                Set<Snp> lowProbSnpsOn1 = new HashSet<>();

                Set<Snp> detectedMutationsOn0 = new HashSet<>();
                Set<Snp> detectedMutationsOn1 = new HashSet<>();

                for (PgxHaplotype pgxHaplotype : pgxGene.getPgxHaplotypes().values()) {
                    if (pgxHaplotype.getSnps().isEmpty()) continue;
                    Set<Snp> snpsOnPgxHaploType = new HashSet<>(pgxHaplotype.getSnpsWithWildType().values());

                    if (snpsOnHaplotype0.equals(snpsOnPgxHaploType)) {
                        lowProbSnpsOn0.addAll(checkProbabilities(snpsOnHaplotype0, pgxHaplotype));
                        detectedMutationsOn0.addAll(getDetectedMutations(snpsOnHaplotype0, pgxHaplotype));
                        if (starAllele0 != null) {
                            logWriter.write(sample.getId() + ": duplicate haplo 0 " + starAllele0 + " " +
                                    pgxHaplotype.getName() + "\n");
                            hasDuplicate = true;
                        }
                        starAllele0 = pgxHaplotype.getCorrected();
                        sample.getGenes().get(pgxGene.getName()).setDisclaimerSubject(pgxHaplotype.getDisclaimer());
                        sample.getGenes().get(pgxGene.getName()).setAllele0(starAllele0);
                    }

                    if (snpsOnHaplotype1.equals(snpsOnPgxHaploType)) {
                        lowProbSnpsOn1.addAll(checkProbabilities(snpsOnHaplotype1, pgxHaplotype));
                        detectedMutationsOn1.addAll(getDetectedMutations(snpsOnHaplotype1, pgxHaplotype));
                        if (starAllele1 != null) {
                            logWriter.write(sample.getId() + ": duplicate haplo 1 " + starAllele1 + " " +
                                    pgxHaplotype.getName() + "\n");
                            hasDuplicate = true;
                        }
                        starAllele1 = pgxHaplotype.getCorrected();
                        sample.getGenes().get(pgxGene.getName()).setDisclaimerSubject(pgxHaplotype.getDisclaimer());
                        sample.getGenes().get(pgxGene.getName()).setAllele1(starAllele1);
                    }

                    if (starAllele0 != null & starAllele1 != null) {
                        //break;
                    }
                }

                addSnpCombination(lowProbSnpsOn0, lowProbSnpsCombinations, sample, 0);
                addSnpCombination(lowProbSnpsOn1, lowProbSnpsCombinations, sample, 1);

                addSnpCombination(detectedMutationsOn0, detectedMutationsCombinations, sample, 0);
                addSnpCombination(detectedMutationsOn1, detectedMutationsCombinations, sample, 1);

                if (starAllele0 == null) {
                    sample.getGenes().get(pgxGene.getName()).setAllele0(pgxGene.getWildType().getName());
                    determineNaHaplotypes(nanUniqueCombinations, snpsOnHaplotype0, sample, 0);
                    naCount += 1;
                }
                if (starAllele1 == null) {
                    sample.getGenes().get(pgxGene.getName()).setAllele1(pgxGene.getWildType().getName());
                    determineNaHaplotypes(nanUniqueCombinations, snpsOnHaplotype1, sample, 1);
                    naCount += 1;
                }
                if (hasDuplicate) duplicateCount += 1;
            }
            System.out.println(pgxGene.getName() + ": " + duplicateCount + " duplicates");
            System.out.println(pgxGene.getName() + ": " + naCount + " nans");

            logNaHaplotype(logWriter, nanUniqueCombinations);
            logLowProbSnps(logWriter, lowProbSnpsCombinations);
            logDetectedMutations(logWriter, detectedMutationsCombinations);

            if (pgxGene.getName().equals("CYP2D6")) {
                new Cyp2d6Caller(samples, pgxGene);
            }

            logWriter.close();
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

    private Set<Snp> checkProbabilities(Set<Snp> snpsOnHaplotype, PgxHaplotype pgxHaplotype) {
        Set<Snp> lowProbSnps = new HashSet<>();
        for (Snp snp : snpsOnHaplotype) {
            snp.setOriginalCall(snp.getVariantAllele());
            if (pgxHaplotype.getSnps().containsKey(snp.getId()) &
                    snp.getProbability() < PROBABILITY_CUT_OFF) {
                lowProbSnps.add(snp);
                if (!snp.getVariantAllele().equals(snp.getReferenceAllele())) {
                    snp.setVariantAllele(snp.getReferenceAllele());
                }

            }
        }
        return lowProbSnps;
    }

    private void determineNaHaplotypes(Map<Set<Snp>, List<String>> nanUniqueCombinations, Set<Snp> snps,
                                       PgxSample sample, int haplotype) {
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

    private void logNaHaplotype(BufferedWriter logWriter, Map<Set<Snp>, List<String>> nanUniqueCombinations)
            throws IOException {
        // Sort the output on occurrence
        TreeMap<Set<Snp>, List<String>> sorter = new TreeMap<>(new ListSizeComparator(nanUniqueCombinations));
        sorter.putAll(nanUniqueCombinations);

        for (Set<Snp> snps : sorter.keySet()) {
            logWriter.write("\nNo translation for SNPs:\n");
            for (Snp snp : snps) {
                logWriter.write(snp.getId() + ", reference: " + snp.getReferenceAllele() + ", variant: " +
                        snp.getVariantAllele() + "\n");
            }
            List<String> sampleHaplotypes = nanUniqueCombinations.get(snps);
            logWriter.write(sampleHaplotypes.size() + " affected haplotype(s):\n");
            for (String sampleHaplotype : sampleHaplotypes) {
                logWriter.write(sampleHaplotype + ", ");
            }
            logWriter.write("\n");
        }
    }

    /**
     * Determines the variant snps that were reliably detected (probability at or above the cut-off) on a haplotype
     * that matched a known star allele. These are the actual mutations underlying the star-allele call, as opposed to
     * the low probability snps that get reverted to reference. Must be called after {@link #checkProbabilities}, so
     * that low probability snps have already been reverted to reference and are excluded here.
     */
    private Set<Snp> getDetectedMutations(Set<Snp> snpsOnHaplotype, PgxHaplotype pgxHaplotype) {
        Set<Snp> detectedMutations = new HashSet<>();
        for (Snp snp : snpsOnHaplotype) {
            if (pgxHaplotype.getSnps().containsKey(snp.getId()) &
                    snp.getProbability() >= PROBABILITY_CUT_OFF &&
                    !snp.getVariantAllele().equals(snp.getReferenceAllele())) {
                detectedMutations.add(snp);
            }
        }
        return detectedMutations;
    }

    private void addSnpCombination(Set<Snp> snps, Map<Set<Snp>, List<SnpLog>> snpCombinations,
                                   PgxSample sample, int haplotype) {
        SnpLog snpLog = new SnpLog(snps, sample, haplotype);
        if (!snps.isEmpty()) {
            if (snpCombinations.containsKey(snps)) {
                List<SnpLog> snpLogs = snpCombinations.get(snps);
                snpLogs.add(snpLog);
            } else {
                List<SnpLog> snpLogs = new ArrayList<>();
                snpLogs.add(snpLog);
                snpCombinations.put(snps, snpLogs);
            }
        }
    }

    private void logLowProbSnps(BufferedWriter logWriter, Map<Set<Snp>, List<SnpLog>> lowProbSnpsCombinations)
            throws IOException {
        TreeMap<Set<Snp>, List<SnpLog>> sorter = new TreeMap<>(new LogSizeComparator(lowProbSnpsCombinations));
        sorter.putAll(lowProbSnpsCombinations);

        for (Set<Snp> snps : sorter.keySet()) {
            logWriter.write("\nSamples with low probability for SNPs:\n");
            for (Snp snp : snps) {
                logWriter.write(snp.getId() + ", reference: " + snp.getReferenceAllele() + ", R2: " + snp.getrSquared() +
                        " " + snp.getType() + "\n");
            }
            List<SnpLog> snpLogs = lowProbSnpsCombinations.get(snps);
            Collections.sort(snpLogs);
            logWriter.write(snpLogs.size() + " affected haplotype(s):\n");
            for (SnpLog snpLog : snpLogs) {
                logWriter.write(snpLog.getSample().getId() + "_hap_" + snpLog.getHaplotype());
                for (Snp snp : snpLog.getSnps()) {
                    logWriter.write(";" + snp.getId() + "_" + snp.getOriginalCall() + "_" + snp.getProbability());
                }
                logWriter.write(",");
            }
            logWriter.write("\n");
        }

    }

    private void logDetectedMutations(BufferedWriter logWriter, Map<Set<Snp>, List<SnpLog>> detectedMutationsCombinations)
            throws IOException {
        TreeMap<Set<Snp>, List<SnpLog>> sorter = new TreeMap<>(new LogSizeComparator(detectedMutationsCombinations));
        sorter.putAll(detectedMutationsCombinations);

        for (Set<Snp> snps : sorter.keySet()) {
            logWriter.write("\nSamples with detected mutation(s):\n");
            for (Snp snp : snps) {
                logWriter.write(snp.getId() + ", reference: " + snp.getReferenceAllele() + ", variant: " +
                        snp.getVariantAllele() + ", R2: " + snp.getrSquared() + " " + snp.getType() + "\n");
            }
            List<SnpLog> snpLogs = detectedMutationsCombinations.get(snps);
            Collections.sort(snpLogs);
            logWriter.write(snpLogs.size() + " affected haplotype(s):\n");
            for (SnpLog snpLog : snpLogs) {
                logWriter.write(snpLog.getSample().getId() + "_hap_" + snpLog.getHaplotype());
                for (Snp snp : snpLog.getSnps()) {
                    logWriter.write(";" + snp.getId() + "_" + snp.getVariantAllele() + "_" + snp.getProbability());
                }
                logWriter.write(",");
            }
            logWriter.write("\n");
        }
    }

    /**
     * Loads the configuration
     */
    private void loadConfig() {
        ConfigProvider configProvider = ConfigProvider.getInstance();
        STAR_ALLELE_OUTPUT_DIR = configProvider.getConfigParam(ConfigConstants.STAR_ALLELE_OUTPUT_DIR);
        SNP_HAPLO_TABLE_DIR = configProvider.getConfigParam(ConfigConstants.SNP_HAPLO_TABLE_DIR);
        HAPLOTYPE_READER = new VcfHaplotypeReader(configProvider.getConfigParam(ConfigConstants.HAPLOTYPE_DIR));
        configProvider.containsConfigParam("");
    }

}
