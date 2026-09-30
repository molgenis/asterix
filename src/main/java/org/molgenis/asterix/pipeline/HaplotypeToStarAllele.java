package org.molgenis.asterix.pipeline;

import org.molgenis.asterix.config.ConfigConstants;
import org.molgenis.asterix.config.ConfigProvider;
import org.molgenis.asterix.io.ApprovalReader;
import org.molgenis.asterix.io.ApprovalWriter;
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
    private static String APPROVAL_INPUT_FILE;
    private static String APPROVAL_OUTPUT_FILE;
    private static HaplotypeReader HAPLOTYPE_READER;

    // Gene identifiers as keys
    private SortedMap<String, PgxGene> genes;
    // Sample identifiers as keys
    private Map<String, PgxSample> samples;
    // The expert decisions of previous review rounds, extended with the cases found in this run
    private ApprovalRegistry approvalRegistry;


    public HaplotypeToStarAllele() throws IOException {
        this.loadConfig();
        this.loadApprovalRegistry();

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

            // Haplotypes Asterix is not certain about, collected here and resolved against the
            // approval registry once the CYP2D6 CNV calls are in
            List<ReviewCase> reviewCases = new ArrayList<>();

            PgxGene pgxGene = genes.get(key);

            int duplicateCount = 0;
            int naCount = 0;

            for (PgxSample sample : samples.values()) {

                sample.getGenes().put(pgxGene.getName(), new PgxGene(pgxGene));

                String starAllele0 = null;
                String starAllele1 = null;

                String disclaimer0 = null;
                String disclaimer1 = null;

                boolean hasDuplicate = false;

                Set<Snp> snpsOnHaplotype0 = sample.getHaplotype0(pgxGene.getVariants());
                Set<Snp> snpsOnHaplotype1 = sample.getHaplotype1(pgxGene.getVariants());

                Set<Snp> lowProbSnpsOn0 = new HashSet<>();
                Set<Snp> lowProbSnpsOn1 = new HashSet<>();

                for (PgxHaplotype pgxHaplotype : pgxGene.getPgxHaplotypes().values()) {
                    if (pgxHaplotype.getSnps().isEmpty()) continue;
                    Set<Snp> snpsOnPgxHaploType = new HashSet<>(pgxHaplotype.getSnpsWithWildType().values());

                    if (snpsOnHaplotype0.equals(snpsOnPgxHaploType)) {
                        lowProbSnpsOn0.addAll(checkProbabilities(snpsOnHaplotype0, pgxHaplotype));
                        if (starAllele0 != null) {
                            logWriter.write(sample.getId() + ": duplicate haplo 0 " + starAllele0 + " " +
                                    pgxHaplotype.getName() + "\n");
                            hasDuplicate = true;
                        }
                        starAllele0 = pgxHaplotype.getCorrected();
                        disclaimer0 = pgxHaplotype.getDisclaimer();
                        sample.getGenes().get(pgxGene.getName()).setDisclaimerSubject(pgxHaplotype.getDisclaimer());
                        sample.getGenes().get(pgxGene.getName()).setAllele0(starAllele0);
                    }

                    if (snpsOnHaplotype1.equals(snpsOnPgxHaploType)) {
                        lowProbSnpsOn1.addAll(checkProbabilities(snpsOnHaplotype1, pgxHaplotype));
                        if (starAllele1 != null) {
                            logWriter.write(sample.getId() + ": duplicate haplo 1 " + starAllele1 + " " +
                                    pgxHaplotype.getName() + "\n");
                            hasDuplicate = true;
                        }
                        starAllele1 = pgxHaplotype.getCorrected();
                        disclaimer1 = pgxHaplotype.getDisclaimer();
                        sample.getGenes().get(pgxGene.getName()).setDisclaimerSubject(pgxHaplotype.getDisclaimer());
                        sample.getGenes().get(pgxGene.getName()).setAllele1(starAllele1);
                    }

                    if (starAllele0 != null & starAllele1 != null) {
                        //break;
                    }
                }

                addLowProbSnps(lowProbSnpsOn0, lowProbSnpsCombinations, sample, 0);
                addLowProbSnps(lowProbSnpsOn1, lowProbSnpsCombinations, sample, 1);

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

                addReviewCase(reviewCases, sample, 0, snpsOnHaplotype0, lowProbSnpsOn0, starAllele0, disclaimer0);
                addReviewCase(reviewCases, sample, 1, snpsOnHaplotype1, lowProbSnpsOn1, starAllele1, disclaimer1);
            }
            System.out.println(pgxGene.getName() + ": " + duplicateCount + " duplicates");
            System.out.println(pgxGene.getName() + ": " + naCount + " nans");

            logNaHaplotype(logWriter, nanUniqueCombinations);
            logLowProbSnps(logWriter, lowProbSnpsCombinations);

            if (pgxGene.getName().equals("CYP2D6")) {
                new Cyp2d6Caller(samples, pgxGene);
            }

            // After the CNV caller, so a CYP2D6 case is reviewed together with its CNV status and
            // an expert decision has the final say over the CNV call
            applyApprovalDecisions(pgxGene, reviewCases, logWriter);

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

    private void addLowProbSnps(Set<Snp> lowProbSnps, Map<Set<Snp>, List<SnpLog>> lowProbSnpsCombinations,
                                PgxSample sample, int haplotype) {
        SnpLog snpLog = new SnpLog(lowProbSnps, sample, haplotype);
        if (!lowProbSnps.isEmpty()) {
            if (lowProbSnpsCombinations.containsKey(lowProbSnps)) {
                List<SnpLog> snpLogs = lowProbSnpsCombinations.get(lowProbSnps);
                snpLogs.add(snpLog);
            } else {
                List<SnpLog> snpLogs = new ArrayList<>();
                snpLogs.add(snpLog);
                lowProbSnpsCombinations.put(lowProbSnps, snpLogs);
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

    /**
     * Record a haplotype for expert review when Asterix is not certain about it: either no
     * translation was found for its SNPs, or the call rests on SNPs that were imputed with a
     * probability below the cut-off. Haplotypes that matched a translation cleanly need no review.
     *
     * @param reviewCases     the cases collected for this gene so far
     * @param sample          the sample the haplotype belongs to
     * @param haplotype       0 or 1
     * @param snpsOnHaplotype the SNPs on the haplotype, with failed calls already reverted to reference
     * @param lowProbSnps     the SNPs that were called with too low a probability to rely on
     * @param starAllele      the star allele that was matched, or null when there was no match
     * @param disclaimer      the disclaimer of the matched star allele, may be null
     */
    private void addReviewCase(List<ReviewCase> reviewCases, PgxSample sample, int haplotype,
                               Set<Snp> snpsOnHaplotype, Set<Snp> lowProbSnps, String starAllele,
                               String disclaimer) {
        if (starAllele != null && lowProbSnps.isEmpty()) return;

        Set<String> failedIds = new HashSet<>();
        for (Snp snp : lowProbSnps) {
            failedIds.add(snp.getId());
        }

        Set<String> variants = new TreeSet<>();
        for (Snp snp : snpsOnHaplotype) {
            if (failedIds.contains(snp.getId())) continue;
            if (!snp.getVariantAllele().equals(snp.getReferenceAllele())) {
                variants.add(HaplotypeCombination.describePresentVariant(snp.getId(), snp.getVariantAllele()));
            }
        }
        for (String failedId : failedIds) {
            variants.add(HaplotypeCombination.describeFailedCall(failedId));
        }

        if (variants.isEmpty()) return;

        reviewCases.add(new ReviewCase(sample, haplotype, variants, disclaimer));
    }

    /**
     * Resolve every uncertain haplotype of a gene against the approval registry.
     *
     * Combinations the registry does not know yet are added to it as unreviewed rows, which is
     * what puts them in the workbook for the next review round; the call Asterix made stands until
     * a reviewer decides otherwise. Combinations that were reviewed get the approved call and
     * disclaimer, overruling what Asterix determined, including the CYP2D6 CNV call.
     *
     * @param pgxGene     the gene being processed
     * @param reviewCases the uncertain haplotypes collected for this gene
     * @param logWriter   the per-gene log, used to record which sample got which decision
     * @throws IOException when the log cannot be written
     */
    private void applyApprovalDecisions(PgxGene pgxGene, List<ReviewCase> reviewCases, BufferedWriter logWriter)
            throws IOException {
        int appliedCount = 0;
        int unreviewedCount = 0;
        StringBuilder appliedLog = new StringBuilder();

        for (ReviewCase reviewCase : reviewCases) {
            PgxGene sampleGene = reviewCase.sample.getGenes().get(pgxGene.getName());

            HaplotypeCombination combination = new HaplotypeCombination(pgxGene.getName(),
                    reviewCase.variants, sampleGene.getCnvStatus());

            String currentCall = reviewCase.haplotype == 0 ? sampleGene.getAllele0() : sampleGene.getAllele1();
            ApprovalDecision decision = approvalRegistry.observe(combination,
                    ApprovalDecision.toBareStarAllele(pgxGene.getName(), currentCall), reviewCase.disclaimer);

            String decidedCall = decision.getDecidedHaplotype();
            if (decidedCall == null) {
                unreviewedCount += 1;
                continue;
            }

            if (reviewCase.haplotype == 0) sampleGene.setAllele0(decidedCall);
            else sampleGene.setAllele1(decidedCall);
            sampleGene.setDisclaimerSubject(decision.getEffectiveDisclaimer());
            appliedCount += 1;

            appliedLog.append(reviewCase.sample.getId()).append("_hap_").append(reviewCase.haplotype)
                    .append(": ").append(currentCall).append(" -> ").append(decidedCall)
                    .append(" (").append(decision.isModified() ? "modified" : "approved")
                    .append(" by ").append(decision.getLoggedBy())
                    .append(" on ").append(decision.getDate()).append(")\n");
        }

        if (appliedCount > 0) {
            logWriter.write("\nApproval decisions applied to " + appliedCount + " haplotype(s):\n");
            logWriter.write(appliedLog.toString());
        }
        if (unreviewedCount > 0) {
            logWriter.write("\n" + unreviewedCount + " haplotype(s) await expert review, see " +
                    APPROVAL_OUTPUT_FILE + "\n");
        }

        System.out.println(pgxGene.getName() + ": " + appliedCount + " approval decision(s) applied, " +
                unreviewedCount + " haplotype(s) awaiting review");
    }

    /**
     * Write the approval (vrijgifte) workbook: the decisions that were reused plus the
     * combinations found in this run that still need one.
     *
     * @throws IOException when the workbook cannot be written
     */
    public void writeApprovalDecisions() throws IOException {
        new ApprovalWriter().write(new File(APPROVAL_OUTPUT_FILE), approvalRegistry, genes);
    }

    public ApprovalRegistry getApprovalRegistry() {
        return approvalRegistry;
    }

    /**
     * Load the decisions of previous review rounds, when an approval file was configured.
     *
     * @throws IOException when the configured file cannot be read
     */
    private void loadApprovalRegistry() throws IOException {
        approvalRegistry = new ApprovalRegistry();

        if (APPROVAL_INPUT_FILE == null || APPROVAL_INPUT_FILE.trim().isEmpty()) {
            System.out.println("No approval decision file configured (use -approval_in), " +
                    "uncertain combinations will only be reported for review");
            return;
        }

        new ApprovalReader().read(new File(APPROVAL_INPUT_FILE), approvalRegistry);
    }

    /**
     * Loads the configuration
     */
    private void loadConfig() {
        ConfigProvider configProvider = ConfigProvider.getInstance();
        STAR_ALLELE_OUTPUT_DIR = configProvider.getConfigParam(ConfigConstants.STAR_ALLELE_OUTPUT_DIR);
        // The haplotype reader already writes its per-gene info files here, so the directory has
        // to exist before it runs, not only by the time writeStarAlleles() is reached
        new File(STAR_ALLELE_OUTPUT_DIR).mkdirs();
        SNP_HAPLO_TABLE_DIR = configProvider.getConfigParam(ConfigConstants.SNP_HAPLO_TABLE_DIR);
        APPROVAL_INPUT_FILE = configProvider.getConfigParam(ConfigConstants.APPROVAL_INPUT_FILE);
        APPROVAL_OUTPUT_FILE = configProvider.getConfigParam(ConfigConstants.APPROVAL_OUTPUT_FILE);
        HAPLOTYPE_READER = new VcfHaplotypeReader(configProvider.getConfigParam(ConfigConstants.HAPLOTYPE_DIR));
        configProvider.containsConfigParam("");
    }

    /**
     * One haplotype of one sample that Asterix is not certain about. Held until the CYP2D6 CNV
     * calls are in, because the CNV status is part of the identity of a case.
     */
    private static class ReviewCase {

        private final PgxSample sample;
        private final int haplotype;
        private final Set<String> variants;
        private final String disclaimer;

        private ReviewCase(PgxSample sample, int haplotype, Set<String> variants, String disclaimer) {
            this.sample = sample;
            this.haplotype = haplotype;
            this.variants = variants;
            this.disclaimer = disclaimer;
        }
    }

}
