package org.molgenis.asterix.io;

import org.apache.commons.io.FileUtils;
import org.molgenis.asterix.config.ConfigConstants;
import org.molgenis.asterix.config.ConfigProvider;
import org.molgenis.asterix.model.PgxGene;
import org.molgenis.asterix.model.PgxHaplotype;
import org.molgenis.asterix.model.PgxSample;
import org.molgenis.asterix.model.Snp;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.vcf.VcfGenotypeData;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * Implementation of HaplotypeReader for phased VCF files.
 *
 * @author Harm Brugge
 */
public class VcfHaplotypeReader implements HaplotypeReader {

    private Map<String, PgxSample> samples;
    private List<String> sampleNames;
    private String haplotypeFolder;

    static final double R_SQUARED_CUT_OFF = 0.97;

    private static String STAR_ALLELE_OUTPUT_DIR;

    static final int REFERENCE_ALLELE_INDEX = 0;
    static final int ALTERNATIVE_ALLELE_INDEX = 1;
    static final int HAPLOTYPE_0_INDEX = 0;
    static final int HAPLOTYPE_1_INDEX = 1;

    public VcfHaplotypeReader(String haplotypeFolder) {
        this.haplotypeFolder = haplotypeFolder;
        samples = new HashMap<>();
        sampleNames = new ArrayList<>();
        loadConfig();
    }

    @Override
    public void readHaplotypes(Collection<PgxGene> genes) throws IOException {
        File[] phasedVcfFiles = FileUtils.listFiles(new File(haplotypeFolder), new String[]{"vcf.gz"}, true).toArray(new File[0]);
        for (File phasedVcfFile : phasedVcfFiles) {
            System.out.println("Reading genotype file " + phasedVcfFile.getName());
            VcfGenotypeData genotypeData = new VcfGenotypeData(phasedVcfFile, 100, 0.0);
            processSampleIds(genotypeData);

            for (PgxGene pgxGene : genes) {
                if (hasVariantsInRange(pgxGene, genotypeData)) {
                    List<Snp> pgxSnpsInGene = new ArrayList<>(pgxGene.getWildType().getSnps().values());
                    int presentCount = 0;
                    int i = 0;
                    for (Snp pgxSnp : pgxSnpsInGene) {
                        i++;
                        System.out.print("Processing PGx gene " + pgxGene.getName() + " SNP " + i + "/" + pgxSnpsInGene.size() + "\r");
                        GeneticVariant snp = genotypeData.getSnpVariantByPos(Integer.toString(pgxSnp.getChr()), pgxSnp.getPos());
                        if (snp != null) {
                            presentCount++;
                            parseSnp(snp, pgxSnp, pgxGene);
                        } else {
                            pgxSnp.setAvailable(false);
//                            removePgxSnpFromGene(pgxGene, pgxSnp);
                        }
                    }
                    System.out.println(presentCount + "/" + pgxGene.getWildType().getSnps().values().size() + " SNPs present");
                }
                writeOutputFile(pgxGene);
            }

        }
    }

    @Override
    public void readHaplotypes() throws IOException {

    }

    @Override
    public Map<String, PgxSample> getSamples() {
        return samples;
    }

    private void processSampleIds(VcfGenotypeData vcfGenotypeData) {
        if (samples.isEmpty()) {
            for (Sample sample : vcfGenotypeData.getSamples()) {
                if (!samples.containsKey(sample.getId())) {
                    samples.put(sample.getId(), new PgxSample(sample.getId()));
                    sampleNames.add(sample.getId());
                } else {
                    throw new RuntimeException("Duplicate sample in genotype file");
                }
            }
        } else if (!isEqualSamples(vcfGenotypeData.getSamples())) {
            throw new RuntimeException("Samples are not equal across genotype files");
        }
    }

    private void parseSnp(GeneticVariant variant, Snp pgxSnp, PgxGene pgxGene) {

        pgxSnp.setrSquared(Double.parseDouble(variant.getAnnotationValues().get("R2").toString()));
        pgxSnp.setMaf(Double.parseDouble(variant.getAnnotationValues().get("MAF").toString()));

        pgxSnp.setAvailable(pgxSnp.getrSquared() > R_SQUARED_CUT_OFF);
        pgxGene.updateSnpInfoOnHaplotypes(pgxSnp);
        //            removePgxSnpFromGene(pgxGene, pgxSnp);

        Snp referenceSnp = pgxSnp.copySnp(pgxSnp);
        referenceSnp.setReferenceAllele(variant.getRefAllele().getAlleleAsString());
        referenceSnp.setMinorAllele(variant.getMinorAllele().getAlleleAsString());
        referenceSnp.setVariantAllele(variant.getRefAllele().getAlleleAsString());

        Snp minorSnp = referenceSnp.copySnp(referenceSnp);
        minorSnp.setVariantAllele(variant.getMinorAllele().getAlleleAsString());

        double[][][] sampleGenotypeProbabilitiesPhased = variant.getSampleGenotypeProbabilitiesPhased();
        for (int i = 0; i < sampleGenotypeProbabilitiesPhased.length; i++) {
            String currentSampleId = sampleNames.get(i);
            PgxSample currentSample = samples.get(currentSampleId);

            double[][] sample = sampleGenotypeProbabilitiesPhased[i];

            // TODO: Should implement min max for probability?
            double probabilityReferenceAlleleHaplotype0 = sample[HAPLOTYPE_0_INDEX][REFERENCE_ALLELE_INDEX];
            if (probabilityReferenceAlleleHaplotype0 > 0.5) {
                referenceSnp.setProbability(probabilityReferenceAlleleHaplotype0);
                currentSample.getHaplotype0().put(referenceSnp.getId(), referenceSnp);
            }
            else {
                minorSnp.setProbability(sample[HAPLOTYPE_0_INDEX][ALTERNATIVE_ALLELE_INDEX]);
                currentSample.getHaplotype0().put(minorSnp.getId(), minorSnp);
            }

            double probabilityReferenceAlleleHaplotype1 = sample[HAPLOTYPE_1_INDEX][REFERENCE_ALLELE_INDEX];
            if (probabilityReferenceAlleleHaplotype1 > 0.5) {
                Snp snp = referenceSnp.copySnp(referenceSnp);
                snp.setProbability(probabilityReferenceAlleleHaplotype1);
                currentSample.getHaplotype1().put(snp.getId(), snp);
            }
            else {
                Snp snp = minorSnp.copySnp(minorSnp);
                snp.setProbability(sample[HAPLOTYPE_1_INDEX][ALTERNATIVE_ALLELE_INDEX]);
                currentSample.getHaplotype1().put(minorSnp.getId(), minorSnp);
            }
        }
    }

    private boolean hasVariantsInRange(PgxGene gene, VcfGenotypeData genotypeData) {
        String chr = Integer.toString(gene.getChr());
        int start = gene.getStartPos();
        int end = gene.getEndPos();
        return genotypeData.getVariantsByRange(chr, start, end).iterator().hasNext();
    }

    private void removePgxSnpFromGene(PgxGene pgxGene, Snp pgxSnp) {
        pgxGene.getWildType().removeSnp(pgxSnp);

        // Remove SNP on all haplotypes and remove haplotype when empty
        for (Iterator<Map.Entry<String, PgxHaplotype>> it = pgxGene.getPgxHaplotypes().entrySet().iterator(); it.hasNext(); ) {
            PgxHaplotype haplotype = it.next().getValue();
            haplotype.removeSnp(pgxSnp);
            if (haplotype.getSnps().size() == 0) it.remove();
        }
    }

    private boolean isEqualSamples(List<Sample> samples) {
        return true;
    }

    private void writeOutputFile(PgxGene pgxGene) throws IOException {
        // Write output file
        File pgxGeneInfoFile = new File(STAR_ALLELE_OUTPUT_DIR + pgxGene.getName() + "_info.txt");
        if (!pgxGeneInfoFile.exists()) pgxGeneInfoFile.createNewFile();

        FileWriter fileWriter = new FileWriter(pgxGeneInfoFile.getAbsoluteFile());
        BufferedWriter bw = new BufferedWriter(fileWriter);

        bw.write("Haplotype\tSNP\tPosition\tReference_allele\tVariant_allele\tStatus\tR_squared\tMAF\n");

        for (PgxHaplotype pgxHaplotype : pgxGene.getPgxHaplotypes().values()) {
            for (Snp haplotypeSnp : pgxHaplotype.getSnps().values()) {
                String isAvailable = "NOT AVAILABLE";
                if (haplotypeSnp.isAvailable()) isAvailable = "AVAILABLE";
                else if (haplotypeSnp.getrSquared() != null &&
                        haplotypeSnp.getrSquared() < R_SQUARED_CUT_OFF)
                    isAvailable = "LOW_IMPUTATION_PROBABILITY";
                bw.write(pgxHaplotype.getName() + "\t" + haplotypeSnp.getId() + "\t" + haplotypeSnp.getPos() + "\t");
                bw.write(haplotypeSnp.getReferenceAllele() + "\t" + haplotypeSnp.getVariantAllele() + "\t");
                bw.write(isAvailable + "\t" + haplotypeSnp.getrSquared() + "\t" + haplotypeSnp.getMaf() + "\n");
            }
        }

        bw.close();
    }

    private void loadConfig() {
        //set config provider
        ConfigProvider configProvider = ConfigProvider.getInstance();
        //load dirs
        STAR_ALLELE_OUTPUT_DIR = configProvider.getConfigParam(ConfigConstants.STAR_ALLELE_OUTPUT_DIR);
        STAR_ALLELE_OUTPUT_DIR = "/groups/umcg-gdio/tmp01/projects/2021001/pgx-pipeline/analyses/impute_pgx_genes/out/20221103/postimpute/";
    }

}
