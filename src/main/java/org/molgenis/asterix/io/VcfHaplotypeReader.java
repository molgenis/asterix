package org.molgenis.asterix.io;

import org.apache.commons.io.FileUtils;
import org.molgenis.asterix.config.ConfigConstants;
import org.molgenis.asterix.config.ConfigProvider;
import org.molgenis.asterix.model.PgxGene;
import org.molgenis.asterix.model.PgxHaplotype;
import org.molgenis.asterix.model.PgxSample;
import org.molgenis.asterix.model.Snp;
import org.molgenis.genotype.Alleles;
import org.molgenis.genotype.Sample;
import org.molgenis.genotype.annotation.Annotation;
import org.molgenis.genotype.variant.GeneticVariant;
import org.molgenis.genotype.vcf.VcfGenotypeData;
import org.molgenis.genotype.vcf.VcfGenotypeField.VcfGenotypeFormat;
import sun.reflect.generics.reflectiveObjects.NotImplementedException;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.Stream;

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
            genotypeData.setPreferredGenotypeFormat(VcfGenotypeFormat.ADS, "HDS");

            Map<String, Annotation> variantAnnotationsMap = genotypeData.getVariantAnnotationsMap();
            processSampleIds(genotypeData);

            for (PgxGene pgxGene : genes) {
                if (hasVariantsInRange(pgxGene, genotypeData)) {
                    List<Snp> pgxSnpsInGene = new ArrayList<>(pgxGene.getVariants().values());
                    int presentCount = 0;
                    int i = 0;
                    for (Snp pgxSnp : pgxSnpsInGene) {
                        i++;
                        System.out.print("Processing PGx gene " + pgxGene.getName() + " SNP " + i + "/" + pgxSnpsInGene.size() + "\r");
                        GeneticVariant genotypesVariant = matchGeneticVariant(genotypeData, pgxSnp);
                        if (genotypesVariant == null) {
                            System.out.println("SNP " + pgxGene.getName() + " " + pgxSnp.getId() + " not available");
                            pgxSnp.setAvailable(false);
                        } else {
                            presentCount++;
                            parseSnp(genotypesVariant, variantAnnotationsMap, pgxSnp, pgxGene);
                        }
                    }
                    System.out.println(presentCount + "/" + pgxGene.getVariants().values().size() + " SNPs present");
                }
                writeOutputFile(pgxGene);
            }

        }
    }

    private GeneticVariant matchGeneticVariant(VcfGenotypeData genotypeData, Snp pgxSnp) {
        GeneticVariant matchedVariant = null;
        Alleles pgxSnpAlleles = Alleles.createBasedOnString(pgxSnp.getReferenceAllele(), pgxSnp.getVariantAllele());
        Iterable<GeneticVariant> positionalMatches = genotypeData.getVariantsByPos(Integer.toString(pgxSnp.getChr()), pgxSnp.getPos());
        for (GeneticVariant candidateVariant : positionalMatches) {
            if (candidateVariant.getVariantAlleles().sameAlleles(pgxSnpAlleles)) {
                if (matchedVariant != null) {
                    System.out.println("pgxSnp = " + pgxSnp);
                    throw new RuntimeException("Multiple SNPs match PGx SNP");
                }
                matchedVariant = candidateVariant;
            }
        }
        return matchedVariant;
    }

    @Override
    public void readHaplotypes() throws IOException {
        throw new NotImplementedException();
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

    private void parseSnp(GeneticVariant variant, Map<String, Annotation> variantAnnotationsMap, Snp pgxSnp, PgxGene pgxGene) {
        Map<String, String> mappedAnnotationValues = new HashMap<>();
        Map<String, ?> annotationValues = variant.getAnnotationValues();
        pgxSnp.setrSquared(((Number) annotationValues.get("R2")).doubleValue());

        for (String annotationId : annotationValues.keySet()) {
            if (variantAnnotationsMap.get(annotationId) != null) {
                if (variantAnnotationsMap.get(annotationId).isList()) {
                    List<?> o = (List<?>) annotationValues.get(annotationId);
                    Object o1 = o.get(0);
                    mappedAnnotationValues.put(annotationId, String.valueOf(o1));
                } else {
                    mappedAnnotationValues.put(annotationId, String.valueOf(annotationValues.get(annotationId)));
                }
            }
        }

        pgxSnp.setAnnotations(mappedAnnotationValues);
        pgxSnp.setHwe(variant.getHwePvalue());
        pgxSnp.setMaf(variant.getMinorAlleleFrequency());

        //pgxSnp.setAvailable(pgxSnp.getrSquared() > R_SQUARED_CUT_OFF);
        pgxSnp.setAvailable(true);
        pgxGene.updateSnpInfo(pgxSnp);
        //removePgxSnpFromGene(pgxGene, pgxSnp);

        Snp referenceSnp = pgxSnp.copySnp(pgxSnp);

        referenceSnp.setReferenceAllele(variant.getRefAllele().getAlleleAsString());
        referenceSnp.setMinorAllele(variant.getMinorAllele().getAlleleAsString());
        Snp alternativeSnp = referenceSnp.copySnp(referenceSnp);

        referenceSnp.setVariantAllele(variant.getRefAllele().getAlleleAsString());
        alternativeSnp.setVariantAllele(variant.getVariantAlleles().get(1).getAlleleAsString());

        double[][][] sampleGenotypeProbabilitiesPhased = variant.getSampleGenotypeProbabilitiesPhased();
        for (int i = 0; i < sampleGenotypeProbabilitiesPhased.length; i++) {
            String currentSampleId = sampleNames.get(i);
            PgxSample currentSample = samples.get(currentSampleId);

            double[][] sample = sampleGenotypeProbabilitiesPhased[i];

            double probabilityReferenceAlleleHaplotype0 = sample[HAPLOTYPE_0_INDEX][REFERENCE_ALLELE_INDEX];
            if (probabilityReferenceAlleleHaplotype0 > 0.5) {
                Snp snp = referenceSnp.copySnp(referenceSnp);
                snp.setProbability(probabilityReferenceAlleleHaplotype0);
                currentSample.getHaplotype0().put(snp.getId(), snp);
            } else {
                Snp snp = alternativeSnp.copySnp(alternativeSnp);
                snp.setProbability(sample[HAPLOTYPE_0_INDEX][ALTERNATIVE_ALLELE_INDEX]);
                currentSample.getHaplotype0().put(snp.getId(), snp);
            }

            double probabilityReferenceAlleleHaplotype1 = sample[HAPLOTYPE_1_INDEX][REFERENCE_ALLELE_INDEX];
            if (probabilityReferenceAlleleHaplotype1 > 0.5) {
                Snp snp = referenceSnp.copySnp(referenceSnp);
                snp.setProbability(probabilityReferenceAlleleHaplotype1);
                currentSample.getHaplotype1().put(snp.getId(), snp);
            } else {
                Snp snp = alternativeSnp.copySnp(alternativeSnp);
                snp.setProbability(sample[HAPLOTYPE_1_INDEX][ALTERNATIVE_ALLELE_INDEX]);
                currentSample.getHaplotype1().put(alternativeSnp.getId(), snp);
            }
        }
    }

    private boolean hasVariantsInRange(PgxGene gene, VcfGenotypeData genotypeData) {
        String chr = Integer.toString(gene.getChr());
        int start = gene.getStartPos() - 1;
        int end = gene.getEndPos() + 1;
        return genotypeData.getVariantsByRange(chr, start, end).iterator().hasNext();
    }

    private void removePgxSnpFromGene(PgxGene pgxGene, Snp pgxSnp) {
        pgxGene.removeVariant(pgxSnp);
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
        //"Haplotype", "SNP", "Position", "Reference_allele", "Variant_allele", "Status", "R_squared", "MAF"
        String[] fixedColumns = new String[]{
                "Haplotype", "SNP", "Position", "Reference_allele", "Variant_allele", "Status", "Calculated_Hwe", "Calculated_Maf"};
        List<String> infoColumns = pgxGene.getAnnotationFields();

        String header = Stream.concat(Stream.of(fixedColumns), infoColumns.stream().map(e -> "Info_" + e))
                .collect(Collectors.joining("\t")).concat("\n");
        bw.write(header);

        for (PgxHaplotype pgxHaplotype : pgxGene.getPgxHaplotypes().values()) {
            Map<String, String> variantAlleles = pgxHaplotype.getVariantAlleles();
            for (Snp haplotypeSnp : pgxHaplotype.getSnps().values()) {
                String variantAllele = variantAlleles.get(haplotypeSnp.getId());
                String isAvailable = "NOT AVAILABLE";
                if (haplotypeSnp.isAvailable()) isAvailable = "AVAILABLE";
                else if (haplotypeSnp.getrSquared() != null &&
                        haplotypeSnp.getrSquared() < R_SQUARED_CUT_OFF)
                    isAvailable = "LOW_IMPUTATION_PROBABILITY";

                String mafFormatted = "na";
                Double maf = haplotypeSnp.getMaf();
                if (maf != null) {
                    mafFormatted = Double.toString(maf);
                }

                List<String> outputValues = new ArrayList<>(Arrays.asList(
                        pgxHaplotype.getName(),
                        haplotypeSnp.getId(),
                        Integer.toString(haplotypeSnp.getPos()),
                        haplotypeSnp.getReferenceAllele(),
                        variantAllele,
                        isAvailable,
                        Double.toString(haplotypeSnp.getHwe()),
                        mafFormatted));

                // Write:
                // - availability
                // - reported MAF
                // - calculated MAF
                // - MAF from reference
                // - calculated HWE
                // - R2
                Map<String, ?> annotations = haplotypeSnp.getAnnotations();
                for (String column : infoColumns) {
                    if (annotations.containsKey(column)) {
                        outputValues.add(String.valueOf(annotations.get(column)));
                    } else {
                        outputValues.add("na");
                    }
                }

                String line = String.join("\t", outputValues).concat("\n");
                bw.write(line);

            }
        }

        bw.close();
    }

    private void loadConfig() {
        //set config provider
        ConfigProvider configProvider = ConfigProvider.getInstance();
        //load dirs
        STAR_ALLELE_OUTPUT_DIR = configProvider.getConfigParam(ConfigConstants.STAR_ALLELE_OUTPUT_DIR);
    }

}
