package org.molgenis.asterix.pipeline;

import org.molgenis.asterix.config.ConfigConstants;
import org.molgenis.asterix.config.ConfigProvider;
import org.molgenis.asterix.model.PGxDiplotype;
import org.molgenis.asterix.model.PgxGene;
import org.molgenis.asterix.model.PgxHaplotype;
import org.molgenis.asterix.model.PgxSample;

import java.io.*;
import java.util.*;

public class StarAlleleToPhenotype {


    //internal value, to switch from having the value names in a vertical column or as header
    private static boolean VALUE_NAME_VERTICAL = false;

    //input dir for converting haplotype to phenotype (star alleles to functions)
    private static String CONVERSION_TABLE_DIR;
    //output dir for converted haplotypes to predicted functions (per gene as file, with rows as persons)
    private static String PREDICTED_PHENOTYPES_OUTPUT_DIR;
    //output file for sample matrix
    private static String SAMPLE_MATRIX_OUT;
    //split sample matrix per sample
    private static boolean SPLIT_SAMPLE_MATRIX_PP;

    //config provider for command line supplied parameters
    private ConfigProvider configProvider = null;

    private SortedMap<String, PgxGene> genes;

    private Map<Set<PgxHaplotype>, PGxDiplotype> diplotypes = new HashMap<>();

    /**
     * @param genes
     */
    public StarAlleleToPhenotype(SortedMap<String, PgxGene> genes) {
        this.genes = genes;
        //set config provider
        this.configProvider = ConfigProvider.getInstance();
        //load dirs
        CONVERSION_TABLE_DIR = this.configProvider.getConfigParam(ConfigConstants.HAPLO_PHENO_TABLE_DIR);
        PREDICTED_PHENOTYPES_OUTPUT_DIR = this.configProvider.getConfigParam(ConfigConstants.PREDICTED_PHENOTYPES_OUTPUT_DIR);
        SAMPLE_MATRIX_OUT = this.configProvider.getConfigParam(ConfigConstants.SAMPLE_MATRIX_OUT);
        SPLIT_SAMPLE_MATRIX_PP = Boolean.parseBoolean(this.configProvider.getConfigParam(ConfigConstants.SPLIT_SAMPLES_PP));
    }

    private void readConversionTables() throws IOException {

        File[] phenotypeTables = new File(CONVERSION_TABLE_DIR).listFiles();

        for (File phenotypeTable : phenotypeTables) {
            try (BufferedReader br = new BufferedReader(new FileReader(phenotypeTable))) {
                br.readLine(); // skip header
                String line = br.readLine();
                readLine:
                while (line != null) {
                    String[] splitLine = line.split("\t");

                    String geneName = splitLine[0];
                    PgxGene pgxGene = genes.get(geneName);

                    String[] diplotypeArray = splitLine[1].split("/");

                    Set<PgxHaplotype> haplotypes = new HashSet<>();
                    for (String haplotypeId : diplotypeArray) {
                        String haplotypeName = pgxGene.getName() + haplotypeId;
                        PgxHaplotype haplotype = pgxGene.getPgxHaplotypes().get(haplotypeName);
                        if (haplotype == null) {
                            System.out.println("Haplotype " + haplotypeName + " in function table, not available in snp table");
                            line = br.readLine();
                            continue readLine;
                        }
                        haplotypes.add(haplotype);
                    }

                    PGxDiplotype diplotype = new PGxDiplotype();
                    diplotype.setPredictedPhenotype(splitLine[3]);
                    diplotype.setContraindication(splitLine[4]);
                    diplotype.setHaplotypes(haplotypes);

                    this.diplotypes.put(haplotypes, diplotype);
                    line = br.readLine();
                }
            }
        }
    }

    /**
     * determine the phenotypes of the samples from their star alleles (values will be set to the
     *
     * @param haplotypeToStarAllele Object to get star alleles from, and set the predicted phenotype for
     * @throws IOException
     */
    public void determinePhenotypes(HaplotypeToStarAllele haplotypeToStarAllele) throws IOException {
        this.readConversionTables();

        Map<String, PgxSample> samples = haplotypeToStarAllele.getSamples();

        for (PgxSample sample : samples.values()) {
            for (PgxGene gene : sample.getGenes().values()) {
                //PgxGene geneWithInfo = this.genes.get(gene.getName());
                if (gene.getAllele0() == null | gene.getAllele1() == null) continue;
                Set<PgxHaplotype> haplotypes = new HashSet<>();
                haplotypes.add(new PgxHaplotype(gene, gene.getAllele0()));
                haplotypes.add(new PgxHaplotype(gene, gene.getAllele1()));

                PGxDiplotype diplotype = diplotypes.get(haplotypes);
                if (diplotype == null) continue;

                gene.setPredictedPhenotype(diplotype.getPredictedPhenotype());
            }
        }


    }

    /**
     * write the phenotypes that were predicted
     *
     * @param haplotypeToStarAllele Object containing the samples and their predicted phenotypes
     * @throws IOException
     */
    public void writePhenotypes(HaplotypeToStarAllele haplotypeToStarAllele) throws IOException {
        new File(PREDICTED_PHENOTYPES_OUTPUT_DIR).mkdirs();

        Map<String, PgxSample> samples = haplotypeToStarAllele.getSamples();

        writeSampleMatrix(haplotypeToStarAllele, samples);

        for (PgxGene pgxGene : haplotypeToStarAllele.getGenes().values()) {

            File predictedPhenotypeOutputFile = new File(PREDICTED_PHENOTYPES_OUTPUT_DIR + pgxGene.getName() + "_predicted_phenotypes.txt");
            if (!predictedPhenotypeOutputFile.exists()) predictedPhenotypeOutputFile.createNewFile();

            FileWriter fileWriter = new FileWriter(predictedPhenotypeOutputFile.getAbsoluteFile());
            BufferedWriter bw = new BufferedWriter(fileWriter);

            for (PgxSample sample : samples.values()) {
                bw.write(sample.getId() + "\t" + sample.getGenes().get(pgxGene.getName()).getPredictedPhenotype() + "\n");
            }

            bw.close();
        }
    }

    /**
     * write the sample matrix
     *
     * @param haplotypeToStarAllele the object containing the phenotype and haplotypes
     * @param samples               the list of samples to write
     * @throws IOException
     */
    private static void writeSampleMatrix(HaplotypeToStarAllele haplotypeToStarAllele, Map<String, PgxSample> samples) throws IOException {
        if (SPLIT_SAMPLE_MATRIX_PP) {
            writeSampleMatrixSplit(haplotypeToStarAllele, samples);
        } else {
            writeSampleMatrixCombined(haplotypeToStarAllele, samples);
        }
    }

    /**
     * write the sample matrix, with a new file for each sample
     *
     * @param haplotypeToStarAllele the object containing the phenotype and haplotypes
     * @param samples               the list of samples to write
     * @throws IOException
     */
    private static void writeSampleMatrixSplit(HaplotypeToStarAllele haplotypeToStarAllele, Map<String, PgxSample> samples) throws IOException {
        //decide on a vertical or horizontal values file
        if (VALUE_NAME_VERTICAL) {
            writeSampleMatrixSplitVertical(haplotypeToStarAllele, samples);
        } else {
            writeSampleMatrixSplitHorizontal(haplotypeToStarAllele, samples);
        }
    }

    /*
    following method \/ outputs in the following format

    ,COMT_hap0,COMT_hap1,CYP1A2_hap0,CYP1A2_hap1
    AA12345,COMTHaplotype low activity,COMTHaplotype normal activity,CYP1A2*1A,CYP1A2*1A

     */

    /**
     * write the sample matrix per sample, with a header and the values, which means a two line file is created
     *
     * @param haplotypeToStarAllele the object containing the phenotype and haplotypes
     * @param samples               the list of samples to write
     * @throws IOException
     */
    private static void writeSampleMatrixSplitHorizontal(HaplotypeToStarAllele haplotypeToStarAllele, Map<String, PgxSample> samples) throws IOException {
        //get the header
        String header = getSampleMatrixHeader(haplotypeToStarAllele);
        //using stringbuilder is more efficient than keep concatenating strings
        StringBuilder stringBuilder = new StringBuilder();
        //go through each sample
        for (PgxSample sample : samples.values()) {

            //get the filepath, and use the sample id in the name of the file
            stringBuilder.append(SAMPLE_MATRIX_OUT).append(sample.getId()).append(".csv");
            String filePath = stringBuilder.toString();

            //create the directory if it does not exist yet
            File dir = new File(SAMPLE_MATRIX_OUT);
            if (!dir.exists()) dir.mkdirs();

            //setup the writer
            File sampleMatrixFile = new File(filePath);
            FileWriter fileWriter0 = new FileWriter(sampleMatrixFile.getAbsoluteFile());
            BufferedWriter bw0 = new BufferedWriter(fileWriter0);

            //write the header
            bw0.write(header);

            //write the id
            bw0.write(sample.getId());

            //write each allele
            for (PgxGene pgxGene : sample.getGenes().values()) {
                bw0.write("," + pgxGene.getAllele0() + "," + pgxGene.getAllele1());
            }

            //write each phenotype
            for (PgxGene pgxGene : sample.getGenes().values()) {
                bw0.write("," + pgxGene.getPredictedPhenotype());
            }

            bw0.newLine();
            bw0.close();

            //empty the stringbuilder so the next filename can be set
            stringBuilder.setLength(0);
        }

    }

    /*
    following method \/ outputs in the following format

    sample_id,AA12345
    COMT_hap0,COMTHaplotype normal activity
    COMT_hap1,COMTHaplotype low activity
    CYP1A2_hap0,CYP1A2*1A
    CYP1A2_hap1,CYP1A2*1A

     */

    /**
     * write the values for one sample to a matrix, with the value name and value on a new line each
     *
     * @param haplotypeToStarAllele the object containing the pgxgenes and their names
     * @param samples               the list of samples to write
     * @throws IOException
     */
    private static void writeSampleMatrixSplitVertical(HaplotypeToStarAllele haplotypeToStarAllele, Map<String, PgxSample> samples) throws IOException {
        //using stringbuilder is more efficient than keep concatenating strings
        StringBuilder stringBuilder = new StringBuilder();
        //go through each sample
        for (PgxSample sample : samples.values()) {

            //get the filepath, and use the sample id in the name of the file
            stringBuilder.append(SAMPLE_MATRIX_OUT).append(sample.getId()).append(".csv");
            String filePath = stringBuilder.toString();

            //create the directory if it does not exist yet
            File dir = new File(SAMPLE_MATRIX_OUT);
            if (!dir.exists()) {
                dir.mkdirs();
            }

            //setup the writer
            File sampleMatrixFile = new File(filePath);
            FileWriter fileWriter0 = new FileWriter(sampleMatrixFile.getAbsoluteFile());
            BufferedWriter bw0 = new BufferedWriter(fileWriter0);

            //write first line
            bw0.write("sample_id,");
            bw0.write(sample.getId());
            bw0.newLine();

            //write each allele
            for (PgxGene pgxGene : sample.getGenes().values()) {
                //first allele
                bw0.write(pgxGene.getName());
                bw0.write("_hap0,");
                bw0.write(pgxGene.getAllele0());
                bw0.newLine();
                //second allele
                bw0.write(pgxGene.getName());
                bw0.write("_hap1,");
                bw0.write(pgxGene.getAllele1());
                bw0.newLine();
            }

            //write each phenotype
            for (PgxGene pgxGene : sample.getGenes().values()) {
                bw0.write(pgxGene.getPredictedPhenotype());
                bw0.write(pgxGene.getName());
                bw0.write("_phenotype,");
                bw0.write(pgxGene.getPredictedPhenotype());
                bw0.newLine();
            }

            bw0.close();

            //empty the stringbuilder so the next filename can be set
            stringBuilder.setLength(0);
        }

    }

    /**
     * write the sample matrix to one file
     *
     * @param haplotypeToStarAllele the object containing the pgxgenes and their names
     * @param samples               the list of samples to write
     * @throws IOException
     */
    private static void writeSampleMatrixCombined(HaplotypeToStarAllele haplotypeToStarAllele, Map<String, PgxSample> samples) throws IOException {
        //File sampleMatrixFile = new File("/Users/harmbrugge/Documents/PGx/data/richtlijn/sample_matrix.csv");
        //File sampleMatrixFile = new File("C:\\molgenis\\asterix_data\\sample_matrix.csv");
        File sampleMatrixFile = new File(SAMPLE_MATRIX_OUT);
        FileWriter fileWriter0 = new FileWriter(sampleMatrixFile.getAbsoluteFile());
        BufferedWriter bw0 = new BufferedWriter(fileWriter0);

        //get the header
        String header = getSampleMatrixHeader(haplotypeToStarAllele);
        bw0.write(header);

        //go through each sample
        for (PgxSample sample : samples.values()) {

            bw0.write(sample.getId());

            //write each allele
            for (PgxGene pgxGene : sample.getGenes().values()) {
                bw0.write(",");
                bw0.write(pgxGene.getAllele0());
                bw0.write(",");
                bw0.write(pgxGene.getAllele1());
            }

            //write each phenotype
            for (PgxGene pgxGene : sample.getGenes().values()) {
                bw0.write(",");
                if (pgxGene.getPredictedPhenotype() == null) {
                    bw0.write("NA");
                } else bw0.write(pgxGene.getPredictedPhenotype());
            }

            bw0.newLine();
        }

        bw0.close();
    }

    /**
     * get the header for the sample matrix
     *
     * @param haplotypeToStarAllele the object containing the pgxgenes and their names
     * @param sampleColumnName      the name to give the column containing the sample identifier
     * @return the header created for the sample matrix
     */
    private static String getSampleMatrixHeader(HaplotypeToStarAllele haplotypeToStarAllele, String sampleColumnName) {
        //set initial value before null-check
        String nullSafeSampleName = "";
        //overwrite if something else than null was supplied
        if (null != sampleColumnName) {
            nullSafeSampleName = sampleColumnName;
        }
        //use stringbuilder object
        StringBuilder stringBuilder = new StringBuilder();
        stringBuilder.append(sampleColumnName);

        //write haplotype column names
        for (PgxGene pgxGene : haplotypeToStarAllele.getGenes().values()) {
            stringBuilder.append(",").append(pgxGene.getName()).append("_hap0,").append(pgxGene.getName()).append("_hap1");
        }

        //write phenotype column names
        for (PgxGene pgxGene : haplotypeToStarAllele.getGenes().values()) {
            stringBuilder.append(",").append(pgxGene.getName()).append("_phenotype");
        }
        stringBuilder.append("\n");
        return stringBuilder.toString();
    }

    /**
     * get the header for the sample matrix
     *
     * @param haplotypeToStarAllele the object containing the pgxgenes and their names
     * @return the header created for the sample matrix
     */
    private static String getSampleMatrixHeader(HaplotypeToStarAllele haplotypeToStarAllele) {
        return getSampleMatrixHeader(haplotypeToStarAllele, "");
    }

    public Map<String, PgxGene> getGenes() {
        if (null == this.genes) {
            return new TreeMap<String, PgxGene>();
        } else {
            return this.genes;
        }
    }
}
