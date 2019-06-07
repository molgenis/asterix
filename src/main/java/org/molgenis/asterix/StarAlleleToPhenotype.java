package org.molgenis.asterix;

import org.apache.commons.io.filefilter.WildcardFileFilter;

import java.io.*;
import java.util.*;

public class StarAlleleToPhenotype {

    private static final String CONVERSION_TABLE_DIR = "/Users/harmbrugge/Documents/PGx/data/richtlijn/haplo_to_pheno/";
    private static final String PREDICTED_PHENOTYPES_OUTPUT_DIR = "/Users/harmbrugge/Documents/PGx/data/richtlijn/predicted_phenotypes/";

    private SortedMap<String, PgxGene> genes;

    public StarAlleleToPhenotype(SortedMap<String, PgxGene> genes) {
        this.genes = genes;
    }

    private void readConversionTables() throws IOException {

        File[] alleleFunctionTables = new File(CONVERSION_TABLE_DIR).listFiles((FileFilter) new WildcardFileFilter("*alleles.txt"));
        File[] phenotypeTables = new File(CONVERSION_TABLE_DIR).listFiles((FileFilter) new WildcardFileFilter("*phenotype.txt"));

        for (File alleleFunctionTable : alleleFunctionTables) {
            try (BufferedReader br = new BufferedReader(new FileReader(alleleFunctionTable))) {
                br.readLine();
                String line = br.readLine();
                String[] splitLine = line.split("\t");

                String geneName = splitLine[0];
                PgxGene pgxGene = genes.get(geneName);

                while (line != null) {
                    splitLine = line.split("\t");
                    String starAlleleName = splitLine[1];
                    String function = splitLine[2];

                    PgxHaplotype pgxHaplotype = pgxGene.getPgxHaplotypes().get(starAlleleName);

                    if (pgxHaplotype == null) System.out.println(starAlleleName);

                    pgxHaplotype.setFunction(function);


                    line = br.readLine();
                }

            }
        }

        for (File phenotypeTable : phenotypeTables) {
            try (BufferedReader br = new BufferedReader(new FileReader(phenotypeTable))) {
                String line = br.readLine();
                String[] splitLine = line.split("\t");

                String geneName = splitLine[0];
                PgxGene pgxGene = genes.get(geneName);

                String[] functions = Arrays.copyOfRange(splitLine, 1, splitLine.length);

                //TODO: remove me
                if (pgxGene == null) System.out.println(geneName);

                Map<List<String>, String> functionToPredictedPhenotype = pgxGene.getFunctionToPredictedPhenotype();

                line = br.readLine();
                int i = 0;
                while (line != null) {
                    splitLine = line.split("\t");

                    for (int j = 0; j < splitLine.length - 1 ; j++) {
                        functionToPredictedPhenotype.put(
                                Collections.unmodifiableList(Arrays.asList(functions[i], functions[j])),
                                splitLine[j+1]);
                    }

                    i++;
                    line = br.readLine();
                }
            }
        }

    }

    public static void main(String[] args) throws IOException {
        HaplotypeToStarAllele haplotypeToStarAllele = new HaplotypeToStarAllele();
        haplotypeToStarAllele.determineStarAlleles();
        haplotypeToStarAllele.writeStarAlleles();

        StarAlleleToPhenotype starAlleleToPhenotype = new StarAlleleToPhenotype(haplotypeToStarAllele.getGenes());
        starAlleleToPhenotype.readConversionTables();

        Map<String, Sample> samples = haplotypeToStarAllele.getSamples();

        for (Sample sample : samples.values()) {
            for (PgxGene gene : sample.getGenes().values()) {

                PgxGene geneWithInfo = starAlleleToPhenotype.genes.get(gene.getName());

                //TODO Fix this Alleles shouldn't be null
                if (gene.getAllele0() == null) continue;
                if (gene.getAllele1() == null) continue;

                String functionAllele0 = geneWithInfo.getPgxHaplotypes().get(gene.getAllele0()).getFunction();
                String functionAllele1 = geneWithInfo.getPgxHaplotypes().get(gene.getAllele1()).getFunction();

                Map<List<String>, String> functionToPredictedPhenotype = starAlleleToPhenotype.genes.get(gene.getName()).getFunctionToPredictedPhenotype();
                String predictedPhenotype = functionToPredictedPhenotype.get(Arrays.asList(functionAllele0, functionAllele1));

                gene.setPredictedPhenotype(predictedPhenotype);
            }
        }




        new File(PREDICTED_PHENOTYPES_OUTPUT_DIR).mkdirs();

        File sampleMatrixFile = new File("/Users/harmbrugge/Documents/PGx/data/richtlijn/sample_matrix.csv");
        FileWriter fileWriter0 = new FileWriter(sampleMatrixFile.getAbsoluteFile());
        BufferedWriter bw0 = new BufferedWriter(fileWriter0);

        for (PgxGene pgxGene : haplotypeToStarAllele.getGenes().values()) {
            bw0.write("," + pgxGene.getName() + "_hap0," + pgxGene.getName() + "_hap1");
        }

        for (PgxGene pgxGene : haplotypeToStarAllele.getGenes().values()) {
            bw0.write("," + pgxGene.getName() + "_phenotype");
        }

        bw0.write("\n");

        for (Sample sample : samples.values()) {

            bw0.write(sample.getId());

            for (PgxGene pgxGene : sample.getGenes().values()) {
                bw0.write("," + pgxGene.getAllele0() + "," + pgxGene.getAllele1());
            }
            for (PgxGene pgxGene : sample.getGenes().values()) {
                bw0.write("," + pgxGene.getPredictedPhenotype());
            }

            bw0.write("\n");
        }

        bw0.close();

        for (PgxGene pgxGene : haplotypeToStarAllele.getGenes().values()) {

            File predictedPhenotypeOutputFile = new File(PREDICTED_PHENOTYPES_OUTPUT_DIR + pgxGene.getName() + "_predicted_phenotypes.txt");
            if (!predictedPhenotypeOutputFile.exists()) predictedPhenotypeOutputFile.createNewFile();

            FileWriter fileWriter = new FileWriter(predictedPhenotypeOutputFile.getAbsoluteFile());
            BufferedWriter bw = new BufferedWriter(fileWriter);

            for (Sample sample : samples.values()) {
                String line = sample.getId() + "\t" + sample.getGenes().get(pgxGene.getName()).getPredictedPhenotype() + "\n";
                bw.write(line);
            }

            bw.close();
        }

    }

}
