package org.molgenis.asterix.utils;

import org.apache.commons.io.FileUtils;
import org.molgenis.asterix.model.PgxGene;
import org.molgenis.asterix.model.PgxHaplotype;
import org.molgenis.asterix.model.Snp;

import java.io.*;
import java.util.HashMap;
import java.util.Map;


public class StarAlleleChecker {

    private static String IMPUTED_SNP_DIR = "/Users/harmbrugge/Documents/PGx/LabGeneeskunde/imputatie/";
    private static String SNP_HAPLO_TALBLE_DIR = "/Users/harmbrugge/Documents/PGx/snp_to_haplo_processed/";

    private Map<String, PgxGene> genes = new HashMap<>();
    private Map<String, Double> rSquaredForSnpLocation = new HashMap<>();

    public StarAlleleChecker() {

    }

    private Double getRsquared(String line) {

        String[] splitLine = line.split("\t");

        int position = Integer.parseInt(splitLine[4]) + 1; // index offset

        String positionString = splitLine[3] + ":" + position;
        return rSquaredForSnpLocation.get(positionString);

    }

    public void readSnpHaploTables() throws IOException {

        File[] snpHaploTables = new File(SNP_HAPLO_TALBLE_DIR).listFiles(file -> !file.isHidden());

        for (File snpHaploTable : snpHaploTables) {

            try (BufferedReader br = new BufferedReader(new FileReader(snpHaploTable))) {
                br.readLine();
                String line = br.readLine();

                String[] splitLine = line.split("\t");
                String geneName = splitLine[1];

                PgxGene pgxGene = new PgxGene(geneName);

                while (line != null) {
                    Double rSquared = getRsquared(line);

//                    pgxGene.addSnpToHaplotypes(line, rSquared);
                    line = br.readLine();
                }

                genes.put(pgxGene.getName(), pgxGene);
            }

        }

    }

    public void writeOutputFiles() throws IOException {

        new File(SNP_HAPLO_TALBLE_DIR + "/output").mkdirs();

        for (PgxGene pgxGene : genes.values()) {

            File outputTableFile = new File(SNP_HAPLO_TALBLE_DIR + "/output/" + pgxGene.getName() + "_alleles_and_snps.txt");
            if (!outputTableFile.exists()) outputTableFile.createNewFile();

            FileWriter fileWriter = new FileWriter(outputTableFile.getAbsoluteFile());
            BufferedWriter bw = new BufferedWriter(fileWriter);

            bw.write("Haplotype.Name\tChr\tPos\tid\tReference.Allele\tVariant.Allele\tType\trSquared\n");

            for (PgxHaplotype pgxHaplotype : pgxGene.getPgxHaplotypes().values()) {
                for (Snp snp :pgxHaplotype.getSnps().values()) {
                    bw.write(pgxHaplotype.getName() + "\t" + snp.getChr() + "\t" + snp.getPos() + "\t" + snp.getId() + "\t" + snp.getReferenceAllele() +
                            "\t" + snp.getVariantAllele() + "\t" + snp.getType() + "\t" + snp.getrSquared() + "\n");
                }
            }

            bw.close();

        }

    }

    public void populateImputedSnpTable() throws IOException {

        File[] snpFiles = FileUtils.listFiles(new File(IMPUTED_SNP_DIR), new String[]{"info"}, true).toArray(new File[0]);

        for (File snpFile: snpFiles) {

            System.out.println("Reading file: " + snpFile.getName());

            try (BufferedReader br = new BufferedReader(new FileReader(snpFile))) {
                br.readLine();
                String line = br.readLine();

                while (line != null) {

                    String[] splitLine = line.split("\t");

                    double rSquared = Double.parseDouble(splitLine[6]);
                    String position = splitLine[0];

//                    if (rSquared > 0.97) {
                        rSquaredForSnpLocation.put(position, rSquared);
//                    }

//                    String[] splitPosition = position.split(":");
//
//                    int chr = Integer.parseInt(splitPosition[0]);
//                    int pos = Integer.parseInt(splitPosition[1]);

                    line = br.readLine();
                }
            }
        }


    }

//    public void start() throws IOException {
//        populateImputedSnpTable();
//        readSnpHaploTables();
//        writeOutputFiles();
//    }
//
//    public static void main(String[] args) throws IOException {
//        StarAlleleChecker starAlleleChecker = new StarAlleleChecker();
//        starAlleleChecker.start();
//    }

}

