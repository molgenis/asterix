package org.molgenis.asterix.utils;

import java.io.*;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class UniqueSnpsPerChromosome {

    private static String SNP_HAPLO_TALBLE_DIR = "/Users/harmbrugge/Documents/PGx/data/snp_to_haplo_with_r2/";
    private static String OUTPUT_DIR = "/Users/harmbrugge/Documents/PGx/data/snps_per_chromosome/";

    private Map<Integer, Set<Integer>> snpsPerChromosome = new HashMap<>();


    public UniqueSnpsPerChromosome() {

    }

    public void readSnpTables() throws IOException {

        File[] snpHaploTables = new File(SNP_HAPLO_TALBLE_DIR).listFiles(file -> !file.isHidden());

        for (File snpHaploTable : snpHaploTables) {


            try (BufferedReader br = new BufferedReader(new FileReader(snpHaploTable))) {
                br.readLine();
                String line = br.readLine();

                String[] splitLine = line.split("\t");
                Integer chromosome = Integer.parseInt(splitLine[1]);

                Set<Integer> snpsForChromosome = snpsPerChromosome.get(chromosome);
                if (snpsForChromosome == null) snpsForChromosome = new HashSet<>();

                while (line != null) {
                    splitLine = line.split("\t");
                    Integer snpPos = Integer.parseInt(splitLine[2]) + 1;

                    if (splitLine[7].equals("null")) {
                        line = br.readLine();
                        continue;
                    }

                    double rSquared = Double.parseDouble(splitLine[7]);

                    if (rSquared > 0.9) snpsForChromosome.add(snpPos);

                    line = br.readLine();
                }

                snpsPerChromosome.put(chromosome, snpsForChromosome);

            }
        }


    }

    public void writeSnpsPerChromosomeFile() throws IOException {

        new File(OUTPUT_DIR).mkdirs();

        for (Map.Entry<Integer, Set<Integer>> entry : snpsPerChromosome.entrySet()) {

            Integer chromosome = entry.getKey();
            Set<Integer> snpPos = entry.getValue();

            File outputFile = new File(OUTPUT_DIR + "chr_" + chromosome + ".txt");
            if (!outputFile.exists()) outputFile.createNewFile();

            FileWriter fileWriter = new FileWriter(outputFile.getAbsoluteFile());
            BufferedWriter bw = new BufferedWriter(fileWriter);

            for (Integer pos : snpPos) {
                bw.write(pos + "\n");
            }

            bw.close();

        }

    }

//    public void start() throws IOException {
//        readSnpTables();
//        writeSnpsPerChromosomeFile();
//    }
//
//    public static void main(String[] args) throws IOException {
//        UniqueSnpsPerChromosome uniqueSnpsPerChromosome = new UniqueSnpsPerChromosome();
//        uniqueSnpsPerChromosome.start();
//    }

}
