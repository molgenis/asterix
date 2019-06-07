package org.molgenis.asterix;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.SortedMap;
import java.util.TreeMap;

public class HaploTableReader {

    public HaploTableReader(String haplotypeTableDir) {
        this.haplotypeTableDir = haplotypeTableDir;
    }

    private String haplotypeTableDir;
    private SortedMap<String, PgxGene> genes = new TreeMap<>();

    public void readSnpHaploTables() throws IOException {

        File[] snpHaploTables = new File(haplotypeTableDir).listFiles(file -> !file.isHidden());

        for (File snpHaploTable : snpHaploTables) {

            try (BufferedReader br = new BufferedReader(new FileReader(snpHaploTable))) {
                br.readLine();
                String line = br.readLine();

                String[] splitLine = line.split("\t");
                String geneName = splitLine[1];

                PgxGene pgxGene = new PgxGene(geneName);

                while (line != null) {
                    pgxGene.addSnpToHaplotypes(line);
                    line = br.readLine();
                }

                genes.put(pgxGene.getName(), pgxGene);
            }

        }

    }

    public SortedMap<String, PgxGene> getGenes() {
        return genes;
    }
}
