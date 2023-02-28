package org.molgenis.asterix.utils;

import org.molgenis.asterix.io.SnpToHaploTableReader;
import org.molgenis.asterix.io.VcfHaplotypeReader;

import java.io.IOException;


public class StarAlleleChecker {

    private String imputedSnpDir;
    private String snpHaploTableDir;

    public StarAlleleChecker(String imputedSnpDir, String SnpHaploTableDir) {
        this.imputedSnpDir = imputedSnpDir;
        this.snpHaploTableDir = SnpHaploTableDir;
    }

    public void readSnpHaploTables() throws IOException {
        SnpToHaploTableReader snpToHaploTableReader = new SnpToHaploTableReader(snpHaploTableDir);
        snpToHaploTableReader.readSnpHaploTables();

        VcfHaplotypeReader vcfHaplotypeReader = new VcfHaplotypeReader(imputedSnpDir);
        vcfHaplotypeReader.readHaplotypes(snpToHaploTableReader.getGenes().values());
    }

    public void start() throws IOException {
        readSnpHaploTables();
        //populateImputedSnpTable();
        //writeOutputFiles();
    }

    public static void main(String[] args) throws IOException {
        StarAlleleChecker starAlleleChecker = new StarAlleleChecker(args[0], args[1]);
        starAlleleChecker.start();
    }

}

