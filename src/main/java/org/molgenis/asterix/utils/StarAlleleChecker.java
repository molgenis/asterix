package org.molgenis.asterix.utils;

import org.molgenis.asterix.io.SnpToHaploTableReader;
import org.molgenis.asterix.io.VcfHaplotypeReader;

import java.io.IOException;


public class StarAlleleChecker {

    private static String IMPUTED_SNP_DIR = "/groups/umcg-gdio/tmp01/projects/2021001/pgx-pipeline/analyses/impute_pgx_genes/out/20221103/postimpute";
    private static String SNP_HAPLO_TABLE_DIR = "/groups/umcg-fg/tmp01/projects/pgx-passport/data/public/pharmvar-5.2.14_haplotypes";
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

