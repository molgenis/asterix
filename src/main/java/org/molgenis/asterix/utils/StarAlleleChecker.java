package org.molgenis.asterix.utils;

import org.molgenis.asterix.io.SnpToHaploTableReader;
import org.molgenis.asterix.io.VcfHaplotypeReader;

import java.io.IOException;


public class StarAlleleChecker {

    private static String IMPUTED_SNP_DIR = "/groups/umcg-gdio/tmp01/projects/2021001/pgx-pipeline/analyses/impute_pgx_genes/out/20221103/postimpute";
    private static String SNP_HAPLO_TABLE_DIR = "/groups/umcg-fg/tmp01/projects/pgx-passport/data/public/pharmvar-5.2.14_haplotypes";

    public StarAlleleChecker() {

    }

    public void readSnpHaploTables() throws IOException {
        SnpToHaploTableReader snpToHaploTableReader = new SnpToHaploTableReader(SNP_HAPLO_TABLE_DIR);
        snpToHaploTableReader.readSnpHaploTables();

        VcfHaplotypeReader vcfHaplotypeReader = new VcfHaplotypeReader(IMPUTED_SNP_DIR);
        vcfHaplotypeReader.readHaplotypes(snpToHaploTableReader.getGenes().values());
    }

    public void start() throws IOException {
        readSnpHaploTables();
        //populateImputedSnpTable();
        //writeOutputFiles();
    }

    public static void main(String[] args) throws IOException {
        StarAlleleChecker starAlleleChecker = new StarAlleleChecker();
        starAlleleChecker.start();
    }

}

