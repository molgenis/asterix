package org.molgenis.asterix.io;

import org.molgenis.asterix.model.PgxGene;
import org.molgenis.asterix.model.PgxHaplotype;
import org.molgenis.asterix.model.Snp;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.SortedMap;
import java.util.TreeMap;

import static java.lang.System.exit;

public class SnpToHaploTableReader {

    public SnpToHaploTableReader(String snpToHaplotypeTableDir) {
        this.haplotypeTableDir = snpToHaplotypeTableDir;
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
                String chrString = splitLine[3];
                String[] chrList = chrString.split("\\.");
                if (chrList.length > 1) {
                    String chr = chrList[0].split("0{4,5}")[1];
                    pgxGene.setChr(Integer.parseInt(chr));
                }
                else pgxGene.setChr(Integer.parseInt(chrString));

                while (line != null) {
                    addSnpToHaplotype(pgxGene, line);
                    line = br.readLine();
                }

                // Add wild type at the end because the reference is only known when all snp are processed
                pgxGene.getPgxHaplotypes().put(pgxGene.getWildType().getName(), pgxGene.getWildType());
                genes.put(pgxGene.getName(), pgxGene);
            }
        }
    }

    private void addSnpToHaplotype(PgxGene gene, String line) {
        String[] splitLine = line.split("\t");
        if (splitLine.length != 9) {
            System.out.println(line);
            throw new IllegalArgumentException("Invalid number of lines in translation table");
        }

        String haplotypeName = splitLine[0];

        Snp snp = new Snp();
        snp.setChr(gene.getChr());
        snp.setPos(Integer.parseInt(splitLine[4]));
        snp.setReferenceAllele(splitLine[6]);
        snp.setType(splitLine[8]);

        String id = splitLine[2];
        if (id.equals("-")) snp.setId(snp.getChr() + ":" + snp.getPos() + ":" + snp.getReferenceAllele());
        else snp.setId(id);


        if (snp.getPos() > gene.getEndPos()) gene.setEndPos(snp.getPos());
        if (snp.getPos() < gene.getStartPos()) gene.setStartPos(snp.getPos());

        if (splitLine[7].equals("-")) snp.setVariantAllele(splitLine[6]); // TODO what about indels
        else snp.setVariantAllele(splitLine[7]);

        if (!gene.getWildType().hasSnp(snp)) {
            Snp wildTypeSnp = snp.copySnp(snp);
            wildTypeSnp.setVariantAllele(wildTypeSnp.getReferenceAllele());
            gene.getWildType().addSnp(wildTypeSnp);
        }
        if (gene.getPgxHaplotypes().containsKey(haplotypeName)) {
            PgxHaplotype pgxHaplotype = gene.getPgxHaplotypes().get(haplotypeName);
            pgxHaplotype.addSnp(snp);
        } else {
            PgxHaplotype pgxHaplotype = new PgxHaplotype(haplotypeName);
            pgxHaplotype.addSnp(snp);
            gene.getPgxHaplotypes().put(haplotypeName, pgxHaplotype);
        }
    }

    public SortedMap<String, PgxGene> getGenes() {
        return genes;
    }
}
