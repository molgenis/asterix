package org.molgenis.asterix;

import java.util.*;

public class PgxGene {

    private String name;
    private SortedMap<String, PgxHaplotype> pgxHaplotypes = new TreeMap<>();

    private Map<List<String>, String> functionToPredictedPhenotype = new HashMap<>();

    private String allele0;
    private String allele1;
    private String predictedPhenotype;

    public PgxGene(String name) {
        this.name = name;
    }

    public void addSnpToHaplotypes(String line, Double rSquared) {

        String[] splitLine = line.split("\t");

        if (splitLine.length < 9) {
            throw new IllegalArgumentException("Invalid number of lines in Haplotype table");
        }

        String haplotypeName = splitLine[0];

        Snp snp = new Snp();
        snp.setId(splitLine[2]);
        snp.setChr(Integer.parseInt(splitLine[3]));
        snp.setPos(Integer.parseInt(splitLine[4]));
        snp.setReferenceAllele(splitLine[6]);

        if (splitLine[7].equals("-")) snp.setVariantAllele(splitLine[6]);
        else snp.setVariantAllele(splitLine[7]);

        snp.setrSquared(rSquared);
        snp.setType(splitLine[8]);

        if (!pgxHaplotypes.containsKey(haplotypeName)) {
            PgxHaplotype pgxHaplotype = new PgxHaplotype(haplotypeName);
            pgxHaplotype.addSnp(snp);
            pgxHaplotypes.put(haplotypeName, pgxHaplotype);
        } else {
            PgxHaplotype pgxHaplotype = pgxHaplotypes.get(haplotypeName);
            pgxHaplotype.addSnp(snp);
        }

    }

    public void addSnpToHaplotypes(String line) {
        addSnpToHaplotypes(line, 0.0);
    }


    public void getDeterminableAlleles() {

    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public SortedMap<String, PgxHaplotype> getPgxHaplotypes() {
        return pgxHaplotypes;
    }

    public void setPgxHaplotypes(SortedMap<String, PgxHaplotype> pgxHaplotypes) {
        this.pgxHaplotypes = pgxHaplotypes;
    }

    public String getAllele0() {
        return allele0;
    }

    public void setAllele0(String allele0) {
        this.allele0 = allele0;
    }

    public String getAllele1() {
        return allele1;
    }

    public void setAllele1(String allele1) {
        this.allele1 = allele1;
    }

    public Map<List<String>, String> getFunctionToPredictedPhenotype() {
        return functionToPredictedPhenotype;
    }

    public void setFunctionToPredictedPhenotype(Map<List<String>, String> functionToPredictedPhenotype) {
        this.functionToPredictedPhenotype = functionToPredictedPhenotype;
    }

    public String getPredictedPhenotype() {
        return predictedPhenotype;
    }

    public void setPredictedPhenotype(String predictedPhenotype) {
        this.predictedPhenotype = predictedPhenotype;
    }

    @Override
    public String toString() {
        return "PgxGene{" +
                "name='" + name + '\'' +
                ", pgxHaplotypes=" + pgxHaplotypes +
                ", functionToPredictedPhenotype=" + functionToPredictedPhenotype +
                ", allele0='" + allele0 + '\'' +
                ", allele1='" + allele1 + '\'' +
                ", predictedPhenotype='" + predictedPhenotype + '\'' +
                '}';
    }
}
