package org.molgenis.asterix.model;

import java.util.HashMap;
import java.util.Map;

public class PgxHaplotype {

    private PgxGene gene;
    private Map<String, String> variantAlleles = new HashMap<>();
    private String name;
    private String function;

    public PgxHaplotype(PgxGene gene, String name) {
        this.name = name;
        this.gene = gene;
    }

//    public void addSnp(Snp snp) {
//        snps.put(snp.getId(), snps);
//    }

//    public boolean hasSnp(Snp snp) {
//        return snps.containsKey(snp.getId());
//    }

//    public void removeSnp(Snp snp) {
//        snps.remove(snp.getId());
//    }

    public Map<String, Snp> getSnps() {
        Map<String, Snp> snps = new HashMap<>();
        for (String id: this.variantAlleles.keySet()) {
            snps.put(id, this.gene.getVariants().get(id));
        }
        return snps;
    }

//    public void setSnps(Map<String, Snp> snps) {
//        this.snps = snps;
//    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getFunction() {
        return function;
    }

    public void setFunction(String function) {
        this.function = function;
    }

    @Override
    public String toString() {
        return "PgxHaplotype{" +
                "snps=" + this.getSnps() +
                ", name='" + name + '\'' +
                ", function='" + function + '\'' +
                '}';
    }

    public Map<String, String> getVariantAlleles() {
        return variantAlleles;
    }

    public void setVariantAlleles(Map<String, String> variantAlleles) {
        this.variantAlleles = variantAlleles;
    }

    public void addVariantAlleles(String id, String variantAllele) {
        this.variantAlleles.put(id, variantAllele);
    }

    public void removeVariantAllele(Snp pgxSnp) {
        variantAlleles.remove(pgxSnp.getId());
    }

    public boolean hasVariantAllele(String id) {
        return variantAlleles.containsKey(id);
    }
}
