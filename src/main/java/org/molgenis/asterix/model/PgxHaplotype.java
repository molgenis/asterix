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

    public Map<String, Snp> getSnps() {
        Map<String, Snp> snps = new HashMap<>();
        for (String id: this.variantAlleles.keySet()) {
            Snp snp = new Snp().copySnp(this.gene.getVariants().get(id));
            snp.setVariantAllele(this.variantAlleles.get(id));
            snps.put(id, snp);
        }
        return snps;
    }

    public Map<String, Snp> getSnpsWithWildType() {
        Map<String, Snp> snps = gene.getWildType().getSnps();
        snps.putAll(this.getSnps());

        return snps;
    }

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
