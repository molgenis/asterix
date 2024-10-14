package org.molgenis.asterix.model;

import java.util.HashMap;
import java.util.Map;

public class PgxHaplotype {

    private PgxGene gene;
    private Map<String, String> variantAlleles = new HashMap<>();
    private String name;
    private String function;
    private String callAs;
    private String disclaimer;

    public PgxHaplotype(PgxGene gene, String name) {
        this.name = name;
        this.gene = gene;
    }

    public PgxHaplotype(PgxGene gene, String name, String callAs) {
        this.name = name;
        this.gene = gene;
        this.callAs = callAs;
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

    public String getCorrected() {
        if (getCallAs() != null) {
            return getCallAs();
        }
        return getName();
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

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        PgxHaplotype haplotype = (PgxHaplotype) o;

        if (!gene.equals(haplotype.gene)) return false;
        return name.equals(haplotype.name);
    }

    @Override
    public int hashCode() {
        int result = gene.hashCode();
        result = 31 * result + name.hashCode();
        return result;
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

    public String getCallAs() {
        return callAs;
    }

    public void setDisclaimer(String disclaimer) {
        this.disclaimer = disclaimer;
    }

    public String getDisclaimer() {
        return disclaimer;
    }
}
