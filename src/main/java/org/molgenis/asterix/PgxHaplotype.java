package org.molgenis.asterix;

import java.util.HashMap;
import java.util.Map;

public class PgxHaplotype {

    private Map<String, Snp> snps = new HashMap<>();
    private String name;
    private String function;

    public PgxHaplotype(String name) {
        this.name = name;
    }

    public void addSnp(Snp snp) {
        snps.put(snp.getId(), snp);
    }

    public Map<String, Snp> getSnps() {
        return snps;
    }

    public void setSnps(Map<String, Snp> snps) {
        this.snps = snps;
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
                "snps=" + snps +
                ", name='" + name + '\'' +
                ", function='" + function + '\'' +
                '}';
    }
}
