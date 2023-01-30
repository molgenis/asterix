package org.molgenis.asterix.model;

import java.util.HashMap;
import java.util.Map;
import java.util.SortedMap;
import java.util.TreeMap;

public class PgxSample {

    private Map<String, Snp> haplotype0;
    private Map<String, Snp> haplotype1;
    private String id;

    private SortedMap<String, PgxGene> genes = new TreeMap<>();

    public PgxSample(String id) {
        this.id = id;
        haplotype0 = new HashMap<>();
        haplotype1 = new HashMap<>();
    }

    public Map<String, Snp> getHaplotype0() {
        return haplotype0;
    }

    public void setHaplotype0(Map<String, Snp> haplotype0) {
        this.haplotype0 = haplotype0;
    }

    public Map<String, Snp> getHaplotype1() {
        return haplotype1;
    }

    public void setHaplotype1(Map<String, Snp> haplotype1) {
        this.haplotype1 = haplotype1;
    }

    public SortedMap<String, PgxGene> getGenes() {
        return genes;
    }

    public void setGenes(SortedMap<String, PgxGene> genes) {
        this.genes = genes;
    }

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }
}
