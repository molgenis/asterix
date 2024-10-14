package org.molgenis.asterix.model;

import java.util.Set;

public class SnpLog implements Comparable<SnpLog> {

    private Set<Snp> snps;
    private PgxSample sample;
    private int haplotype;

    public SnpLog(Set<Snp> snps, PgxSample sample, int haplotype) {
        this.snps = snps;
        this.sample = sample;
        this.haplotype = haplotype;
    }

    public Set<Snp> getSnps() {
        return snps;
    }

    public void setSnps(Set<Snp> snps) {
        this.snps = snps;
    }

    public PgxSample getSample() {
        return sample;
    }

    public void setSample(PgxSample sample) {
        this.sample = sample;
    }

    public int getHaplotype() {
        return haplotype;
    }

    public void setHaplotype(int haplotype) {
        this.haplotype = haplotype;
    }

    public double getFirstProbability() {
        return getSnps().iterator().next().getProbability();
    }

    @Override
    public int compareTo(SnpLog o) {
        return Double.compare(o.getFirstProbability(), this.getFirstProbability());
    }

}
