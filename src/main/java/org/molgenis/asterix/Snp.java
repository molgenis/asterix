package org.molgenis.asterix;

import java.util.Objects;

public class Snp {

    private int chr;
    private int pos;
    private String id;
    private String referenceAllele;
    private String variantAllele;
    private String type;
    private Double rSquared;
    private String minorAllele;

    public int getChr() {
        return chr;
    }

    public void setChr(int chr) {
        this.chr = chr;
    }

    public int getPos() {
        return pos;
    }

    public void setPos(int pos) {
        this.pos = pos;
    }

    public String getId() {
        return id;
    }

    public void setId(String id) {
        this.id = id;
    }

    public String getReferenceAllele() {
        return referenceAllele;
    }

    public void setReferenceAllele(String referenceAllele) {
        this.referenceAllele = referenceAllele;
    }

    public String getVariantAllele() {
        return variantAllele;
    }

    public void setVariantAllele(String variantAllele) {
        this.variantAllele = variantAllele;
    }

    public String getType() {
        return type;
    }

    public void setType(String type) {
        this.type = type;
    }

    public Double getrSquared() {
        return rSquared;
    }

    public void setrSquared(Double rSquared) {
        this.rSquared = rSquared;
    }

    public String getMinorAllele() {
        return minorAllele;
    }

    public void setMinorAllele(String minorAllele) {
        this.minorAllele = minorAllele;
    }

    public Snp copySnp(Snp snpToCopy) {
        Snp snp = new Snp();
        snp.setVariantAllele(snpToCopy.getVariantAllele());
        snp.setReferenceAllele(snpToCopy.getVariantAllele());
        snp.setPos(snpToCopy.getPos());
        snp.setChr(snpToCopy.getChr());
        snp.setId(snpToCopy.getId());
        snp.setMinorAllele(snpToCopy.getMinorAllele());
        snp.setType(snpToCopy.getType());

        return snp;
    }

    @Override
    public String toString() {
        return "Snp{" +
                "chr=" + chr +
                ", pos=" + pos +
                ", id='" + id + '\'' +
                ", referenceAllele='" + referenceAllele + '\'' +
                ", variantAllele='" + variantAllele + '\'' +
                ", type='" + type + '\'' +
                ", rSquared=" + rSquared +
                '}';
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Snp snp = (Snp) o;
        return Objects.equals(id, snp.id) &&
                Objects.equals(variantAllele, snp.variantAllele);
    }

    @Override
    public int hashCode() {
        return Objects.hash(id, variantAllele);
    }
}
