package org.molgenis.asterix.model;

import org.molgenis.genotype.annotation.Annotation;

import java.util.HashMap;
import java.util.Map;
import java.util.Objects;

public class Snp {

    private int chr;
    private int pos;
    private String id;
    private String referenceAllele;
    private String variantAllele;
    private String type;
    private Double rSquared;
    private Double probability;
    private String minorAllele;
    private Double maf;
    private boolean isAvailable;
    private double hwePvalue;
    private Map<String, String> annotationValues;

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

    public Double getMaf() {
        return maf;
    }

    public void setMaf(Double maf) {
        this.maf = maf;
    }

    public boolean isAvailable() {
        return isAvailable;
    }

    public void setAvailable(boolean available) {
        isAvailable = available;
    }

    public Double getProbability() {
        return probability;
    }

    public void setProbability(Double probability) {
        this.probability = probability;
    }

    public Snp copySnp(Snp snpToCopy) {
        Snp snp = new Snp();
        snp.setVariantAllele(snpToCopy.getVariantAllele());
        snp.setReferenceAllele(snpToCopy.getReferenceAllele());
        snp.setPos(snpToCopy.getPos());
        snp.setChr(snpToCopy.getChr());
        snp.setId(snpToCopy.getId());
        snp.setMinorAllele(snpToCopy.getMinorAllele());
        snp.setType(snpToCopy.getType());
        snp.setMaf(snpToCopy.getMaf());
        snp.setrSquared(snpToCopy.getrSquared());
        snp.setAvailable(snpToCopy.isAvailable());
        snp.setAnnotations(snpToCopy.getAnnotations());

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

/*    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Snp snp = (Snp) o;
        return Objects.equals(id, snp.id) &&
                Objects.equals(variantAllele, snp.variantAllele);
    }*/

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        Snp snp = (Snp) o;

        if (chr != snp.chr) return false;
        if (pos != snp.pos) return false;
        return variantAllele.equals(snp.variantAllele);
    }

    @Override
    public int hashCode() {
        return Objects.hash(id, variantAllele);
    }

    public void setHwe(double hwePvalue) {
        this.hwePvalue = hwePvalue;
    }

    public void setAnnotations(Map<String, String> annotationValues) {
        this.annotationValues = annotationValues;
    }

    public Map<String, String> getAnnotations() {
        if (annotationValues == null) {
            return new HashMap<>();
        }
        return annotationValues;
    }

    public double getHwe() {
        return hwePvalue;
    }
}
