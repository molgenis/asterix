package org.molgenis.asterix.model;

import java.util.Set;

public class PGxDiplotype {

    private Set<PgxHaplotype> haplotypes;

    private String predictedPhenotype;
    private String Contraindication;

    public Set<PgxHaplotype> getHaplotypes() {
        return haplotypes;
    }

    public void setHaplotypes(Set<PgxHaplotype> haplotypes) {
        this.haplotypes = haplotypes;
    }

    public String getPredictedPhenotype() {
        return predictedPhenotype;
    }

    public void setPredictedPhenotype(String predictedPhenotype) {
        this.predictedPhenotype = predictedPhenotype;
    }

    public String getContraindication() {
        return Contraindication;
    }

    public void setContraindication(String contraindication) {
        Contraindication = contraindication;
    }

}
