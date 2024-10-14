package org.molgenis.asterix.model;

import java.util.Iterator;
import java.util.Set;

public class PGxDiplotype {

    private Set<PgxHaplotype> haplotypes;

    private String predictedPhenotype;
    private String contraindication;

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
        return contraindication;
    }

    public void setContraindication(String contraindication) {
        this.contraindication = contraindication;
    }

    public String getDiplotypeString() {
        StringBuilder diplotypeString = new StringBuilder();
        for (Iterator<PgxHaplotype> iterator = getHaplotypes().iterator(); iterator.hasNext(); ) {
            PgxHaplotype pgxHaplotype = iterator.next();
            diplotypeString.append(pgxHaplotype.getCorrected());
            if (iterator.hasNext()) diplotypeString.append("/");
        }
        if (getHaplotypes().size() == 1) {
            diplotypeString.append("/");
            diplotypeString.append(getHaplotypes().iterator().next().getCorrected());
        }
        return diplotypeString.toString();
    }

}
