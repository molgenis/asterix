package org.molgenis.asterix.pipeline;

import org.molgenis.asterix.io.Hl7Writer;

import java.io.IOException;

 /**
 * Class that contains the pipeline that needs to be gone through
 *
 * @author OelenR & BruggeH
 */
public class Pipeline {

    /**
     * Do the actual pipeline
     * @throws IOException
     */
    public void doPipeline() throws IOException {
        HaplotypeToStarAllele haplotypeToStarAllele = this.getStarAlleleFromHaploType();
        haplotypeToStarAllele.writeStarAlleles();
        StarAlleleToPhenotype starAlleleToPhenotype = getPhenotypeFromStarAllele(haplotypeToStarAllele);
        starAlleleToPhenotype.writePhenotypes(haplotypeToStarAllele);
        Hl7Writer hl7Writer = new Hl7Writer();
        hl7Writer.writeJson(haplotypeToStarAllele.getSamples());
    }

    /**
     * Get the star alleles from the haplotypes
     * @return the haplotypes in an object, with the star alleles determined
     * @throws IOException
     */
    private HaplotypeToStarAllele getStarAlleleFromHaploType() throws IOException {
        HaplotypeToStarAllele haplotypeToStarAllele = new HaplotypeToStarAllele();
        haplotypeToStarAllele.determineStarAlleles();
        return haplotypeToStarAllele;
    }

    /**
     * Get the phenotypes determined from the star alleles
     * @param haplotypeToStarAllele the object containing the star alleles
     * @return the object that determined the phenotypes
     * @throws IOException
     */
    private StarAlleleToPhenotype getPhenotypeFromStarAllele(HaplotypeToStarAllele haplotypeToStarAllele) throws IOException {
        StarAlleleToPhenotype starAlleleToPhenotype = new StarAlleleToPhenotype(haplotypeToStarAllele.getGenes());
        starAlleleToPhenotype.determinePhenotypes(haplotypeToStarAllele);
        return starAlleleToPhenotype;
    }
}
