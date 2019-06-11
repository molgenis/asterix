package org.molgenis.asterix.pipeline;

import org.molgenis.asterix.HaplotypeToStarAllele;
import org.molgenis.asterix.StarAlleleToPhenotype;

import java.io.IOException;
import java.util.Arrays;
import java.util.List;
import java.util.Map;

/**
 * class that contains the pipeline that needs to be gone through
 *
 * @author OelenR
 */
public class Pipeline {


    /**
     * do the actual pipeline
     * @throws IOException
     */
    public void doPipeline() throws IOException{
        //get the haplotypes
        HaplotypeToStarAllele haplotypeToStarAllele = this.getStarAlleleFromHaploType();
        //write the haplotypes
        haplotypeToStarAllele.writeStarAlleles();
        //get the phenotpyes
        StarAlleleToPhenotype starAlleleToPhenotype = getPhenotypeFromStarAllele(haplotypeToStarAllele);
        //write the phenotypes
        starAlleleToPhenotype.writePhenotypes(haplotypeToStarAllele);
    }

    /**
     * get the star alleles from the haplotypes
     * @return the haplotypes in an object, with the star alleles determined
     * @throws IOException
     */
    private HaplotypeToStarAllele getStarAlleleFromHaploType() throws IOException {
        HaplotypeToStarAllele haplotypeToStarAllele = new HaplotypeToStarAllele();
        haplotypeToStarAllele.determineStarAlleles();
        return haplotypeToStarAllele;
    }

    /**
     * get the phenotypes determined from the star alleles
     * @param haplotypeToStarAllele the object containing the star alleles
     * @return the object that determined the phenotypes
     * @throws IOException
     */
    private StarAlleleToPhenotype getPhenotypeFromStarAllele(HaplotypeToStarAllele haplotypeToStarAllele) throws IOException{
        StarAlleleToPhenotype starAlleleToPhenotype = new StarAlleleToPhenotype(haplotypeToStarAllele.getGenes());
        starAlleleToPhenotype.determinePhenotypes(haplotypeToStarAllele);
        return starAlleleToPhenotype;
    }
}
