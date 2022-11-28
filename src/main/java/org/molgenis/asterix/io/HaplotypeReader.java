package org.molgenis.asterix.io;

import org.molgenis.asterix.model.PgxGene;
import org.molgenis.asterix.model.PgxSample;

import java.io.IOException;
import java.util.Collection;
import java.util.Map;

/**
 *
 * Implementations of this interface read a haplotype file containing phased genotypes
 * and collect them on a Map with PgxSamples as values and identifiers as keys.
 * The PgxSample objects contain the PgxHaplotypes.
 *
 * @author Harm Brugge
 *
 */
public interface HaplotypeReader {

    /**
     * Starts reading the haplotype file.
     * @throws IOException when file is not found or not readable.
     */
    void readHaplotypes() throws IOException;

    /**
     * Starts reading the haplotype file.
     * @param genes Collection of genes of which the haplotypes need to be extracted.
     * @throws IOException when file is not found or not readable.
     */
    void readHaplotypes(Collection<PgxGene> genes) throws IOException;

    /**
     *
     * @return the samples
     */
    Map<String, PgxSample> getSamples();
}
