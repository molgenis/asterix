package org.molgenis.asterix.io;

import org.molgenis.asterix.model.PgxGene;
import org.molgenis.asterix.model.PgxSample;

import java.io.IOException;
import java.util.Collection;
import java.util.Map;

/**
 *
 */
public interface HaplotypeReader {

    void readHaplotypes() throws IOException;

    void readHaplotypes(Collection<PgxGene> genes) throws IOException;

    Map<String, PgxSample> getSamples();
}
