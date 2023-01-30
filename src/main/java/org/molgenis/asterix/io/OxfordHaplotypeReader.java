package org.molgenis.asterix.io;

import org.apache.commons.io.FileUtils;
import org.molgenis.asterix.model.PgxGene;
import org.molgenis.asterix.model.PgxSample;
import org.molgenis.asterix.model.Snp;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

/**
 * Implementation of HaplotypeReader for the Oxford haplotype file type (.haps)
 */
public class OxfordHaplotypeReader implements HaplotypeReader {

    private Map<String, PgxSample> samples;
    private String haplotypeFolder;
    private List<String> sampleNames;

    public OxfordHaplotypeReader(String haplotypeFolder) {
        this.haplotypeFolder = haplotypeFolder;
        this.sampleNames = new ArrayList<>();
        this.samples = new HashMap<>();
    }

    @Override
    public Map<String, PgxSample> getSamples() {
        return samples;
    }

    @Override
    public void readHaplotypes() throws IOException {
        this.parseSampleNames();

        File[] haplotypeFiles = FileUtils.listFiles(new File(haplotypeFolder), new String[]{"haps"}, true).toArray(new File[0]);

        for (File snpHaploTable : haplotypeFiles) {
            try (BufferedReader br = new BufferedReader(new FileReader(snpHaploTable))) {
                String line = br.readLine();
                while (line != null) {

                    String[] splitLine = line.split("\\s");
                    Snp referenceSnp = new Snp();
                    //TODO Chr is not always available in haps file
//                    referenceSnp.setChr(Integer.parseInt(splitLine[0]));
                    referenceSnp.setId(splitLine[1]);
                    referenceSnp.setPos(Integer.parseInt(splitLine[2]));
                    referenceSnp.setReferenceAllele(splitLine[3]);
                    referenceSnp.setMinorAllele(splitLine[4]);
                    referenceSnp.setVariantAllele(splitLine[3]);

                    Snp minorSnp = referenceSnp.copySnp(referenceSnp);
                    minorSnp.setVariantAllele(splitLine[4]);

                    int sampleCount = 0;

                    for (int i = 5; i < splitLine.length; i = i + 2) {
                        String currentSampleId = sampleNames.get(sampleCount);
                        PgxSample currentSample = samples.get(currentSampleId);

                        int haplotype0 = Integer.parseInt(splitLine[i]);
                        int haplotype1 = Integer.parseInt(splitLine[i + 1]);

                        if (haplotype0 == 0) currentSample.getHaplotype0().put(referenceSnp.getId(), referenceSnp);
                        else currentSample.getHaplotype0().put(minorSnp.getId(), minorSnp);

                        if (haplotype1 == 0) currentSample.getHaplotype1().put(referenceSnp.getId(), referenceSnp);
                        else currentSample.getHaplotype1().put(minorSnp.getId(), minorSnp);

                        sampleCount++;

                    }

                    line = br.readLine();
                }
            }
        }
    }

    @Override
    public void readHaplotypes(Collection<PgxGene> genes) throws IOException {

    }

    /**
     * Parses the sample names from the first .sample file in the haplotype
     * input directory.
     *
     * .haps files have separate sample files
     *
     * @throws IOException when the sample file is missing or is not readable
     */
    private void parseSampleNames() throws IOException {
        File sampleFile = FileUtils.listFiles(new File(haplotypeFolder), new String[]{"sample"}, true).toArray(new File[0])[0];

        try (BufferedReader br = new BufferedReader(new FileReader(sampleFile))) {
            br.readLine();
            br.readLine();
            String line = br.readLine();

            while (line != null) {
                String[] splitLine = line.split("\\s");
                String sampleName = splitLine[1];

                sampleNames.add(sampleName);
                samples.put(sampleName, new PgxSample(sampleName));

                line = br.readLine();
            }
        }
    }
}
