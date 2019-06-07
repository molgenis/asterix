package org.molgenis.asterix;

import org.apache.commons.io.FileUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class HaplotypeReader {

    private Map<String, Sample> samples;
    private String haplotypeFolder;
    private List<String> sampleNames;

    public HaplotypeReader(String haplotypeFolder) {
        this.haplotypeFolder = haplotypeFolder;
        this.sampleNames = new ArrayList<>();
        this.samples = new HashMap<>();
    }

    public Map<String, Sample> getSamples() {
        return samples;
    }

    public void readHaplotypes() throws IOException {

        File sampleFile = FileUtils.listFiles(new File(haplotypeFolder), new String[]{"sample"}, true).toArray(new File[0])[0];

        try (BufferedReader br = new BufferedReader(new FileReader(sampleFile))) {
            br.readLine();
            br.readLine();
            String line = br.readLine();


            while (line != null) {
                String[] splitLine = line.split("\\s");
                String sampleName = splitLine[1];

                sampleNames.add(sampleName);
                samples.put(sampleName, new Sample(sampleName));

                line = br.readLine();
            }
        }

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
                        Sample currentSample = samples.get(currentSampleId);

                        int haplotype0 = Integer.parseInt(splitLine[i]);
                        int haplotype1 = Integer.parseInt(splitLine[i+1]);

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

}
