package org.molgenis.asterix.pipeline;

import org.molgenis.asterix.config.ConfigConstants;
import org.molgenis.asterix.config.ConfigProvider;
import org.molgenis.asterix.model.PgxGene;
import org.molgenis.asterix.model.PgxHaplotype;
import org.molgenis.asterix.model.PgxSample;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.Map;

public class Cyp2d6Caller {

    private String CYP2D6_CNV_STATUS_FILE;
    private PgxGene cyp2d6;
    private Map<String, PgxSample> samples;


    public Cyp2d6Caller(Map<String, PgxSample> samples, PgxGene cyp2d6) throws IOException {
        this.loadConfig();
        this.cyp2d6 = cyp2d6;
        this.samples = samples;
        this.readCnvStatusFile();
    }

    private void readCnvStatusFile() throws IOException {
        cyp2d6.getPgxHaplotypes().put("CYP2D6*5", new PgxHaplotype(cyp2d6, "CYP2D6*5"));

        File cyp2d6File = new File(CYP2D6_CNV_STATUS_FILE);

        if (!cyp2d6File.exists()) {
            throw new IllegalArgumentException("CYP2D6 cnv status file not found: " + cyp2d6File.getAbsolutePath());
        }

        try (BufferedReader br = new BufferedReader(new FileReader(cyp2d6File))) {
            br.readLine();

            int i = 0;
            String line = br.readLine();
            double[] probabilities = new double[4];
            while (line != null) {
                String[] splitLine = line.split("\\s");
                double probability = Double.parseDouble(splitLine[2]);
                probabilities[i] = probability;
                line = br.readLine();
                i++;
                if (i % 4 == 0) {
                    i = 0;
                    String sampleId = splitLine[0];
                    procesSample(sampleId, probabilities);
                }
            }
        }
    }

    private void procesSample(String sampleId, double[] probabilities) {
        PgxSample sample = samples.get(sampleId);
        String cnvStatus = cnvStatusFromProbabilities(sample, probabilities);
    }

    private String cnvStatusFromProbabilities(PgxSample sample, double[] probabilities) {
        PgxGene cyp2d6 = sample.getGenes().get("CYP2D6");

        String cnvStatus = null;
        for (int i = 0; i < probabilities.length; i++) {
            double prob = probabilities[i];
            switch (i) {
                case 0:
                    if (prob > 0.97) {
                        cnvStatus = "HOMOZYGOUS_DELETION";
                        cyp2d6.setAllele0("CYP2D6*5");
                        cyp2d6.setAllele1("CYP2D6*5");
                    }
                    break;
                case 1:
                    if (prob > 0.97) {
                        cnvStatus = "HETEROZYGOUS_DELETION";
                        cyp2d6.setAllele0("CYP2D6*5");
                    }
                    break;
                case 2:
                    if (prob > 0.97) cnvStatus = "NO_CNV";
                    break;
                case 3:
                    if (prob > 0.97) cnvStatus = "DUPLICATION";
                    break;
            }
        }
        if (cnvStatus == null) {
            cnvStatus = "LOW_PROBABILITY";
            cyp2d6.setAllele0("CYP2D6*NA");
            cyp2d6.setAllele1("CYP2D6*NA");
        }
        return cnvStatus;
    }

    private void loadConfig() {
        ConfigProvider configProvider = ConfigProvider.getInstance();
        CYP2D6_CNV_STATUS_FILE = configProvider.getConfigParam(ConfigConstants.CYP2D6_CNV_STATUS_FILE);
    }

}
