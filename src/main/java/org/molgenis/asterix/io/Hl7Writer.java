package org.molgenis.asterix.io;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import org.molgenis.asterix.config.ConfigConstants;
import org.molgenis.asterix.config.ConfigProvider;
import org.molgenis.asterix.model.PgxGene;
import org.molgenis.asterix.model.PgxSample;
import org.molgenis.asterix.model.hl7.DiagnosticReport;
import org.molgenis.asterix.model.hl7.Result;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class Hl7Writer {

    private static String HL7_OUTPUT_FILE;
    private static String HL7_INPUT_FILE;

    private static String DATE_FORMAT = "yyyy-MM-dd'T'HH:mm:ssXXX";

    private DiagnosticReport dummyDiagnosticReport;
    private Result dummyResult;

    public Hl7Writer() throws FileNotFoundException {
        loadConfig();
        readDummyHl7Report();
    }

    public void writeJson(Map<String, PgxSample> samples) throws IOException {
        FileWriter fileWriter = new FileWriter(HL7_OUTPUT_FILE);
        Gson gson = new GsonBuilder().setPrettyPrinting().disableHtmlEscaping()
                .setDateFormat(DATE_FORMAT).create();

        List<DiagnosticReport> diagnosticReports = new ArrayList<>();

        for (PgxSample sample : samples.values()) {
            DiagnosticReport diagnosticReport = new DiagnosticReport(dummyDiagnosticReport);
            diagnosticReport.getIdentifier().get(0).setValue("report_" + sample.getId());
            diagnosticReport.getSubject().getIdentifier().get(0).setValue(sample.getId());

            List<Result> results = new ArrayList<>();
            for (PgxGene pgxGene : sample.getGenes().values()) {
                if (pgxGene.getContraindication() == null) continue;
                Result result = new Result(dummyResult);
                Map<String, String> codingClinvar = result.getValueCodeableConcept().getCoding().get(0);
                codingClinvar.put("code", pgxGene.getDiplotype());
                codingClinvar.put("display", pgxGene.getDiplotype());
                codingClinvar.put("disclaimer_general", pgxGene.getDisclaimerGeneral());
                codingClinvar.put("disclaimer_subject", pgxGene.getDisclaimerSubject());

                Map<String, String> codingContraindication = result.getValueCodeableConcept().getCoding().get(1);
                codingContraindication.put("code", pgxGene.getContraindication().replaceAll("\\s+", "_"));
                codingContraindication.put("display", pgxGene.getContraindication());

                results.add(result);
            }
            diagnosticReport.setResult(results);
            diagnosticReports.add(diagnosticReport);

        }
        gson.toJson(diagnosticReports, fileWriter);
        fileWriter.close();
    }

    private void readDummyHl7Report() throws FileNotFoundException {
        Gson gson = new Gson();
        this.dummyDiagnosticReport = gson.fromJson(new FileReader(HL7_INPUT_FILE), DiagnosticReport.class);
        this.dummyResult = dummyDiagnosticReport.getResult().get(0);
    }

    private void loadConfig() {
        ConfigProvider configProvider = ConfigProvider.getInstance();
        HL7_INPUT_FILE = configProvider.getConfigParam(ConfigConstants.HL7_INPUT_FILE);
        HL7_OUTPUT_FILE = configProvider.getConfigParam(ConfigConstants.HL7_OUTPUT_FILE);
    }

}
