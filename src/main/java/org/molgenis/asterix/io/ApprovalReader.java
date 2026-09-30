package org.molgenis.asterix.io;

import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.CellType;
import org.apache.poi.ss.usermodel.DataFormatter;
import org.apache.poi.ss.usermodel.DateUtil;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.ss.usermodel.WorkbookFactory;
import org.molgenis.asterix.model.ApprovalDecision;
import org.molgenis.asterix.model.HaplotypeCombination;
import org.molgenis.asterix.pipeline.ApprovalRegistry;

import java.io.File;
import java.io.IOException;
import java.text.SimpleDateFormat;
import java.util.HashMap;
import java.util.Map;

/**
 * Reads the approval (vrijgifte) decision registry from an Excel workbook, so decisions a reviewer
 * made in an earlier round are reused instead of being asked again.
 *
 * Every sheet is read, and rows are keyed on the gene and variants columns rather than on the
 * sheet they live in, so an existing single-sheet decision log is read just as well as the
 * one-tab-per-gene workbook Asterix writes.
 *
 * @author BruggeH
 */
public class ApprovalReader {

    private static final String DATE_FORMAT = "yyyy-MM-dd";

    private final DataFormatter dataFormatter = new DataFormatter();

    /**
     * Read all decisions from the workbook into the registry.
     *
     * @param approvalFile the workbook to read
     * @param registry     the registry to fill
     * @throws IOException              when the file cannot be read
     * @throws IllegalArgumentException when the file does not exist
     */
    public void read(File approvalFile, ApprovalRegistry registry) throws IOException {
        if (!approvalFile.exists()) {
            // Failing here is deliberate: silently continuing would drop expert decisions that
            // were already made. Leave the option off for a first run, there is nothing to load yet.
            throw new IllegalArgumentException("Approval decision file not found: " +
                    approvalFile.getAbsolutePath());
        }

        int rowCount = 0;
        try (Workbook workbook = WorkbookFactory.create(approvalFile, null, true)) {
            for (int i = 0; i < workbook.getNumberOfSheets(); i++) {
                rowCount += readSheet(workbook.getSheetAt(i), registry);
            }
        }
        System.out.println("Read " + rowCount + " approval decision(s) from " + approvalFile.getAbsolutePath());
    }

    private int readSheet(Sheet sheet, ApprovalRegistry registry) {
        Row header = sheet.getRow(sheet.getFirstRowNum());
        if (header == null) return 0;

        Map<String, Integer> columns = readHeader(header);
        if (!columns.containsKey(ApprovalWorkbook.COL_GENE)
                || !columns.containsKey(ApprovalWorkbook.COL_VARIANTS)) {
            // Not a decision sheet, for instance a notes tab a reviewer added
            System.out.println("Skipping sheet without approval columns: " + sheet.getSheetName());
            return 0;
        }

        int rowCount = 0;
        for (int i = header.getRowNum() + 1; i <= sheet.getLastRowNum(); i++) {
            Row row = sheet.getRow(i);
            if (row == null) continue;

            String gene = value(row, columns, ApprovalWorkbook.COL_GENE);
            if (gene == null || gene.isEmpty()) continue;

            HaplotypeCombination combination = HaplotypeCombination.fromStrings(
                    gene,
                    value(row, columns, ApprovalWorkbook.COL_VARIANTS),
                    value(row, columns, ApprovalWorkbook.COL_CNV_STATUS));

            ApprovalDecision decision = new ApprovalDecision(combination);
            decision.setCallInitial(value(row, columns, ApprovalWorkbook.COL_CALL_INITIAL));
            decision.setDisclaimerInitial(value(row, columns, ApprovalWorkbook.COL_DISCLAIMER_INITIAL));
            decision.setDecision(value(row, columns, ApprovalWorkbook.COL_DECISION));
            decision.setDisclaimer(value(row, columns, ApprovalWorkbook.COL_DISCLAIMER));
            decision.setLoggedBy(value(row, columns, ApprovalWorkbook.COL_LOGGED_BY));
            decision.setDate(value(row, columns, ApprovalWorkbook.COL_DATE));
            decision.setPharmacistConsultationRequired(
                    value(row, columns, ApprovalWorkbook.COL_PHARMACIST_CONSULTATION));
            decision.setPharmacistName(value(row, columns, ApprovalWorkbook.COL_PHARMACIST_NAME));

            registry.addExistingDecision(decision);
            rowCount += 1;
        }
        return rowCount;
    }

    private Map<String, Integer> readHeader(Row header) {
        Map<String, Integer> columns = new HashMap<>();
        for (int i = header.getFirstCellNum(); i < header.getLastCellNum(); i++) {
            String name = cellValue(header.getCell(i));
            if (name != null && !name.isEmpty()) columns.put(name.trim().toLowerCase(), i);
        }
        return columns;
    }

    private String value(Row row, Map<String, Integer> columns, String column) {
        Integer index = columns.get(column);
        if (index == null) return null;
        return cellValue(row.getCell(index));
    }

    private String cellValue(Cell cell) {
        if (cell == null) return null;
        // Existing logs keep the date as a real Excel date, which formats as a serial number
        // unless it is recognised as one
        if (cell.getCellType() == CellType.NUMERIC && DateUtil.isCellDateFormatted(cell)) {
            return new SimpleDateFormat(DATE_FORMAT).format(cell.getDateCellValue());
        }
        return dataFormatter.formatCellValue(cell).trim();
    }
}
