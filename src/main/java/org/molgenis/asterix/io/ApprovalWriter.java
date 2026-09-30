package org.molgenis.asterix.io;

import org.apache.poi.ss.usermodel.Cell;
import org.apache.poi.ss.usermodel.CellStyle;
import org.apache.poi.ss.usermodel.FillPatternType;
import org.apache.poi.ss.usermodel.Font;
import org.apache.poi.ss.usermodel.HorizontalAlignment;
import org.apache.poi.ss.usermodel.IndexedColors;
import org.apache.poi.ss.usermodel.Row;
import org.apache.poi.ss.usermodel.Sheet;
import org.apache.poi.ss.usermodel.Workbook;
import org.apache.poi.ss.util.WorkbookUtil;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
import org.molgenis.asterix.model.ApprovalDecision;
import org.molgenis.asterix.model.HaplotypeCombination;
import org.molgenis.asterix.model.PgxGene;
import org.molgenis.asterix.model.Snp;
import org.molgenis.asterix.pipeline.ApprovalRegistry;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;

/**
 * Writes the approval (vrijgifte) decision registry to an Excel workbook, one tab per gene.
 *
 * The workbook is both the review document and the registry: it holds the decisions that were
 * already made together with the combinations that turned up in this run and still need one. The
 * whole registry is rewritten every run, so pointing the output at the same file that was read as
 * input keeps a single, growing decision file.
 *
 * The layout follows the pruning-decision workbooks the reviewers already work with: the identity
 * of the case first, then the decision columns, then one column per variant of the gene marking
 * which SNPs are involved.
 *
 * @author BruggeH
 */
public class ApprovalWriter {

    private static final int WIDE_COLUMN_CHARS = 60;
    private static final int NORMAL_COLUMN_CHARS = 22;
    private static final int SNP_COLUMN_CHARS = 16;

    private CellStyle headerStyle;
    private CellStyle reviewerHeaderStyle;
    private CellStyle wrappedStyle;
    private CellStyle snpStyle;

    /**
     * Write the registry.
     *
     * @param approvalFile the workbook to write, its parent directories are created when missing
     * @param registry     the registry to write
     * @param genes        the genes as read from the translation tables, used for the per-SNP columns
     * @throws IOException when the file cannot be written
     */
    public void write(File approvalFile, ApprovalRegistry registry, SortedMap<String, PgxGene> genes)
            throws IOException {
        File parent = approvalFile.getAbsoluteFile().getParentFile();
        if (parent != null) parent.mkdirs();

        try (Workbook workbook = new XSSFWorkbook()) {
            createStyles(workbook);

            for (String geneName : geneNames(registry, genes)) {
                List<ApprovalDecision> decisions = registry.getDecisionsForGene(geneName);
                writeSheet(workbook, geneName, decisions, genes.get(geneName));
            }

            try (OutputStream outputStream = new FileOutputStream(approvalFile)) {
                workbook.write(outputStream);
            }
        }

        System.out.println("Wrote " + registry.getDecisions().size() + " approval decision(s) to " +
                approvalFile.getAbsolutePath() + " (" + registry.getUnreviewedCount() + " awaiting review)");
    }

    /**
     * All genes that need a tab: the genes of this run first, then any gene that only exists in
     * previously recorded decisions, so no decision is ever silently dropped.
     */
    private Set<String> geneNames(ApprovalRegistry registry, SortedMap<String, PgxGene> genes) {
        Set<String> geneNames = new LinkedHashSet<>(genes.keySet());
        for (ApprovalDecision decision : registry.getDecisions()) {
            geneNames.add(decision.getCombination().getGeneName());
        }
        return geneNames;
    }

    private void writeSheet(Workbook workbook, String geneName, List<ApprovalDecision> decisions, PgxGene gene) {
        // The gene is read back from the gene column, not from the tab name, so sanitising here is safe
        Sheet sheet = workbook.createSheet(WorkbookUtil.createSafeSheetName(geneName));
        List<Snp> variants = variantsByPosition(gene);

        writeHeader(sheet, variants);

        int rowNum = 1;
        int comboNum = 1;
        for (ApprovalDecision decision : decisions) {
            writeDecision(sheet.createRow(rowNum), decision, comboNum, variants);
            rowNum += 1;
            comboNum += 1;
        }

        sheet.createFreezePane(0, 1);
        setColumnWidths(sheet, variants);
    }

    private void writeHeader(Sheet sheet, List<Snp> variants) {
        Row header = sheet.createRow(0);
        for (int i = 0; i < ApprovalWorkbook.COLUMNS.length; i++) {
            String column = ApprovalWorkbook.COLUMNS[i];
            Cell cell = header.createCell(i);
            cell.setCellValue(column);
            cell.setCellStyle(isReviewerColumn(column) ? reviewerHeaderStyle : headerStyle);
        }
        for (int i = 0; i < variants.size(); i++) {
            Snp variant = variants.get(i);
            Cell cell = header.createCell(ApprovalWorkbook.COLUMNS.length + i);
            cell.setCellValue(variant.getId() + "_" + variant.getVariantAllele());
            cell.setCellStyle(headerStyle);
        }
    }

    private void writeDecision(Row row, ApprovalDecision decision, int comboNum, List<Snp> variants) {
        HaplotypeCombination combination = decision.getCombination();

        row.createCell(0).setCellValue(comboNum);
        setString(row, 1, combination.getGeneName());
        setWrapped(row, 2, combination.getVariantsString());
        setString(row, 3, combination.getCnvStatus());
        row.createCell(4).setCellValue(decision.getObservedHaplotypeCount());
        setString(row, 5, decision.getCallInitial());
        setWrapped(row, 6, decision.getDisclaimerInitial());
        setString(row, 7, decision.getDecision());
        setWrapped(row, 8, decision.getDisclaimer());
        setString(row, 9, decision.getLoggedBy());
        setString(row, 10, decision.getDate());
        setString(row, 11, decision.getPharmacistConsultationRequired());
        setString(row, 12, decision.getPharmacistName());

        Map<String, String> marksBySnpId = new LinkedHashMap<>();
        for (String variant : combination.getVariants()) {
            marksBySnpId.put(HaplotypeCombination.snpId(variant),
                    HaplotypeCombination.isFailedCall(variant)
                            ? ApprovalWorkbook.MARK_FAILED
                            : ApprovalWorkbook.MARK_PRESENT);
        }

        for (int i = 0; i < variants.size(); i++) {
            String mark = marksBySnpId.get(variants.get(i).getId());
            Cell cell = row.createCell(ApprovalWorkbook.COLUMNS.length + i);
            cell.setCellValue(mark == null ? ApprovalWorkbook.MARK_ABSENT : mark);
            cell.setCellStyle(snpStyle);
        }
    }

    private List<Snp> variantsByPosition(PgxGene gene) {
        List<Snp> variants = new ArrayList<>();
        if (gene == null) return variants;

        for (Map.Entry<String, Snp> variant : gene.getVariants().entrySet()) {
            variants.add(variant.getValue());
        }
        variants.sort(Comparator.comparingInt(Snp::getPos).thenComparing(Snp::getId));
        return variants;
    }

    private void setString(Row row, int column, String value) {
        if (value == null) return;
        row.createCell(column).setCellValue(value);
    }

    private void setWrapped(Row row, int column, String value) {
        if (value == null || value.isEmpty()) return;
        Cell cell = row.createCell(column);
        cell.setCellValue(value);
        cell.setCellStyle(wrappedStyle);
    }

    private boolean isReviewerColumn(String column) {
        for (String reviewerColumn : ApprovalWorkbook.REVIEWER_COLUMNS) {
            if (reviewerColumn.equals(column)) return true;
        }
        return false;
    }

    private void setColumnWidths(Sheet sheet, List<Snp> variants) {
        for (int i = 0; i < ApprovalWorkbook.COLUMNS.length; i++) {
            String column = ApprovalWorkbook.COLUMNS[i];
            boolean wide = ApprovalWorkbook.COL_VARIANTS.equals(column)
                    || ApprovalWorkbook.COL_DISCLAIMER_INITIAL.equals(column)
                    || ApprovalWorkbook.COL_DISCLAIMER.equals(column);
            sheet.setColumnWidth(i, (wide ? WIDE_COLUMN_CHARS : NORMAL_COLUMN_CHARS) * 256);
        }
        for (int i = 0; i < variants.size(); i++) {
            sheet.setColumnWidth(ApprovalWorkbook.COLUMNS.length + i, SNP_COLUMN_CHARS * 256);
        }
    }

    private void createStyles(Workbook workbook) {
        Font boldFont = workbook.createFont();
        boldFont.setBold(true);

        headerStyle = workbook.createCellStyle();
        headerStyle.setFont(boldFont);

        // The columns the reviewer fills in are tinted, so it is obvious what is editable
        reviewerHeaderStyle = workbook.createCellStyle();
        reviewerHeaderStyle.setFont(boldFont);
        reviewerHeaderStyle.setFillForegroundColor(IndexedColors.LIGHT_YELLOW.getIndex());
        reviewerHeaderStyle.setFillPattern(FillPatternType.SOLID_FOREGROUND);

        wrappedStyle = workbook.createCellStyle();
        wrappedStyle.setWrapText(true);

        snpStyle = workbook.createCellStyle();
        snpStyle.setAlignment(HorizontalAlignment.CENTER);
    }
}
