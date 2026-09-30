package org.molgenis.asterix.io;

/**
 * Column names shared by {@link ApprovalReader} and {@link ApprovalWriter}.
 *
 * The names follow the decision logs the reviewers already keep, so an existing log can be read
 * without converting it. The reader resolves columns by header name rather than by position, so a
 * reviewer may reorder columns or add columns of their own without breaking the round trip.
 *
 * @author BruggeH
 */
public final class ApprovalWorkbook {

    // Identity of the case, written by Asterix and not to be edited by the reviewer
    public static final String COL_COMBO_NUM = "combo_num";
    public static final String COL_GENE = "gene";
    public static final String COL_VARIANTS = "variants";
    public static final String COL_CNV_STATUS = "cnv_status";
    public static final String COL_N_HAPLOTYPES = "n_haplotypes";

    // The call Asterix made
    public static final String COL_CALL_INITIAL = "call_initial";
    public static final String COL_DISCLAIMER_INITIAL = "disclaimer_initial";

    // Filled in by the reviewer. A decision equal to call_initial approves it unchanged, a
    // different one overrides it, and NA means the haplotype cannot be released.
    public static final String COL_DECISION = "decision";
    public static final String COL_DISCLAIMER = "disclaimer";
    public static final String COL_LOGGED_BY = "logged_by";
    public static final String COL_DATE = "date";
    public static final String COL_PHARMACIST_CONSULTATION = "pharmacist_consultation_required";
    public static final String COL_PHARMACIST_NAME = "pharmacist_name";

    /** Columns in the order they are written, before the per-SNP columns. */
    public static final String[] COLUMNS = {
            COL_COMBO_NUM,
            COL_GENE,
            COL_VARIANTS,
            COL_CNV_STATUS,
            COL_N_HAPLOTYPES,
            COL_CALL_INITIAL,
            COL_DISCLAIMER_INITIAL,
            COL_DECISION,
            COL_DISCLAIMER,
            COL_LOGGED_BY,
            COL_DATE,
            COL_PHARMACIST_CONSULTATION,
            COL_PHARMACIST_NAME
    };

    /** Columns the reviewer fills in; tinted in the workbook so it is obvious what is editable. */
    public static final String[] REVIEWER_COLUMNS = {
            COL_DECISION,
            COL_DISCLAIMER,
            COL_LOGGED_BY,
            COL_DATE,
            COL_PHARMACIST_CONSULTATION,
            COL_PHARMACIST_NAME
    };

    /** Marks a SNP whose variant allele is present on the haplotype. */
    public static final String MARK_PRESENT = "X";
    /** Marks a SNP that could not be called reliably. */
    public static final String MARK_FAILED = "F";
    /** Marks a SNP that plays no role in a combination. */
    public static final String MARK_ABSENT = ".";

    private ApprovalWorkbook() {

    }
}
