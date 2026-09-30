package org.molgenis.asterix.model;

/**
 * One row of the approval (vrijgifte) decision registry: the call Asterix made for a
 * {@link HaplotypeCombination}, and what the reviewer decided about it.
 *
 * The columns follow the decision logs the reviewers already keep: {@code call_initial} and
 * {@code disclaimer_initial} are what Asterix produced, {@code decision} and {@code disclaimer}
 * are what the reviewer concluded. There is no separate "approved or modified" column; a decision
 * equal to the initial call is an approval, a different one is an override, and {@code NA} means
 * the haplotype cannot be released.
 *
 * Star alleles are held the way the logs write them, without the gene prefix ({@code *4}, not
 * {@code CYP2D6*4}); {@link #getDecidedHaplotype()} returns the qualified name Asterix uses.
 *
 * @author BruggeH
 */
public class ApprovalDecision {

    /** A decision of {@code NA} means the reviewer could not release a call for this haplotype. */
    public static final String NO_CALL = "NA";

    private static final String NO_CALL_HAPLOTYPE_SUFFIX = "*NA";

    private final HaplotypeCombination combination;

    private String callInitial;
    private String disclaimerInitial;
    private String decision;
    private String disclaimer;
    private String loggedBy;
    private String date;
    private String pharmacistConsultationRequired;
    private String pharmacistName;

    /** Number of sample haplotypes in the current run that hit this combination. Not a decision. */
    private int observedHaplotypeCount = 0;

    public ApprovalDecision(HaplotypeCombination combination) {
        this.combination = combination;
    }

    /**
     * Strip the gene prefix off a haplotype name, so {@code CYP2D6*4} is stored as {@code *4}.
     */
    public static String toBareStarAllele(String geneName, String haplotypeName) {
        if (haplotypeName == null) return null;
        if (haplotypeName.startsWith(geneName)) return haplotypeName.substring(geneName.length());
        return haplotypeName;
    }

    /**
     * @return whether a reviewer has recorded a decision on this row
     */
    public boolean isReviewed() {
        return emptyToNull(decision) != null;
    }

    /**
     * @return whether the reviewer decided this haplotype cannot be released
     */
    public boolean isNoCall() {
        return NO_CALL.equalsIgnoreCase(emptyToNull(decision));
    }

    /**
     * @return whether the reviewer overruled the call Asterix made
     */
    public boolean isModified() {
        if (!isReviewed()) return false;
        return !decision.trim().equals(emptyToNull(callInitial));
    }

    /**
     * @return the haplotype that should be called for this combination, qualified with the gene
     * name, or null when the row has not been reviewed. A no-call decision yields the gene's
     * {@code *NA} haplotype, which has no phenotype and is therefore left out of the report.
     */
    public String getDecidedHaplotype() {
        if (!isReviewed()) return null;
        if (isNoCall()) return combination.getGeneName() + NO_CALL_HAPLOTYPE_SUFFIX;

        String decided = decision.trim();
        if (decided.startsWith(combination.getGeneName())) return decided;
        return combination.getGeneName() + decided;
    }

    /**
     * @return the disclaimer to attach to the call: the reviewer's, falling back to the one
     * Asterix produced when the reviewer left it empty, or null when there is neither
     */
    public String getEffectiveDisclaimer() {
        if (!isReviewed()) return null;
        String reviewed = emptyToNull(disclaimer);
        return reviewed != null ? reviewed : emptyToNull(disclaimerInitial);
    }

    private static String emptyToNull(String value) {
        if (value == null) return null;
        String trimmed = value.trim();
        return trimmed.isEmpty() || ".".equals(trimmed) ? null : trimmed;
    }

    public HaplotypeCombination getCombination() {
        return combination;
    }

    public String getCallInitial() {
        return callInitial;
    }

    public void setCallInitial(String callInitial) {
        this.callInitial = callInitial;
    }

    public String getDisclaimerInitial() {
        return disclaimerInitial;
    }

    public void setDisclaimerInitial(String disclaimerInitial) {
        this.disclaimerInitial = disclaimerInitial;
    }

    public String getDecision() {
        return decision;
    }

    public void setDecision(String decision) {
        this.decision = decision;
    }

    public String getDisclaimer() {
        return disclaimer;
    }

    public void setDisclaimer(String disclaimer) {
        this.disclaimer = disclaimer;
    }

    public String getLoggedBy() {
        return loggedBy;
    }

    public void setLoggedBy(String loggedBy) {
        this.loggedBy = loggedBy;
    }

    public String getDate() {
        return date;
    }

    public void setDate(String date) {
        this.date = date;
    }

    public String getPharmacistConsultationRequired() {
        return pharmacistConsultationRequired;
    }

    public void setPharmacistConsultationRequired(String pharmacistConsultationRequired) {
        this.pharmacistConsultationRequired = pharmacistConsultationRequired;
    }

    public String getPharmacistName() {
        return pharmacistName;
    }

    public void setPharmacistName(String pharmacistName) {
        this.pharmacistName = pharmacistName;
    }

    public int getObservedHaplotypeCount() {
        return observedHaplotypeCount;
    }

    public void incrementObservedHaplotypeCount() {
        this.observedHaplotypeCount += 1;
    }

    @Override
    public String toString() {
        return "ApprovalDecision{" +
                "combination=" + combination +
                ", callInitial='" + callInitial + '\'' +
                ", decision='" + decision + '\'' +
                '}';
    }
}
