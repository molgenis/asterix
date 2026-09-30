package org.molgenis.asterix.model;

import java.util.Collection;
import java.util.Objects;
import java.util.SortedSet;
import java.util.TreeSet;

/**
 * Identifies one reviewable case: the variants on a single haplotype of a single gene that
 * Asterix could not resolve with certainty.
 *
 * This is the key of the approval (vrijgifte) decision registry, so it has to be stable across
 * runs and round-trip through the Excel file without loss. The variants are held the way the
 * decision logs write them, as sorted tokens in a single string:
 *
 * <ul>
 *   <li>{@code rs3892097_T_P} - the variant allele is present on this haplotype</li>
 *   <li>{@code rs776746_F} - the SNP could not be called reliably (imputed below the cut-off)</li>
 * </ul>
 *
 * For CYP2D6 the CNV status is part of the identity as well: the same variants mean something
 * different on a deleted or duplicated gene, so those cases are reviewed separately.
 *
 * @author BruggeH
 */
public class HaplotypeCombination {

    private static final String VARIANT_SEPARATOR = ";";
    private static final String PRESENT_SUFFIX = "_P";
    private static final String FAILED_SUFFIX = "_F";

    private final String geneName;
    private final SortedSet<String> variants;
    private final String cnvStatus;

    public HaplotypeCombination(String geneName, Collection<String> variants, String cnvStatus) {
        this.geneName = geneName;
        this.variants = new TreeSet<>(variants);
        this.cnvStatus = emptyToNull(cnvStatus);
    }

    /**
     * Rebuild a combination from the strings as they are stored in the decision registry.
     */
    public static HaplotypeCombination fromStrings(String geneName, String variants, String cnvStatus) {
        return new HaplotypeCombination(geneName, splitVariants(variants), cnvStatus);
    }

    /**
     * Describe a variant allele that is present on the haplotype, for example
     * {@code rs3892097_T_P}.
     */
    public static String describePresentVariant(String id, String allele) {
        return id + "_" + allele + PRESENT_SUFFIX;
    }

    /**
     * Describe a SNP that could not be called reliably, for example {@code rs776746_F}. The
     * allele is deliberately left out: what matters to the reviewer is that the call failed, not
     * which way it happened to lean.
     */
    public static String describeFailedCall(String id) {
        return id + FAILED_SUFFIX;
    }

    public static boolean isFailedCall(String variant) {
        return variant.endsWith(FAILED_SUFFIX);
    }

    /**
     * @return the SNP identifier a variant token refers to, so the per-SNP columns of the workbook
     * can be marked
     */
    public static String snpId(String variant) {
        if (isFailedCall(variant)) {
            return variant.substring(0, variant.length() - FAILED_SUFFIX.length());
        }
        String withoutSuffix = variant.endsWith(PRESENT_SUFFIX)
                ? variant.substring(0, variant.length() - PRESENT_SUFFIX.length())
                : variant;
        int separator = withoutSuffix.lastIndexOf('_');
        return separator < 0 ? withoutSuffix : withoutSuffix.substring(0, separator);
    }

    public String getGeneName() {
        return geneName;
    }

    public SortedSet<String> getVariants() {
        return variants;
    }

    public String getCnvStatus() {
        return cnvStatus;
    }

    public String getVariantsString() {
        StringBuilder builder = new StringBuilder();
        for (String variant : variants) {
            if (builder.length() > 0) builder.append(VARIANT_SEPARATOR);
            builder.append(variant);
        }
        return builder.toString();
    }

    /**
     * @return whether this combination has nothing to review, i.e. a plain wild-type haplotype
     * with no failed calls
     */
    public boolean isEmpty() {
        return variants.isEmpty();
    }

    private static SortedSet<String> splitVariants(String value) {
        SortedSet<String> variants = new TreeSet<>();
        if (value == null) return variants;
        for (String variant : value.split(VARIANT_SEPARATOR)) {
            String trimmed = variant.trim();
            if (!trimmed.isEmpty()) variants.add(trimmed);
        }
        return variants;
    }

    private static String emptyToNull(String value) {
        if (value == null) return null;
        String trimmed = value.trim();
        return trimmed.isEmpty() ? null : trimmed;
    }

    @Override
    public String toString() {
        return geneName + "{" +
                "variants=" + getVariantsString() +
                ", cnvStatus=" + cnvStatus +
                '}';
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        HaplotypeCombination that = (HaplotypeCombination) o;

        if (!geneName.equals(that.geneName)) return false;
        if (!variants.equals(that.variants)) return false;
        return Objects.equals(cnvStatus, that.cnvStatus);
    }

    @Override
    public int hashCode() {
        return Objects.hash(geneName, variants, cnvStatus);
    }
}
