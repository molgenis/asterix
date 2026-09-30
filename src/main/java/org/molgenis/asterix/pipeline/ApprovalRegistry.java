package org.molgenis.asterix.pipeline;

import org.molgenis.asterix.model.ApprovalDecision;
import org.molgenis.asterix.model.HaplotypeCombination;

import java.util.ArrayList;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

/**
 * The approval (vrijgifte) decision registry: every combination of mutations and unreliable calls
 * that has been seen, together with the decision a reviewer made about it.
 *
 * Rows read from a previously reviewed workbook are kept verbatim, so decisions survive across
 * runs. Combinations that turn up in the current run and are not in the registry yet are added as
 * unreviewed rows, which is what puts them in front of a reviewer in the next round.
 *
 * @author BruggeH
 */
public class ApprovalRegistry {

    // Insertion ordered, so previously reviewed rows keep their place and new cases are appended
    private final Map<HaplotypeCombination, ApprovalDecision> decisions = new LinkedHashMap<>();

    private int newCombinationCount = 0;

    /**
     * Add a row that was read from an existing decision workbook.
     */
    public void addExistingDecision(ApprovalDecision decision) {
        decisions.put(decision.getCombination(), decision);
    }

    /**
     * Look the combination up, adding it as a new unreviewed row when the registry does not know
     * it yet, and count this observation.
     *
     * @param combination the combination observed on a sample haplotype
     * @param callInitial the call Asterix made for it without the gene prefix, recorded on new
     *                    rows so the reviewer can approve it unchanged
     * @param disclaimer  the disclaimer Asterix produced for it, may be null
     * @return the registry row for this combination, never null
     */
    public ApprovalDecision observe(HaplotypeCombination combination, String callInitial, String disclaimer) {
        ApprovalDecision decision = decisions.get(combination);
        if (decision == null) {
            decision = new ApprovalDecision(combination);
            decision.setCallInitial(callInitial);
            decision.setDisclaimerInitial(disclaimer);
            decisions.put(combination, decision);
            newCombinationCount += 1;
        }
        decision.incrementObservedHaplotypeCount();
        return decision;
    }

    public ApprovalDecision get(HaplotypeCombination combination) {
        return decisions.get(combination);
    }

    public List<ApprovalDecision> getDecisions() {
        return new ArrayList<>(decisions.values());
    }

    /**
     * @return the rows belonging to one gene, in registry order
     */
    public List<ApprovalDecision> getDecisionsForGene(String geneName) {
        List<ApprovalDecision> geneDecisions = new ArrayList<>();
        for (ApprovalDecision decision : decisions.values()) {
            if (decision.getCombination().getGeneName().equals(geneName)) {
                geneDecisions.add(decision);
            }
        }
        return geneDecisions;
    }

    /**
     * @return the number of combinations added during this run, i.e. the cases awaiting review
     */
    public int getNewCombinationCount() {
        return newCombinationCount;
    }

    /**
     * @return the number of registry rows that are still waiting for a reviewer
     */
    public int getUnreviewedCount() {
        int unreviewed = 0;
        for (ApprovalDecision decision : decisions.values()) {
            if (!decision.isReviewed()) unreviewed += 1;
        }
        return unreviewed;
    }
}
