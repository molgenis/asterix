package org.molgenis.asterix.model.hl7;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class Result {

    private ValueCodeableConcept valueCodeableConcept = new ValueCodeableConcept();

    public Result(Result other) {
        this.valueCodeableConcept = new ValueCodeableConcept(other.valueCodeableConcept);
    }

    public class ValueCodeableConcept {

        private List<Map<String, String>> coding;

        public ValueCodeableConcept() {
        }

        public ValueCodeableConcept(ValueCodeableConcept other) {
            this.coding = new ArrayList<>();
            for (Map<String, String> codingItem : other.coding) {
                this.coding.add(new HashMap<>(codingItem));
            }
        }

        public List<Map<String, String>> getCoding() {
            return coding;
        }

        public void setCoding(List<Map<String, String>> coding) {
            this.coding = coding;
        }
    }

    public ValueCodeableConcept getValueCodeableConcept() {
        return valueCodeableConcept;
    }

    public void setValueCodeableConcept(ValueCodeableConcept valueCodeableConcept) {
        this.valueCodeableConcept = valueCodeableConcept;
    }

}
