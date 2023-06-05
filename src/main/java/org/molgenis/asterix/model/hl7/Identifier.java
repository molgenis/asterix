package org.molgenis.asterix.model.hl7;

import java.util.Map;

public class Identifier {

    private String use;
    private Map<String, String> type;
    private String system;
    private String value;
    private String assigner;

    public Identifier(Identifier other) {
        this.use = other.use;
        this.type = other.type;
        this.system = other.system;
        this.value = other.value;
        this.assigner = other.assigner;
    }

    public String getUse() {
        return use;
    }

    public void setUse(String use) {
        this.use = use;
    }

    public Map<String, String> getType() {
        return type;
    }

    public void setType(Map<String, String> type) {
        this.type = type;
    }

    public String getSystem() {
        return system;
    }

    public void setSystem(String system) {
        this.system = system;
    }

    public String getValue() {
        return value;
    }

    public void setValue(String value) {
        this.value = value;
    }

    public String getAssigner() {
        return assigner;
    }

    public void setAssigner(String assigner) {
        this.assigner = assigner;
    }
}
