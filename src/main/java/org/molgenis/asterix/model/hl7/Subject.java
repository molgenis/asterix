package org.molgenis.asterix.model.hl7;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class Subject {

    private String resourceType;
    private List<Identifier> identifier;
    private List<Map<String, String>> name;
    private String gender;

    public Subject() {
    }

    public Subject(Subject other) {
        this.resourceType = other.resourceType;
        this.identifier = new ArrayList<>();
        for (Identifier identifierOther : other.identifier) this.identifier.add(new Identifier(identifierOther));
        this.name = other.name;
        this.gender = other.gender;
    }

    public String getResourceType() {
        return resourceType;
    }

    public void setResourceType(String resourceType) {
        this.resourceType = resourceType;
    }

    public List<Identifier> getIdentifier() {
        return identifier;
    }

    public void setIdentifier(List<Identifier> identifier) {
        this.identifier = identifier;
    }

    public List<Map<String, String>> getName() {
        return name;
    }

    public void setName(List<Map<String, String>> name) {
        this.name = name;
    }

    public String getGender() {
        return gender;
    }

    public void setGender(String gender) {
        this.gender = gender;
    }
}
