package org.molgenis.asterix.model.hl7;

import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.Map;

public class DiagnosticReport {

    private String resourceType;
    private List<Identifier> identifier;
    private String status;
    private Subject subject;
    private Date effectiveDateTime;
    private Map<String, Date> effectivePeriod;
    private Date issued = new Date(System.currentTimeMillis());
    private List<Result> result;

    public DiagnosticReport() {
    }

    public DiagnosticReport(DiagnosticReport other) {
        this.resourceType = other.resourceType;
        this.identifier = new ArrayList<>();
        for (Identifier identifierOther : other.identifier) this.identifier.add(new Identifier(identifierOther));
        this.status = other.status;
        this.subject = new Subject(other.subject);
        this.effectiveDateTime = other.effectiveDateTime;
        this.effectivePeriod = other.effectivePeriod;
        this.issued = other.issued;
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

    public String getStatus() {
        return status;
    }

    public void setStatus(String status) {
        this.status = status;
    }

    public Subject getSubject() {
        return subject;
    }

    public void setSubject(Subject subject) {
        this.subject = subject;
    }

    public Date getEffectiveDateTime() {
        return effectiveDateTime;
    }

    public void setEffectiveDateTime(Date effectiveDateTime) {
        this.effectiveDateTime = effectiveDateTime;
    }

    public Map<String, Date> getEffectivePeriod() {
        return effectivePeriod;
    }

    public void setEffectivePeriod(Map<String, Date> effectivePeriod) {
        this.effectivePeriod = effectivePeriod;
    }

    public Date getIssued() {
        return issued;
    }

    public void setIssued(Date issued) {
        this.issued = issued;
    }

    public List<Result> getResult() {
        return result;
    }

    public void setResult(List<Result> result) {
        this.result = result;
    }

}
