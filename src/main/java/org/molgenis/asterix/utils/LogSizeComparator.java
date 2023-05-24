package org.molgenis.asterix.utils;

import org.molgenis.asterix.model.Snp;
import org.molgenis.asterix.model.SnpLog;

import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class LogSizeComparator implements Comparator<Set<Snp>> {

    private final Map<Set<Snp>, List<SnpLog>> map;

    public LogSizeComparator(final Map<Set<Snp>, List<SnpLog>> map) {
        this.map = map;
    }

    @Override
    public int compare(Set<Snp> s1, Set<Snp> s2) {
        List<?> list1 = this.map.get(s1);
        List<?> list2 = this.map.get(s2);
        Integer length1 = list1.size();
        Integer length2 = list2.size();
        return length2.compareTo(length1);
    }
}