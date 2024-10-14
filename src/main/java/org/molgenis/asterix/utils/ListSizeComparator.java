package org.molgenis.asterix.utils;

import org.molgenis.asterix.model.Snp;

import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class ListSizeComparator implements Comparator<Set<Snp>> {

    private final Map<Set<Snp>, List<String>> map;

    public ListSizeComparator(final Map<Set<Snp>, List<String>> map) {
        this.map = map;
    }

    @Override
    public int compare(Set<Snp> s1, Set<Snp> s2) {
        List<String> list1 = this.map.get(s1);
        List<String> list2 = this.map.get(s2);
        Integer length1 = list1.size();
        Integer length2 = list2.size();
        return length2.compareTo(length1);
    }
}