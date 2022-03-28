package org.Retrosynthesis.models;

import java.util.Set;

public class Enzyme {
    private final String name;
    private final Set<Protein> subunits;

    public Enzyme(String name, Set<Protein> subunits) {
        this.name = name;
        this.subunits = subunits;
    }

    public String getName() {
        return name;
    }

    public Set<Protein> getSubunits() {
        return subunits;
    }
}

