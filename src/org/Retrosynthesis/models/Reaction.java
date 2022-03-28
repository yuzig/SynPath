package org.Retrosynthesis.models;

import java.util.Set;

public class Reaction {
    private final String ecnum;
    private final Set<Chemical> substrates;
    private final Set<Chemical> products;

    public Reaction(String ecnum, Set<Chemical> substrates, Set<Chemical> products) {
        this.ecnum = ecnum;
        this.substrates = substrates;
        this.products = products;
    }

    public String getecnum() {
        return ecnum;
    }

    public Set<Chemical> getSubstrates() {
        return substrates;
    }

    public Set<Chemical> getProducts() {
        return products;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append(ecnum).append("\t");

        for (Chemical c : substrates) {
            sb.append(c.getId()).append("\t");
        }
        sb.append("-->");
        for (Chemical c : products) {
            sb.append(c.getId()).append("\t");
        }
        String out = sb.toString();
        return out;
    }
}
