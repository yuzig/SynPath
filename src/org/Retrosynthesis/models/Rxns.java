package org.Retrosynthesis.models;

import java.util.Set;
/**
 * A model for reactions
 * @author Y. C. Gao
 */
public class Rxns {
    private final String ecnum;
    private final Set<Chems> substrates;
    private final Set<Chems> products;
    private final String name;
    private final Set<String> pathways;

    public Rxns(String ecnum, Set<Chems> substrates, Set<Chems> products, String name, Set<String> pathways) {
        this.ecnum = ecnum;
        this.substrates = substrates;
        this.products = products;
        this.name = name;
        this.pathways = pathways;
    }

    public String getecnum() {
        return ecnum;
    }

    public Set<Chems> getSubstrates() {
        return substrates;
    }

    public Set<Chems> getProducts() {
        return products;
    }

    public boolean inPathway() {
        if (pathways.isEmpty()){
            return false;
        } else {
            return true;
        }
    }
}
