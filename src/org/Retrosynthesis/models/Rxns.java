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
    private final Double gibbs;

    public Rxns(String ecnum, Set<Chems> substrates, Set<Chems> products, String name, Set<String> pathways, Double Gibbs) {
        this.ecnum = ecnum;
        this.substrates = substrates;
        this.products = products;
        this.name = name;
        this.pathways = pathways;
        this.gibbs = Gibbs;
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

    public Set<String> getPathways() {
        return pathways;
    }

    public String getName() {
        return name;
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (Chems sub : substrates) {
            sb.append(sub.getID());
            sb.append(" ");
        }
        sb.append(" --> ");
        for (Chems pro : products) {
            sb.append(pro.getID());
            sb.append(" ");
        }
        String out = sb.toString();
        return out;
    }

    public Double getGibbs() {
        return gibbs;
    }
}

