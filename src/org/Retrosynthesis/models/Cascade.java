package org.Retrosynthesis.models;

import java.util.Set;

public class Cascade {
    private final Chemical product;  //This Cascade represents all routes to this chemical
    private final Set<Reaction> rxnsThatFormPdt;

    public Cascade(Chemical product, Set<Reaction> rxnsThatFormPdt) {
        this.product = product;
        this.rxnsThatFormPdt = rxnsThatFormPdt;
    }

    public Chemical getProduct() {
        return product;
    }

    public Set<Reaction> getRxnsThatFormPdt() {
        return rxnsThatFormPdt;
    }

    public void addRxn(Reaction rxn){
        rxnsThatFormPdt.add(rxn);
    }
    
    //void addReaction(Reaction rxn) {
     //   rxnsThatFormPdt.add(rxn);
    //}
}
