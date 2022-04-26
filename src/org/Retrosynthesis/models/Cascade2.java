package org.Retrosynthesis.models;

import java.util.Set;

public class Cascade2 {

    private final Chems product;  //This Cascade represents all routes to this chemical
    private final Set<Rxns> rxnsThatFormPdt;

    public Cascade2(Chems product, Set<Rxns> rxnsThatFormPdt) {
        this.product = product;
        this.rxnsThatFormPdt = rxnsThatFormPdt;
    }

    public Chems getProduct() {
        return product;
    }

    public Set<Rxns> getRxnsThatFormPdt() {
        return rxnsThatFormPdt;
    }

    public void addRxn(Rxns rxn){
        rxnsThatFormPdt.add(rxn);
    }
}
