package org.Retrosynthesis.models;

import java.util.Set;
/**
 * model for cascades
 * includes a product, and all reactions that form the product
 * @author Y. C. Gao
 */

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
