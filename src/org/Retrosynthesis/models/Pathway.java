package org.Retrosynthesis.models;


import java.util.List;


public class Pathway {
    private final List<Reaction> Reactions;

    public Pathway(List<Reaction> Reactions) {
        this.Reactions = Reactions;
    }

//    public void printPathway() {
//        ReactionPrinter RP = new ReactionPrinter();
//        for (Reaction c : Reactions){
//            RP.run(c);
//        }
//    }


    public List<Reaction> getReactions() {
        return Reactions;
    }
    public void addRxn(Reaction rxn){
        this.Reactions.add(rxn);
    }
    public void remove(Reaction rxn){
        this.Reactions.remove(rxn);
    }

    public String toString() {
        StringBuilder sb = new StringBuilder();
        for (Reaction r : Reactions){
            sb.append(r.toString()).append("\n");
        }
        String out = sb.toString();
        return out;
    }

    }
