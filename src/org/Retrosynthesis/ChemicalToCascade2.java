package org.Retrosynthesis;

import org.Retrosynthesis.models.*;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

public class ChemicalToCascade2 {
    private HashMap<Chems, Cascade2> out;

    public void initiate(){
        out = new HashMap<>();
    }

    public HashMap<Chems, Cascade2> run(List<Rxns> rxn, List<Chems> allchems) throws Exception {
        for (Chems c : allchems){
            if (c.getInchi() == null){
                continue;
            }
            Cascade2 newCascade = new Cascade2(c, new HashSet<>());
            out.put(c, newCascade);
        }
        for (Rxns r : rxn){
            for (Chems c : r.getProducts()) {
                if (c.getInchi() == null){
                    continue;
                }
                out.get(c).addRxn(r);}
        }
        System.out.println("done populating cascademap");
        return out;
    }
}
