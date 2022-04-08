package org.Retrosynthesis;

import org.Retrosynthesis.models.Cascade;
import org.Retrosynthesis.models.Chemical;
import org.Retrosynthesis.models.Reaction;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class ChemicalToCascade {
    private HashMap<Chemical, Cascade> out;

    public void initiate(){
        out = new HashMap<>();
    }

    public HashMap<Chemical, Cascade> run(List<Reaction> rxn, List<Chemical> allchems) throws Exception {
        for (Chemical c : allchems){
            Cascade newCascade = new Cascade(c, new HashSet<>());
            out.put(c, newCascade);
        }
        for (Reaction r : rxn){
            for (Chemical c : r.getProducts()) { out.get(c).addRxn(r);}
        }
        return out;
    }
}
