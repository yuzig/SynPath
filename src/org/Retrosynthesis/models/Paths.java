package org.Retrosynthesis.models;

import java.util.Set;
/**
 * This extracts pathway data from metacyc reaction database
 * @author Y. C. Gao
 */

public class Paths {

    private final String uniqueID;
    private final String name;
    private final Set<String> Species;

    public Paths(String uniqueID, String name, Set<String> species) {
        this.name = name;
        this.Species = species;
        this.uniqueID = uniqueID;
    }

    public Set<String> getSpecies() {
        return Species;
    }
    public Boolean containSpecies(String speciesName){
        if (Species.contains(speciesName)){
            return true;
        }
        return false;
    }
}
