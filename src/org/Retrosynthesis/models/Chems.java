package org.Retrosynthesis.models;

/**
 * A model for chemicals
 * @author Y. C. Gao
 */
public class Chems {

    private final String inchi;
    private final String name;
    private final String id;
    private final String Smiles;

    public Chems(String id, String inchi, String name, String Smiles) {
        this.inchi = inchi;
        this.name = name;
        this.id = id;
        this.Smiles = Smiles;
    }


    public String getInchi() {
        return inchi;
    }

    public String getName() {
        return name;
    }

    public String getID() {
        return id;
    }

    public String getSmiles() {
        return Smiles;
    }
}

