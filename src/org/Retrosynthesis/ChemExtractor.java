package org.Retrosynthesis;

import org.Retrosynthesis.models.Chemical;
import org.Utils.FileUtils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class ChemExtractor {
    private HashMap<String, Chemical> chemicalHashMap;
    public void initiate() {
        chemicalHashMap = new HashMap<>();
    }

    public List<Chemical> run(String chempath) throws Exception {
        //Read in all the chemicals
        List<Chemical> allChemicals = new ArrayList<Chemical>();
        String chemdata = FileUtils.readFile(chempath);
        chemdata = chemdata.replaceAll("\"", "");
        String[] lines = chemdata.trim().split("\\r|\\r?\\n");
        chemdata = null;

        //Each line of the file is a chemical after a header
        for (int i = 1; i < lines.length; i++) {
            String aline = lines[i];
            String[] tabs = aline.trim().split("\t");

            Long id = Long.parseLong(tabs[0]);
            String name = tabs[1];
            String inchi = tabs[2];
            String smiles = tabs[3];

            Chemical achem = new Chemical(id, inchi, smiles, name);
            chemicalHashMap.put(achem.getName(), achem);
            allChemicals.add(achem);

        }
        System.out.println("done populating chemicals");
        return allChemicals;
    }

    public HashMap<String, Chemical> getChemicalHashMap() {
        return chemicalHashMap;
    }
}
