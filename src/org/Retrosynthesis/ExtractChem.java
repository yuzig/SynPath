package org.Retrosynthesis;

import org.Retrosynthesis.models.Chems;
import org.Utils.FileUtils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

/**
 * This extracts chemical compounds from metacyc compound database
 * @author Y. C. Gao
 */
public class ExtractChem {
    private HashMap<String, Chems> chemsHashMap;
    private List<String> AAlists;

    public void initiate(){
        chemsHashMap = new HashMap<>();
        AAlists = new ArrayList<>();
    }

    public List<Chems> run(String chempath) throws Exception {
        //Read in all the chemicals
        List<Chems> allChemicals = new ArrayList<>();
        String chemdata = FileUtils.readFile(chempath);
        String[] lines = chemdata.trim().split("//");
        chemdata = null;
        String commonName = null;
        String uniqueID = null;
        String Inchi = null;
        Boolean isAA = false;

        for (int i = 2; i < lines.length; i++) {
            String aCompound = lines[i];
            String[] aline = aCompound.trim().split("\\r|\\r?\\n");
            String[] tabs;

            for (String str : aline){
                if (str.startsWith("INCHI-KEY")){
                    continue;
                }
                if (str.startsWith("UNIQUE-ID")){
                    tabs = str.split(" - ");
                    uniqueID = tabs[1];
                    continue;
                }
                if (str.startsWith("INCHI")){
                    tabs = str.split(" - ");
                    Inchi = tabs[1];
                    continue;
                }
                if (str.startsWith("COMMON-NAME")) {
                    tabs = str.split(" - ");
                    commonName = tabs[1];
                    continue;
                }

                if (str.startsWith("TYPES")) {
                    tabs = str.split(" - ");
                    if (tabs[1].contains("Inorganic-Compounds") || tabs[1].contains("Electorin-Carriers") || tabs[1].contains("Amino-Acids") || tabs[1].contains("Amino-Acid") || tabs[1].contains("Glutathione") || tabs[1].contains("Polyphosphates") || tabs[1].contains("Metal-Cations") || tabs[1].contains("Inorganic-Anions")){
                        isAA = true;
                        continue;
                    }
                }

                if (str.startsWith("COMMON-NAME")) {
                    tabs = str.split(" - ");
                    commonName = tabs[1];
                    continue;
                }
            }

            Chems achem = new Chems(uniqueID, Inchi,commonName);
            if (isAA == true){
                AAlists.add(achem.getID());
            }
            allChemicals.add(achem);
            chemsHashMap.put(uniqueID, achem);
            commonName = null;
            uniqueID = null;
            Inchi = null;
            isAA = false;
        }
        System.out.println("done populating chemicals");
        return allChemicals;
    }
    public HashMap<String, Chems> getChemsHashMap() {
        return chemsHashMap; }

    public List<String> getAAlists() {
        return AAlists;
    }
}
