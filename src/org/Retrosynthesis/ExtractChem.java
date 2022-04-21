package org.Retrosynthesis;

import org.Retrosynthesis.models.Chems;
import org.Utils.FileUtils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class ExtractChem {

    public List<Chems> run(String chempath) throws Exception {
        //Read in all the chemicals
        List<Chems> allChemicals = new ArrayList<>();
        String chemdata = FileUtils.readFile(chempath);
        String[] lines = chemdata.trim().split("//");
        chemdata = null;
        String commonName = null;
        String uniqueID = null;
        String Inchi = null;

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
            }

            Chems achem = new Chems(uniqueID, Inchi,commonName);
            allChemicals.add(achem);

        }
        System.out.println("done populating chemicals");
        return allChemicals;
    }

}
