package org.Retrosynthesis;

import org.Retrosynthesis.models.Chems;
import org.Retrosynthesis.models.Rxns;
import org.Utils.FileUtils;

import java.util.*;
/**
 * This extracts reaction data from metacyc reaction database
 * @author Y. C. Gao
 */

public class ExtractRxns {
    public void initiate(){

    }
    public List<Rxns> run(String filepath, HashMap<String, Chems> allchems) throws Exception {
        List<Rxns> allRxn = new ArrayList<>();
        String rxndata = FileUtils.readFile(filepath);
        String[] lines = rxndata.trim().split("//");
        rxndata = null;
        String uniqueID = null;
        String ECnum = null;
        Set<Chems> substrates = new HashSet<>();
        Set<Chems> products = new HashSet<>();
        Set<String> pathways = new HashSet<>();

        for (int i = 2; i < lines.length; i++) {
            String aRxn = lines[i];
            String[] aline = aRxn.trim().split("\\r|\\r?\\n");
            String[] tabs;

            for (String str : aline) {
                if (str.startsWith("UNIQUE-ID")) {
                    tabs = str.split(" - ");
                    uniqueID = tabs[1];
                    continue;
                }
                if (str.startsWith("EC-NUMBER")) {
                    tabs = str.split(" - ");
                    ECnum = tabs[1];
                    continue;
                }
                if (str.startsWith("LEFT")) {
                    tabs = str.split(" - ");
                    Chems substrate = allchems.get(tabs[1]);
                    if (substrate == null){
                        substrate = new Chems(tabs[1],null, null);
                    }
                    substrates.add(substrate);
                    continue;
                }
                if (str.startsWith("RIGHT")) {
                    tabs = str.split(" - ");
                    Chems product = allchems.get(tabs[1]);
                    if (product == null){
                        product = new Chems(tabs[1],null, null);
                    }
                    products.add(product);
                    continue;
                }
                if (str.startsWith("IN-PATHWAY")){
                    tabs = str.split(" - ");
                    pathways.add(tabs[1]);
                }
            }
            Rxns rxn = new Rxns(ECnum, substrates, products, uniqueID,pathways);
            allRxn.add(rxn);
            substrates = new HashSet<>();
            products = new HashSet<>();
            pathways = new HashSet<>();
        }
        return allRxn;
    }

}

