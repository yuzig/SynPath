package org.Retrosynthesis;

import org.C5;
import org.Retrosynthesis.models.*;

import java.util.*;
/**
 * function for enumerating all pathways leading to the production of target metabolite
 * @author Y. C. Gao
 */

public class PathEnum {

    private HashMap<String, Chems> chems;
    private List<Rxns> reactions;
    private HashMap<String, Paths> PathMap;
    private List<String> natives;
    private HashMap<Chems, Cascade2> chemToCascadeMap;
    private HashMap<String, Rxns> rxnsHashMap;

    public void initiate() throws Exception{
        ExtractPaths ep = new ExtractPaths();
        ExtractChem ec = new ExtractChem();
        ExtractRxns er = new ExtractRxns();
        ExtractCoreMet ecm = new ExtractCoreMet();
        ChemicalToCascade2 ctc2 = new ChemicalToCascade2();

        ep.initiate();
        ec.initiate();
        er.initiate();
        ecm.initiate();
        ctc2.initiate();

        String chempath = new C5().getClass().getResource("Data" + "/" + "compounds.dat").getFile();
        String rxnpath = new C5().getClass().getResource("Data" + "/" + "reactions.dat").getFile();
        String CoreMetpath = new C5().getClass().getResource("Data" + "/" + "e_coli_core_metabolites.csv").getFile();
        String PathPath = new C5().getClass().getResource("Data" + "/" + "pathways.dat").getFile();

        List<Chems> listofchems = ec.run(chempath);
        List<String> aminoacids = ec.getAAlists();
        ep.run(PathPath);
        rxnsHashMap = er.getRxnsHashmap();
        chems = ec.getChemsHashMap();
        reactions = er.run(rxnpath,chems);
        PathMap = ep.getPathMap();
        natives = ecm.run(CoreMetpath);
        natives.addAll(aminoacids);
        chemToCascadeMap = ctc2.run(reactions, listofchems);
    }

    public List<List<Rxns>> run(Cascade2 cascade){
        Chems product = cascade.getProduct();
        List<Rxns> CurrPath = new ArrayList<>();
        List<List<Rxns>> enumPath = new ArrayList<>();
        Set<String> visited = new HashSet<>();
        Set<String> visited_rxn = new HashSet<>();
        visited.add(product.getID());
        depthSearch(product,enumPath,CurrPath,visited,visited_rxn, 0);
        return enumPath;
    }

    private void depthSearch(Chems chem, List<List<Rxns>> allPaths,List<Rxns> CurrPath,Set<String> visitedChem, Set<String> visited_rxn, int layer){
        Cascade2 cascade = chemToCascadeMap.get(chem);
        if (layer > 6 || cascade == null){
            return;
        }

        if (cascade.getRxnsThatFormPdt().size() == 0){
            allPaths.add(CurrPath);
            return;
        }

        for (Rxns r : cascade.getRxnsThatFormPdt()) {
            List<Rxns> helperlist = new ArrayList<>(CurrPath);

            boolean visited = false;
            for (Chems substrate : r.getSubstrates()) {
                if (visitedChem.contains(substrate.getID())) {
                    visited = true;
                }
            }
            if (visited) {
                continue;
            }
            helperlist.add(r);

            if (allNatives(r.getSubstrates()) || isNativeRxn(r)) {
                List<Rxns> newPath = new ArrayList<>(helperlist);
                allPaths.add(newPath);
                continue;
            } else {
                for (Chems c : r.getSubstrates()){
                    if (isNatives(c)){
                        continue;
                    }
                    if (visitedChem.contains(c.getID())) {
                        helperlist.remove(r);
                        break;
                    }
                    if (c.getInchi() == null){
                        continue;
                    }
                    else {
                        visitedChem.add(c.getID());
                        depthSearch(c, allPaths, helperlist, visitedChem, visited_rxn, layer + 1);
                        visitedChem.remove(c.getID());
                        continue;
                    }

                }
            }
        }
    }

    private boolean isNatives(Chems chem) {
        Boolean ret = false;
        if (natives.contains(chem.getID())){
            ret = true;
        }
        return ret;
    }

    private boolean allNatives(Set<Chems> chems){
        Boolean ret = true;
        for (Chems c : chems){
            if (!natives.contains(c.getID())){
                ret = false;
            }
        }
        return ret;
    }

    private Boolean isNativeRxn(Rxns r){
        boolean ret = false;
        if (!r.inPathway()){
            return ret;
        } else {
            for (String str : r.getPathways()){
                Paths p = PathMap.get(str);
                if (p.containSpecies("Escherichia coli K-12 substr. MG1655")) {
                    ret = true;
                    break;
                }
            }
        }
        return ret;
    }

    public HashMap<Chems, Cascade2> getChemToCascadeMap() {
        return chemToCascadeMap;
    }

    public HashMap<String, Chems> getChems() {
        return chems;
    }

}
