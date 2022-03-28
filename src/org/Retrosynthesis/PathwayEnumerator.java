package org.Retrosynthesis;

import org.C5;
import org.Retrosynthesis.models.*;

import java.util.*;

/**
 * Enumarate all pathways that produces a product from native metabolites
 * @author carol_gyz Carol Gao
 */

public class PathwayEnumerator {
    private HashMap<Chemical, Integer> chemicalToShell;
    private HashMap<Chemical, Cascade> chemicalToCascade;
    private List<Chemical> allChems ;

    public void initiate() throws Exception {
        ChemExtractor ce = new ChemExtractor();
        ce.initiate();
        RxnExtractor re = new RxnExtractor();
        re.initiate();
        MetaboliteExtractor me = new MetaboliteExtractor();
        me.initiate();
        Synthesizer s = new Synthesizer();
        s.initiate();

        String chempath = new C5().getClass().getResource("data" + "/" + "good_chems.txt").getFile();
        String rxnpath = new C5().getClass().getResource("data" + "/" + "good_reactions.txt").getFile();
        List<String> metPaths = new ArrayList<>();
        metPaths.add(new C5().getClass().getResource("data" + "/" + "minimal_metabolites.txt").getFile());
        metPaths.add(new C5().getClass().getResource("data" + "/" + "universal_metabolites.txt").getFile());
        allChems = ce.run(chempath);
        List<Reaction> allRxns = re.run(rxnpath, allChems);
        Set<String> natives = me.run(metPaths);
        HyperGraph hg = s.run(allRxns, allChems, natives);

        chemicalToShell = hg.getChemicalToShell();
        chemicalToCascade = hg.getChemicalToCascade();
    }

    public List<Pathway> run(Cascade cascade) throws Exception {
        Chemical product = cascade.getProduct();
        List<Pathway> enumPath = new ArrayList<>();
        List<Reaction> CurrPath = new ArrayList<>();
        Set<Chemical> visited = new HashSet<>();
        depthSearch(product,enumPath, CurrPath, visited);
        return enumPath;
    }

    private void depthSearch(Chemical chem, List<Pathway> allPaths, List<Reaction> CurrPath, Set<Chemical> visited) {
        visited.add(chem);
        Cascade cascade = chemicalToCascade.get(chem);
            for (Reaction r : cascade.getRxnsThatFormPdt()) {

                List<Reaction> helperlist = new ArrayList<>(CurrPath);
                helperlist.add(r);
                if (allNatives(r.getSubstrates())){
                    Pathway path = new Pathway(helperlist);
                    allPaths.add(path);
                    return;

                } else {
                    for (Chemical c : r.getSubstrates()) {

                        if(isNatives(c)){
                            continue;
                        }
                        if (visited.contains(c)){
                            continue;
                        }
                        depthSearch(c, allPaths, helperlist, visited);

                    }
                }
            }
        }

    private boolean allNatives(Set<Chemical> chems){
        Boolean ret = true;
        for (Chemical c : chems){
            if (chemicalToShell.get(c) != 0){
                ret = false;
            }
        }
        return ret;
    }

    private boolean isNatives(Chemical chem) {
        Boolean ret = false;
        if (chemicalToShell.get(chem) == 0){
            ret = true;
        }
        return ret;
    }


    public HashMap<Chemical, Cascade>  getChemicalToCascade() {
        return chemicalToCascade;
    }
    public List<Chemical> getAllChems() {
        return allChems;
    }


//    public static void main(String[] args) throws Exception {
//        PathwayEnumerator PE = new PathwayEnumerator();
//        PE.initiate();
//
//        Chemical butanol =  PE.allChems.get(5133);
//        Cascade butanolCascade = PE.chemicalToCascade.get(butanol);
//
//        List<Pathway> output = PE.run(butanolCascade);
//        for (Pathway path : output){
//            System.out.println(">>>");
//            path.printPathway();
//            System.out.println(">>>");
//        }
//    }
}


