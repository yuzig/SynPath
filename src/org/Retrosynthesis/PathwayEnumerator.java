package org.Retrosynthesis;

import org.C5;
import org.Retrosynthesis.models.*;

import java.sql.Array;
import java.util.*;

/**
 * Enumarate all pathways that produces a product from native metabolites
 * @author carol_gyz Carol Gao
 */

public class PathwayEnumerator {
    private HashMap<Chemical, Integer> chemicalToShell;
    private HashMap<Chemical, Cascade> chemicalToCascade;
    private List<Chemical> allChems ;
    private HashMap<Chemical, Cascade> chemtoCascade;
    private HashMap<String, Chemical> chemMap;
    private Set<String> natives;

    public void initiate() throws Exception {
        ChemExtractor ce = new ChemExtractor();
        ce.initiate();
        RxnExtractor re = new RxnExtractor();
        re.initiate();
        MetaboliteExtractor me = new MetaboliteExtractor();
        me.initiate();
        Synthesizer s = new Synthesizer();
        s.initiate();
        ChemicalToCascade CC = new ChemicalToCascade();
        CC.initiate();

        String chempath = new C5().getClass().getResource("data" + "/" + "good_chems.txt").getFile();
        String rxnpath = new C5().getClass().getResource("data" + "/" + "good_reactions.txt").getFile();
        List<String> metPaths = new ArrayList<>();
        metPaths.add(new C5().getClass().getResource("data" + "/" + "minimal_metabolites.txt").getFile());
        metPaths.add(new C5().getClass().getResource("data" + "/" + "universal_metabolites.txt").getFile());
        allChems = ce.run(chempath);
        chemMap = ce.getChemicalHashMap();
        List<Reaction> allRxns = re.run(rxnpath, allChems);
        natives = me.run(metPaths);
        HyperGraph hg = s.run(allRxns, allChems, natives);
        chemtoCascade = CC.run(allRxns,allChems);

        chemicalToShell = hg.getChemicalToShell();
        chemicalToCascade = hg.getChemicalToCascade();
    }

    public List<Pathway> run(Cascade cascade) throws Exception {
        Chemical product = cascade.getProduct();
        List<Pathway> enumPath = new ArrayList<>();
        List<Reaction> CurrPath = new ArrayList<>();
        Set<String> visited = new HashSet<String>();
        visited.add(product.getInchi());
        depthSearch(product,enumPath, CurrPath, visited, 0);
        System.out.println(enumPath);
        return enumPath;
    }

    private void depthSearch(Chemical chem, List<Pathway> allPaths, List<Reaction> CurrPath, Set<String> visitedChem, int layer) {
        Cascade cascade = chemtoCascade.get(chem);
        if (layer > 5) {
            return;
        }
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
                    if (visitedChem.contains(c.getInchi())){
                        continue;
                    } else {
                        visitedChem.add(c.getInchi());
                        depthSearch(c, allPaths, helperlist, visitedChem, layer + 1);
                        visitedChem.remove(c.getInchi());
                    }
                }
            }
        }
    }

    private boolean allNatives(Set<Chemical> chems){
        Boolean ret = true;
        for (Chemical c : chems){
            if (!natives.contains(c.getInchi())){
                ret = false;
            }
        }
        return ret;
    }

    private boolean isNatives(Chemical chem) {
        Boolean ret = false;
        if (natives.contains(chem.getInchi())){
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

    public Map<String, Chemical> getChemMaps() {
        return chemMap;
    }
}


