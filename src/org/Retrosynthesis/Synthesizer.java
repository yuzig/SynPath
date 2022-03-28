package org.Retrosynthesis;
import org.Retrosynthesis.models.Chemical;
import org.Retrosynthesis.models.HyperGraph;
import org.Retrosynthesis.models.Reaction;

import java.util.HashMap;
import java.util.List;
import java.util.Set;

public class Synthesizer {

    public void initiate() throws Exception {}

    public HyperGraph run(List<Reaction> allRxn, List<Chemical> allChem, Set<String> natives) throws Exception {
        HashMap<Chemical, Integer> chemicalToShell = new HashMap<Chemical, Integer>();
        HashMap<Reaction, Integer> reactionToShell = new HashMap<Reaction, Integer>();

        int currShell = 1;

        for (Chemical chem : allChem) {
            if (natives.contains(chem.getInchi())) {
                chemicalToShell.put(chem, 0);
            }
        }

        while (ExpandOnce(currShell, allRxn, chemicalToShell, reactionToShell)) {
            currShell++;
        }

        return new HyperGraph(reactionToShell, chemicalToShell);
    }

    private boolean ExpandOnce(int currshell, List<Reaction> allReactions, HashMap<Chemical, Integer> chemicalToShell,
                               HashMap<Reaction, Integer> reactionToShell) throws Exception {
        //Increment the current shell
        boolean isExpanded = false;

        //Iterate through reactions
        outer:
        for (Reaction rxn : allReactions) {
            //If the reaction has already been put in the expansion, skip this reaction
            if (reactionToShell.containsKey(rxn)) {
                continue outer;
            }

            //If any of the substates are not enabled, skip this reaction
            for (Chemical achem : rxn.getSubstrates()) {
                if (!chemicalToShell.containsKey(achem)) {
                    continue outer;
                }
            }

            //If gets this far, the Reaction is enabled and new, thus expansion will occur
            isExpanded = true;

            //Log the reaction into the expansion at the current shell
            reactionToShell.put(rxn, currshell);

            //For each product, enable it with current shell (if it isn't already)
            for (Chemical chemid : rxn.getProducts()) {
                if (!chemicalToShell.containsKey(chemid)) {
                    chemicalToShell.put(chemid, currshell);
                }
            }
        }

        System.out.println("Expanded shell: " + currshell + " result " + isExpanded + " with " + chemicalToShell.size() + " reachables");

        return isExpanded;
    }

}
