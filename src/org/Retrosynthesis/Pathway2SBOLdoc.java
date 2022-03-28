package org.Retrosynthesis;
import org.C5;
import org.Retrosynthesis.models.*;
import org.sbolstandard.core2.SBOLConversionException;
import org.sbolstandard.core2.SBOLDocument;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

public class Pathway2SBOLdoc {

    private List<Chemical> allChems;
    private List<Reaction> allRxns;

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
        allRxns = re.run(rxnpath, allChems);
    }

    private SBOLDocument run(Pathway pathway) throws SBOLConversionException{

        return null;
    }
}



//    private final SBOLDocument sboldoc;
//    private final List<Reaction> pathwayIn;
//    private final List<Chemical> Chemicals;
//
//    public Pathway2SBOLdoc (Pathway pathway){
//        pathwayIn = pathway.getReactions();
//
//        }
//
//    public static void main(String[] args) {
//
//    }
//    }
//
//

