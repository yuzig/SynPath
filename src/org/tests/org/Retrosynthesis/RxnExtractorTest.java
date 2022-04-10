package org.Retrosynthesis;

import org.C5;
import org.Retrosynthesis.models.Chemical;
import org.Retrosynthesis.models.Reaction;
import org.junit.Test;

import java.util.ArrayList;
import java.util.List;


class RxnExtractorTest {

    private List<Chemical> allChems;
    private List<Reaction> allRxns;

    @Test
    void reactionFormatTest() throws Exception{
        ChemExtractor ce = new ChemExtractor();
        ce.initiate();
        RxnExtractor re = new RxnExtractor();
        re.initiate();
        MetaboliteExtractor me = new MetaboliteExtractor();
        me.initiate();


        String chempath = new C5().getClass().getResource("data" + "/" + "good_chems.txt").getFile();
        String rxnpath = new C5().getClass().getResource("data" + "/" + "good_reactions.txt").getFile();
        List<String> metPaths = new ArrayList<>();
        metPaths.add(new C5().getClass().getResource("data" + "/" + "minimal_metabolites.txt").getFile());
        metPaths.add(new C5().getClass().getResource("data" + "/" + "universal_metabolites.txt").getFile());
        allChems = ce.run(chempath);
        List<Reaction> allRxns = re.run(rxnpath, allChems);
        System.out.println(allRxns.get(0).getecnum());
//        assertTrue(true);
    }
}

