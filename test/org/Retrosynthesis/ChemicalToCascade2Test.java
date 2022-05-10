package org.Retrosynthesis;

import org.C5;
import org.Retrosynthesis.models.Cascade2;
import org.Retrosynthesis.models.Chems;
import org.Retrosynthesis.models.Rxns;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.HashMap;
import java.util.List;

import static org.junit.jupiter.api.Assertions.*;

class ChemicalToCascade2Test {
    private ExtractChem ec;
    private ExtractRxns er;
    private List<Rxns> output;
    private List<Chems> chems;

    @BeforeEach
    void setUp() throws Exception {
        String chempath = new C5().getClass().getResource("Data" + "/" + "compounds.dat").getFile();
        String rxnpath = new C5().getClass().getResource("Data" + "/" + "reactions.dat").getFile();
        ec = new ExtractChem();
        ec.initiate();
        chems = ec.run(chempath);
        HashMap<String, Chems> chems = ec.getChemsHashMap();

        er = new ExtractRxns();
        er.initiate();
        output = er.run(rxnpath,chems);
    }

    @Test
    void run() throws Exception{
        ChemicalToCascade2 ctc2 = new ChemicalToCascade2();
        ctc2.initiate();
        HashMap<Chems, Cascade2> cascades = ctc2.run(output, chems);
        Chems cannabidiol = ec.getChemsHashMap().get("CPD-7173");
        Cascade2 cannabidiol_cascade = cascades.get(cannabidiol);
        assertTrue(cannabidiol_cascade.getRxnsThatFormPdt().size() == 1);

        Chems hydroxybutrylcoa = ec.getChemsHashMap().get("4-HYDROXY-BUTYRYL-COA");
        Cascade2 hydroxybutylcoa_cascade = cascades.get(hydroxybutrylcoa);
        assertTrue(hydroxybutylcoa_cascade.getRxnsThatFormPdt().size() == 3);


    }
}