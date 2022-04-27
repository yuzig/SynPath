package org.Retrosynthesis;

import org.C5;
import org.Retrosynthesis.models.Cascade2;
import org.Retrosynthesis.models.Chems;
import org.Retrosynthesis.models.Rxns;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.HashMap;
import java.util.List;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.*;

class ExtractRxnsTest {
    private static ExtractChem ec;
    private static List<Chems> chems;
    private static ExtractRxns er;
    private static List<Rxns> output;
    private static ChemicalToCascade2 ctc2;
    private static HashMap<Chems, Cascade2> chems2cascades ;

    @BeforeEach
    void setUp() throws Exception {
        String chempath = new C5().getClass().getResource("Data" + "/" + "compounds.dat").getFile();
        String rxnpath = new C5().getClass().getResource("Data" + "/" + "reactions.dat").getFile();
        ec = new ExtractChem();
        ec.initiate();
        chems = ec.run(chempath);
        HashMap<String, Chems> chemsHashMap = ec.getChemsHashMap();

        er = new ExtractRxns();
        output = er.run(rxnpath,chemsHashMap);

        ctc2 = new ChemicalToCascade2();
        ctc2.initiate();
        chems2cascades = ctc2.run(output,chems);

    }

    @Test
    void run() {
        HashMap<String, Chems> map = ec.getChemsHashMap();
        Chems butanediol = map.get("CPD-13560");
        Cascade2 cascade = chems2cascades.get(butanediol);
        assertTrue(cascade.getRxnsThatFormPdt().size() == 1);
    }
}