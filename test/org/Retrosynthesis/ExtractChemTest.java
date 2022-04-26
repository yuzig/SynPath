package org.Retrosynthesis;

import org.C5;
import org.Retrosynthesis.models.Chems;
import org.junit.BeforeClass;
import org.junit.jupiter.api.Test;

import java.util.HashMap;

import static org.junit.jupiter.api.Assertions.*;

class ExtractChemTest {
    private static ExtractChem ec;
    private static HashMap<String, Chems> chems;

    @BeforeClass
    public static void setUpClass() throws Exception {
//        ec = new ExtractChem();
//        ec.initiate();
//        String chempath = new C5().getClass().getResource("Data" + "/" + "compounds.dat").getFile();
//        ec.run(chempath);
//        chems = ec.getChemsHashMap();
    }
    @Test
    void run() throws Exception {
        ec = new ExtractChem();
        ec.initiate();
        String chempath = new C5().getClass().getResource("Data" + "/" + "compounds.dat").getFile();
        ec.run(chempath);
        chems = ec.getChemsHashMap();

        Chems butanediol = new Chems("CPD-13560","InChI=1S/C4H10O2/c5-3-1-2-4-6/h5-6H,1-4H2","1,4-butanediol");
        Chems cannabidiol = new Chems("CPD-7173", "InChI=1S/C21H30O2/c1-5-6-7-8-16-12-19(22)21(20(23)13-16)18-11-15(4)9-10-17(18)14(2)3/h11-13,17-18,22-23H,2,5-10H2,1,3-4H3/t17-,18+/m0/s1","cannabidiol");
        Chems pyruvate = new Chems("PYRUVATE", "InChI=1S/C3H4O3/c1-2(4)3(5)6/h1H3,(H,5,6)/p-1","pyruvate");

        Chems but = chems.get(butanediol.getID());
        assertEquals(but.getInchi(),butanediol.getInchi());
        assertEquals(but.getName(), butanediol.getName());

        Chems cann = chems.get(cannabidiol.getID());
        assertEquals(cann.getInchi(),cannabidiol.getInchi());
        assertEquals(cann.getName(), cannabidiol.getName());

        Chems pyr = chems.get(pyruvate.getID());
        assertEquals(pyr.getInchi(),pyruvate.getInchi());
        assertEquals(pyr.getName(), pyruvate.getName());

    }

    @Test
    void getChemsHashMap() {
    }
}