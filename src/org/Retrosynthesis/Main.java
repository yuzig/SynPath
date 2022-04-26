package org.Retrosynthesis;

import org.C5;
import org.Retrosynthesis.models.*;
import java.io.FileWriter;


import java.nio.file.Path;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

public class Main {
    private static ExtractChem ec;
    private static ExtractRxns er;
    private static ExtractCoreMet ECM;
    private static ExtractPaths ep;

    public static void main(String[] args) throws Exception {
        String chempath = new C5().getClass().getResource("Data" + "/" + "compounds.dat").getFile();
        String rxnpath = new C5().getClass().getResource("Data" + "/" + "reactions.dat").getFile();
        String CoreMetpath = new C5().getClass().getResource("Data" + "/" + "e_coli_core_metabolites.csv").getFile();
        String PathPath = new C5().getClass().getResource("Data" + "/" + "pathways.dat").getFile();


        ec = new ExtractChem();
        ec.initiate();
        ec.run(chempath);
        HashMap<String, Chems> chems = ec.getChemsHashMap();

        er = new ExtractRxns();
        List<Rxns> output = er.run(rxnpath,chems);

        ECM = new ExtractCoreMet();
        List<String> natives = ECM.run(CoreMetpath);

        ep = new ExtractPaths();
        ep.initiate();
        List<Paths> ListofPathways = ep.run(PathPath);


        int a = 1;


//        PathwayEnumerator pathwayEnumerator = new PathwayEnumerator();
//        pathwayEnumerator.initiate();
//
//        Pathway2SBOLdoc converter = new Pathway2SBOLdoc();
//        converter.initiate();
//
//        Chemical butanol = pathwayEnumerator.getChemMaps().get("cannabidiol");
//        Cascade butanolCascade = pathwayEnumerator.getChemicalToCascade().get(butanol);
//        List<Pathway> output = pathwayEnumerator.run(butanolCascade);
//
//
////         for (int i = 0; i < 10; i ++) {
////            converter.run(output.get(i), i);
////        }
//
//        FileWriter writer = new FileWriter("output_cannabidiol.txt");
//        for(Pathway p: output) {
//            String str = p.toString();
//            writer.write(str + System.lineSeparator());
//        }
//        writer.close();

//        converter.out("test.xml");






    }

}

