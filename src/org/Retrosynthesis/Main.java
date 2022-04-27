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

    private static PathEnum pe;

    public static void main(String[] args) throws Exception {

        pe = new PathEnum();
        pe.initiate();
        HashMap<Chems, Cascade2> cascadeMap = pe.getChemToCascadeMap();
        HashMap<String, Chems> chemMaps = pe.getChems();

        Chems butanediol = chemMaps.get("CPD-13560");
        Set<List<Rxns>> output = pe.run(cascadeMap.get(butanediol));
        StringBuilder sb = new StringBuilder();
        for(List<Rxns> list : output) {
            for (Rxns r : list) {
                sb.append(r.toString());
                sb.append('\n');
            }
            sb.append('\n');
        }
        String out = sb.toString();
        FileWriter writer1 = new FileWriter("output_butanediol.txt");
        writer1.write(out);
        writer1.close();




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

