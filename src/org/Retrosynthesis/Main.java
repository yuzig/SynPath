package org.Retrosynthesis;

import org.C5;
import org.Retrosynthesis.models.*;
import java.io.FileWriter;

import java.io.IOException;
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

        Chems butanediol = chemMaps.get("BUTANOL");
        Set<List<Rxns>> output = pe.run(cascadeMap.get(butanediol));
        Pathway2Files(output);

//        Pathway2SBOLdoc pathway2SBOLdoc = new Pathway2SBOLdoc();
//        pathway2SBOLdoc.initiate();
//        List<List<Rxns>> output_list = output.stream().toList();
//
//        for (int i = 0; i < 10; i ++) {
//            pathway2SBOLdoc.run(output_list.get(i), i);
//        }
//
//        pathway2SBOLdoc.out("test.xml");



    }
    public static void Pathway2Files(Set<List<Rxns>> output) throws IOException {
        StringBuilder sb = new StringBuilder();
        for(List<Rxns> list : output) {
            for (Rxns r : list) {
                sb.append(r.getName() + " ");
                sb.append(r.toString());
                sb.append('\n');
            }
            sb.append('\n');
        }
        String out = sb.toString();
        FileWriter writer1 = new FileWriter("output_butanol.txt");
        writer1.write(out);
        writer1.close();
    }

}

