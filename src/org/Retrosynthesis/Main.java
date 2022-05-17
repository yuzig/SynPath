package org.Retrosynthesis;

import org.C5;
import org.Retrosynthesis.models.*;
import java.io.FileWriter;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

public class Main {

    private static PathEnum pe;
    private static PathRanker pr;

    public static void main(String[] args) throws Exception {

        pe = new PathEnum();
        pe.initiate();
        pr = new PathRanker();
        pr.initiate();
        HashMap<Chems, Cascade2> cascadeMap = pe.getChemToCascadeMap();
        HashMap<String, Chems> chemMaps = pe.getChems();

        Chems butanediol = chemMaps.get("CPD-13560");
        List<List<Rxns>> output = pe.run(cascadeMap.get(butanediol));
        String pathways = Pathway2String(output);
        HashMap<Integer, Double> out = pr.runCobraPy("/Users/carol_gyz/IdeaProjects/SBOLmetPathDesign/cobrapyconverter/GenomeScaleModels/iJO1366.xml",pathways);
        for (Integer key : out.keySet()) {
            List<Rxns> path = output.get(key);
            String pathway = OnePathway2String(path);
            Double theoretical_yield = out.get(key);

            System.out.println("Theoretical Yield:" + theoretical_yield);
            System.out.println(pathway);
        }


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
    public static String OnePathway2String(List<Rxns> output) {
        StringBuilder sb = new StringBuilder();
        for(Rxns r : output) {
            sb.append(r.getName() + "\t");
            sb.append(r.toString());
            sb.append('\n');
        }
        return sb.toString();
    }
    public static String Pathway2String(List<List<Rxns>> output) {
        StringBuilder sb = new StringBuilder();
        for(List<Rxns> list : output) {
            for (Rxns r : list) {
                sb.append(r.getName() + "\t");
                sb.append(r.toString());
                sb.append('\n');
            }
            sb.append("//");
        }
        return sb.toString();
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

