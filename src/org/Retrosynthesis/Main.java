package org.Retrosynthesis;

import org.Retrosynthesis.models.*;

import java.io.FileWriter;

import java.io.IOException;
import java.util.HashMap;
import java.util.List;
import java.util.Scanner;
import java.util.Set;


public class Main {

    private static PathEnum pe;
    private static PathRanker pr;

    static{
        Scanner readData = new Scanner(System.in);
        System.out.println("Enter target ID, precursor ID, and path to chassis separated by space :");
        String str = readData.nextLine();
        String [] str1 = str.split(" ");
        try {
            Main.main(str1);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static void main(String[] args) throws Exception {
        String chemID = args[0];
        String precursor = null;
        String chassismodelPath = null;
        if (args.length == 3) {
            precursor = args[1];
            chassismodelPath = args[2];
        } else {
            precursor = null;
            chassismodelPath = args[1];
        }
        pe = new PathEnum();
        pe.initiate();
        pr = new PathRanker();
        pr.initiate();
//        Pathway2SBOLdoc pathway2SBOLdoc = new Pathway2SBOLdoc();
//        pathway2SBOLdoc.initiate(pe.getRxnsHashMap());

        HashMap<Chems, Cascade2> cascadeMap = pe.getChemToCascadeMap();
        HashMap<String, Chems> chemMaps = pe.getChems();


        Chems chemical = chemMaps.get(chemID);
        List<List<Rxns>> output = pe.run(cascadeMap.get(chemical), precursor);
        String pathways = Pathway2String(output);Set<RankedPathsObj> out = pr.runCobraPy(chassismodelPath,pathways);
        for (RankedPathsObj obj : out) {
            System.out.println("Theoretical Yield: " + obj.getTheoretical_yield());
            System.out.println("atp consumption: " + obj.getAtpm_generation());
            System.out.println("nadph_consumption: " + obj.getNadph_consumption());
            System.out.println("nadh_consumption: " + obj.getNadh_consumption());
            System.out.println("eng_yield: " + obj.getFva_dif());
            System.out.println("pathway: " + obj.getPathway());
        }
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

