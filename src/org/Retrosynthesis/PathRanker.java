package org.Retrosynthesis;


import java.io.*;
import java.util.*;
import java.util.stream.Collectors;
import org.Retrosynthesis.ArrayIndexComparator;
import org.Retrosynthesis.models.RankedPathsObj;
import org.Retrosynthesis.models.Rxns;

public class PathRanker  {

    public void initiate(){

    }

    public Set<RankedPathsObj> runCobraPy(String chasismodelPath, String listOfPaths) throws Exception {
        String file = "file.txt";
        Pathway2Files(listOfPaths, file);
        String pathwayFile = "/Users/carol_gyz/IdeaProjects/SBOLmetPathDesign/file.txt";

        ProcessBuilder processBuilder = new ProcessBuilder("python3", "/Users/carol_gyz/IdeaProjects/SBOLmetPathDesign/cobrapyconverter/GenomeScaleModels/CobraConverter.py",  chasismodelPath, pathwayFile);
        processBuilder.redirectErrorStream(true);
        String[] paths = listOfPaths.split("//");
        Process process = processBuilder.start();
        Set<RankedPathsObj> ret = new HashSet<>();
        List<String> results = readProcessOutput(process.getInputStream());

        for (int i = 3; i < results.size(); i ++) {
            String[] tabs = results.get(i).split("\t");
            RankedPathsObj rankedPathsObj = new RankedPathsObj(Double.parseDouble(tabs[0]), Double.parseDouble(tabs[1]), Double.parseDouble(tabs[2]), Double.parseDouble(tabs[3]), Double.parseDouble(tabs[4]), paths[Integer.parseInt(tabs[5])]);
            ret.add(rankedPathsObj);
        }
        return ret;
    }

    private List<String> readProcessOutput(InputStream inputStream) throws IOException {
        try (BufferedReader output = new BufferedReader(new InputStreamReader(inputStream))) {
            return output.lines()
                    .collect(Collectors.toList());
        }
    }

    public static void Pathway2Files(String Pathways, String fileName) throws IOException {
        FileWriter writer1 = new FileWriter(fileName);
        writer1.write(Pathways);
        writer1.close();
    }
}
