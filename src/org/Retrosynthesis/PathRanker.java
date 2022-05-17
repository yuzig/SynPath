package org.Retrosynthesis;


import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;
import org.Retrosynthesis.ArrayIndexComparator;

public class PathRanker  {

    public void initiate(){

    }

    public HashMap<Integer, Double>  runCobraPy(String chasismodelPath, String listOfPaths) throws Exception {
        ProcessBuilder processBuilder = new ProcessBuilder("python3", "/Users/carol_gyz/IdeaProjects/SBOLmetPathDesign/cobrapyconverter/GenomeScaleModels/CobraConverter.py",  chasismodelPath, listOfPaths);

        Process process = processBuilder.start();
        List<String> results = readProcessOutput(process.getInputStream());
        Double[] ranked_pwy = new Double[results.size()-3];
        int c = 0;
        for (int i = 3; i < results.size(); i ++) {
           ranked_pwy[c] = (Double.parseDouble(results.get(i)));
           c++;
        }

        ArrayIndexComparator comparator = new ArrayIndexComparator(ranked_pwy);
        Integer[] indexes = comparator.createIndexArray();
        Arrays.sort(indexes, comparator);
        HashMap<Integer, Double> ret = new HashMap<>();
        for (int i = indexes.length - 1; i > indexes.length - 10; i--) {
            ret.put(indexes[i], ranked_pwy[indexes[i]]);
        }
        return ret;
    }

    private List<String> readProcessOutput(InputStream inputStream) throws IOException {
        try (BufferedReader output = new BufferedReader(new InputStreamReader(inputStream))) {
            return output.lines()
                    .collect(Collectors.toList());
        }
    }
}
