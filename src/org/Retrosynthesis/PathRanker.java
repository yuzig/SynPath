package org.Retrosynthesis;

import org.Retrosynthesis.models.Rxns;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.List;
import java.util.stream.Collectors;

public class PathRanker  {

    public void initiate(){

    }

    public List<List<Rxns>> runCobraPy(String chasismodelPath, String listOfPaths) throws Exception {
        ProcessBuilder processBuilder = new ProcessBuilder("python3", "/Users/carol_gyz/IdeaProjects/SBOLmetPathDesign/cobrapyconverter/GenomeScaleModels/CobraConverter.py",  chasismodelPath, listOfPaths);
        Process process = processBuilder.start();
//        List<String> results = readProcessOutput(process.getInputStream());
        BufferedReader reader = new BufferedReader((new InputStreamReader(process.getInputStream())));
        BufferedReader errorReader = new BufferedReader(new InputStreamReader(process.getErrorStream()));
        String lines = null;
        while ((lines = reader.readLine()) != null) {
            System.out.println(lines);
        }
//        while ((lines = errorReader.readLine()) != null) {
//            System.out.println(lines);
//        }
        return null;
    }
//        List<String> results = readProcessOutput(process.getInputStream());
//        return  null;

    private List<String> readProcessOutput(InputStream inputStream) throws IOException {
        try (BufferedReader output = new BufferedReader(new InputStreamReader(inputStream))) {
            return output.lines()
                    .collect(Collectors.toList());
        }
    }
}
