package org.Retrosynthesis;

import org.Retrosynthesis.models.Chems;
import org.Retrosynthesis.models.Paths;
import org.Retrosynthesis.models.Pathway;
import org.Retrosynthesis.models.Rxns;
import org.Utils.FileUtils;

import java.util.*;

public class ExtractPaths {
    private HashMap<String, Paths> PathMap;

    public void initiate(){
        PathMap = new HashMap<>();
    }

    public List<Paths> run(String filepath ) throws Exception {
        List<Paths> allPaths = new ArrayList<>();
        String pathdata = FileUtils.readFile(filepath);
        String[] lines = pathdata.trim().split("//");
        pathdata = null;
        String uniqueID = null;
        String name = null;
        Set<String> species = new HashSet<>();

        for (int i = 2; i < lines.length; i++) {
            String apath = lines[i];
            String[] aline = apath.trim().split("\\r|\\r?\\n");
            String[] tabs;

            for (String str : aline) {
                if (str.startsWith("UNIQUE-ID")) {
                    tabs = str.split(" - ");
                    uniqueID = tabs[1];
                    continue;
                }
                if (str.startsWith("COMMON-NAME")) {
                    tabs = str.split(" - ");
                    name = tabs[1];
                    continue;
                }
                if (str.startsWith("SPECIES")){
                    tabs = str.split(" - ");
                    species.add(tabs[1]);
                }
            }
            Paths path = new Paths(uniqueID,name,species);
            allPaths.add(path);
            PathMap.put(uniqueID,path);
            species = new HashSet<>();
        }
        return allPaths;
    }

    public HashMap<String, Paths> getPathMap() {
        return PathMap;
    }
}
