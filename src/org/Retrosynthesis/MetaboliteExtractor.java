package org.Retrosynthesis;


import org.Utils.FileUtils;

import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class MetaboliteExtractor {

    public void initiate() {
    }

    public Set<String> run(List<String> metabolitePaths) throws Exception {

        Set<String> nativeInchis = new HashSet<String>();

        for (String path : metabolitePaths) {
            String nativedata = FileUtils.readFile(path);
            nativedata = nativedata.replaceAll("\"", "");
            String[] lines = nativedata.trim().split("\\r|\\r?\\n");
            nativedata = null;

            //Each line of the file is a chemical, add it to the list
            for (int i = 1; i < lines.length; i++) {
                String aline = lines[i];
                String[] tabs = aline.split("\t");
                String inchi = tabs[1];
                nativeInchis.add(inchi);
            }
        }

        return nativeInchis;
    }

}
