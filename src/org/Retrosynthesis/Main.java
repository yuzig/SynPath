package org.Retrosynthesis;

import org.C5;
import org.Retrosynthesis.models.*;


import java.util.ArrayList;
import java.util.List;
import java.util.Set;

public class Main {
    public static void main(String[] args) throws Exception {
        PathwayEnumerator pathwayEnumerator = new PathwayEnumerator();
        pathwayEnumerator.initiate();

        Pathway2SBOLdoc converter = new Pathway2SBOLdoc();
        converter.initiate();

        Chemical butanol = pathwayEnumerator.getChemMaps().get("1-butanol");
        Cascade butanolCascade = pathwayEnumerator.getChemicalToCascade().get(butanol);
        List<Pathway> output = pathwayEnumerator.run(butanolCascade);

//         for (int i = 0; i < 10; i ++) {
//            converter.run(output.get(i), i);
//        }

        converter.out("test.xml");






    }

}

