package org.Retrosynthesis;

import org.Retrosynthesis.models.Cascade2;
import org.Retrosynthesis.models.Chems;
import org.Retrosynthesis.models.Rxns;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import static org.junit.jupiter.api.Assertions.*;

class PathEnumTest {
    private PathEnum PE;
    private HashMap<Chems, Cascade2> cascadeMap;
    private HashMap<String, Chems> chemMaps;

    @BeforeEach
    void setUp() throws Exception {
        PE = new PathEnum();
        PE.initiate();
        cascadeMap = PE.getChemToCascadeMap();
        chemMaps = PE.getChems();
    }

    @Test
    void run() {
        Chems c = chemMaps.get("CPD-13560");
        List<List<Rxns>> output = PE.run(cascadeMap.get(c), null);
        for (List<Rxns> pathways : output) {
            StringBuilder sb = new StringBuilder();
            for(Rxns r: pathways) {
                for (Chems sub : r.getSubstrates()){
                    sb.append(sub.getName());
                    sb.append(" ");
                }
                sb.append(" -->");
                for (Chems pro : r.getProducts()){
                    sb.append(pro.getName());
                    sb.append(" ");
                }
                sb.append("\n");
            }
            String print = sb.toString();
            System.out.println(print);
        }
        assertTrue(true);
    }
}