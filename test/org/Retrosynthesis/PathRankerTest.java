package org.Retrosynthesis;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

class PathRankerTest {
    private PathRanker pathranker;

    @BeforeEach
    void setUp() {
        pathranker = new PathRanker();
        pathranker.initiate();
    }

    @Test
    void runCobraPy() throws Exception {
        String modelpath = "/Users/carol_gyz/IdeaProjects/SBOLmetPathDesign/cobrapyconverter/GenomeScaleModels/iJO1366.xml";
        String apath = "RXN-12595\tPROTON CPD-13555 NADH --> CPD-13560 NAD\n"
                + "RXN-12594\tPROTON 4-HYDROXY-BUTYRYL-COA NADH --> CPD-13555 NAD CO-A\n"
                + "RXN-9092\t4-HYDROXY-BUTYRATE ATP CO-A --> PPI 4-HYDROXY-BUTYRYL-COA AMP\n"
                + "4-HYDROXYBUTYRATE-DEHYDROGENASE-RXN\tSUCC-S-ALD PROTON NADH --> 4-HYDROXY-BUTYRATE NAD\n"
                + "RXN-13328\tGLYOX 4-AMINO-BUTYRATE --> SUCC-S-ALD GLY\n"
                + "GUANIDINOBUTYRASE-RXN\tWATER CPD-592 --> UREA 4-AMINO-BUTYRATE\n"
                + "GUANIDINOBUTANAMIDE-NH3-RXN\t4-GUANIDO-BUTYRAMIDE WATER --> CPD-592 AMMONIUM\n"
                + "ARGININE-2-MONOOXYGENASE-RXN\tOXYGEN-MOLECULE ARG --> 4-GUANIDO-BUTYRAMIDE WATER CARBON-DIOXIDE\n";

        pathranker.runCobraPy(modelpath, apath);
        assertTrue(true);

    }
}