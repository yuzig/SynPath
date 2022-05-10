package org.Retrosynthesis;
import org.C5;
import org.Retrosynthesis.models.*;
import org.sbolstandard.core.datatree.Document;
import org.sbolstandard.core2.*;

import javax.xml.namespace.QName;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URI;
import java.util.*;
import java.lang.*;
import java.util.stream.Collectors;

public class Pathway2SBOLdoc {

    private HashMap<String, Chems> allChems;
    private Set<String> natives;
    private URI BIOPAX_BiochemicalRxn_uri;
    private URI InteractionBiochemicalRxnuri;
    private URI BIOPAX_SmallMolecule_uri;
    private URI ParticipationReactantUri;
    private URI ParticipationProductUri;
    private URI RoleProcessuri;
    private URI RoleChemicalUri;
    private URI ParticipationSideSubstrateUri;
    private URI ParticipationSideProductUri;
    private SBOLDocument document;
    private URI SequenceSmilesURI;
    private Set<Chems> visitedChem;
    private HashMap<String, ComponentDefinition> visitedComponents;
    private PathEnum pe;

    public void initiate() throws Exception {
        String uriPrefix = "http://sbols.org/dummy_uri";
        document = new SBOLDocument ();
        document.setDefaultURIprefix (uriPrefix);
        document.setComplete(true);
        document.setCreateDefaults (true);
        visitedChem = new HashSet<>();
        visitedComponents = new HashMap<>();

        BIOPAX_BiochemicalRxn_uri = URI.create("http://www.biopax.org/release/biopax-level3.owl#BiochemicalReaction");
        RoleProcessuri = URI.create("http://identifiers.org/so/SBO:0000375");
        RoleChemicalUri = URI.create("http://identifiers.org/so/SBO:0000247");
        InteractionBiochemicalRxnuri = URI.create("http://identifiers.org/biomodels.sbo/SBO:0000176");
        BIOPAX_SmallMolecule_uri = URI.create("http://www.biopax.org/release/biopax-level3.owl#SmallMolecule");
        ParticipationReactantUri = URI.create("http://identifiers.org/biomodels.sbo/SBO:0000010");
        ParticipationProductUri = URI.create("http://identifiers.org/biomodels.sbo/SBO:0000011");
        ParticipationSideSubstrateUri = URI.create("http://identifiers.org/biomodels.sbo/SBO:0000604");
        ParticipationSideProductUri = URI.create("http://identifiers.org/biomodels.sbo/SBO:0000603");
        SequenceSmilesURI = URI.create("http://www.opensmiles.org/opensmiles.html");

        pe = new PathEnum();
        pe.initiate();

        allChems = pe.getChems();
    }
    void run(List<Rxns> pathway, int i) throws SBOLValidationException, SBOLConversionException, IOException{
        List<Chems> metabolites = new ArrayList<>();
        ModuleDefinition moduleDef = document.createModuleDefinition("generic_pathway" + String.valueOf(i));
        HashMap<String, ComponentDefinition> componentDefinitionHashMap = new HashMap<>();
        ComponentDefinition Curr;

        for (Rxns r : pathway) {
            String str = r.getName().replaceAll("-","_");
            str = str.replaceAll("\\.","_");
            try {
                componentDefinitionHashMap.put("rxn_" + str,document.createComponentDefinition("rxn_" + str,BIOPAX_BiochemicalRxn_uri));
            } catch (SBOLValidationException ex) {
                break;
            }
            Curr = componentDefinitionHashMap.get("rxn_" + str);
            Curr.addRole(RoleProcessuri);
            Curr.createAnnotation(new QName("http://sbols.org/dummy_uri/annotation/ecnum#","ecnum","id"), r.getecnum());
            Curr.createAnnotation(new QName("http://sbols.org/dummy_uri/annotation/ecnum#","ecnum","biocyc"), r.getName());
            Curr.createAnnotation(new  QName("http://sbols.org/dummy_uri/annotation/ecnum#","ecnum","formula"), r.toString());
            for (Chems c: r.getSubstrates()) {
                String chem = c.getID().replaceAll("-","_");
                if (!componentDefinitionHashMap.containsKey("c_" + chem)){
                    componentDefinitionHashMap.put("c_" + chem, document.createComponentDefinition("c_" + chem,BIOPAX_SmallMolecule_uri));
                    Curr.addRole(RoleChemicalUri);
                    Curr.createAnnotation(new QName("http://sbols.org/dummy_uri/annotation/compound#", "Compound","biocyc"), c.getID());
                    Curr.addSequence(document.createSequence("seq_" + chem, "1", c.getSmiles(), SequenceSmilesURI));
                }
            }
            for (Chems c: r.getProducts()) {
                String chem = c.getID().replaceAll("-","_");
                if (!componentDefinitionHashMap.containsKey("c_" + chem)) {
                    componentDefinitionHashMap.put("c_" + chem, document.createComponentDefinition("c_" + chem, BIOPAX_SmallMolecule_uri));
                    Curr.addRole(RoleChemicalUri);
                    Curr.createAnnotation(new QName("http://sbols.org/dummy_uri/annotation/compound#", "Compound", "biocyc"), c.getID());
                    Curr.addSequence(document.createSequence( "_seq" + chem, "1", c.getSmiles(), SequenceSmilesURI));
                }
            }
        }
    }

    public void out(String filename) throws SBOLConversionException, IOException {
        document.write(filename);
        System.out.println("Done converting SBOLdoc to" + filename);
    }

}
