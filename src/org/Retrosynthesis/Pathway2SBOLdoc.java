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
    private URI roleIntegrationURI;
    private URI SequenceSmilesURI;
    private Set<Chems> visitedChem;
    private HashMap<String, ComponentDefinition> visitedComponents;
    private PathEnum pe;
    private HashMap<String, Rxns> rxnsHashMap;

    public void initiate(HashMap<String, Rxns> rxnsHashmap) throws Exception {
        String uriPrefix = "http://sbols.org/dummy_uri";
        document = new SBOLDocument ();
        document.setDefaultURIprefix (uriPrefix);
        document.setComplete(true);
        document.setCreateDefaults (true);
        visitedChem = new HashSet<>();
        visitedComponents = new HashMap<>();
        this.rxnsHashMap = rxnsHashmap;

        BIOPAX_BiochemicalRxn_uri = URI.create("http://www.biopax.org/release/biopax-level3.owl#BiochemicalReaction");
        RoleProcessuri = URI.create("http://identifiers.org/so/SBO:0000375");
        RoleChemicalUri = URI.create("http://identifiers.org/so/SBO:0000247");
        roleIntegrationURI = URI.create("https://sbols.org/v2#overrideRoles");
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
    void run(RankedPathsObj pathway, int i) throws SBOLValidationException, SBOLConversionException, IOException{
        List<Chems> metabolites = new ArrayList<>();
        ModuleDefinition moduleDef = document.createModuleDefinition("generic_pathway" + String.valueOf(i));
        moduleDef.createAnnotation(new QName("http://sbols.org/dummy_uri/annotation/metrics#", "Theoretical_Yield", "cobrapy"), pathway.getTheoretical_yield());
        moduleDef.createAnnotation(new QName("http://sbols.org/dummy_uri/annotation/metrics#", "Eng_yield", "cobrapy"), 1);
        moduleDef.createAnnotation(new QName("http://sbols.org/dummy_uri/annotation/metrics#", "nadph_consumption", "cofactor"), 1);
        moduleDef.createAnnotation(new QName("http://sbols.org/dummy_uri/annotation/metrics#", "nadh_consumption", "cofactor"), 1);
        ComponentDefinition Curr;

        String path = pathway.getPathway();
        String[] tabs = path.split("\n");
        int a = 1;
        for (String r : tabs) {
            String rxnID = r.split("\t")[0];
            Rxns reaction = rxnsHashMap.get(rxnID);
            String str = reaction.getName().replaceAll("-","_");
            str = str.replaceAll("\\.","_");
            try {
                visitedComponents.put("rxn_" + str,document.createComponentDefinition("rxn_" + str,BIOPAX_BiochemicalRxn_uri));
            } catch (SBOLValidationException ex) {
                Curr = visitedComponents.get("rxn_" + str);
                moduleDef.createFunctionalComponent("Pathway_Step_" + String.valueOf(a), AccessType.PUBLIC, Curr.getPersistentIdentity(),DirectionType.NONE);
                continue;
            }
            Curr = visitedComponents.get("rxn_" + str);
            Curr.addRole(RoleProcessuri);
            QName order = new  QName("http://sbols.org/dummy_uri/annotation/reaction#", "reaction", "order");
            Curr.createAnnotation(order, a);
            Curr.createAnnotation(new QName("http://sbols.org/dummy_uri/annotation/reaction#","reaction","ecnum"), reaction.getecnum());
            Curr.createAnnotation(new QName("http://sbols.org/dummy_uri/annotation/reaction#","reaction","biocyc"), reaction.getName());
            Curr.createAnnotation(new  QName("http://sbols.org/dummy_uri/annotation/reaction#","reaction","formula"), r.toString());
            moduleDef.createFunctionalComponent("Pathway_Step_" + String.valueOf(a), AccessType.PUBLIC, Curr.getPersistentIdentity(),DirectionType.NONE);

            for (Chems c: reaction.getSubstrates()) {
                String chem = c.getID().replaceAll("-","_");
                ComponentDefinition cd = null;
                try {
                    visitedComponents.put("c_" + chem, document.createComponentDefinition("c_" + chem,BIOPAX_SmallMolecule_uri));
                } catch(SBOLValidationException e) {
                    cd = visitedComponents.get("c_" + chem);
                    Component component = Curr.createComponent("c_" + chem, AccessType.PUBLIC, cd.getPersistentIdentity());
                    component.setRoleIntegration(RoleIntegrationType.OVERRIDEROLES);
                    component.addRole(ParticipationSideSubstrateUri);
                    continue;
                }
                cd = visitedComponents.get("c_" + chem);
                cd.addRole(RoleChemicalUri);
                cd.createAnnotation(new QName("http://sbols.org/dummy_uri/annotation/compound#", "Compound","biocyc"), c.getID());
                if (c.getSmiles() != null) {
                    cd.addSequence(document.createSequence("seq_" + chem, "1", c.getSmiles(), SequenceSmilesURI));
                }
                Component component = Curr.createComponent("c_" + chem, AccessType.PUBLIC,cd.getPersistentIdentity());
                component.setRoleIntegration(RoleIntegrationType.OVERRIDEROLES);
                component.addRole(ParticipationReactantUri);

            }
            for (Chems c: reaction.getProducts()) {
                String chem = c.getID().replaceAll("-","_");
                if (!visitedComponents.containsKey("c_" + chem)) {
                    try {
                        visitedComponents.put("c_" + chem, document.createComponentDefinition("c_" + chem,BIOPAX_SmallMolecule_uri));
                    } catch(SBOLValidationException e) {
                        continue;
                    }
                    ComponentDefinition cd = visitedComponents.get("c_" + chem);
                    cd.addRole(RoleChemicalUri);
                    cd.createAnnotation(new QName("http://sbols.org/dummy_uri/annotation/compound#", "Compound", "biocyc"), c.getID());
                    if (c.getSmiles() != null) {
                        cd.addSequence(document.createSequence( "seq_" + chem, "1", c.getSmiles(), SequenceSmilesURI));
                    }
                    Component component = Curr.createComponent("c_" + chem, AccessType.PUBLIC,cd.getPersistentIdentity());
                    component.setRoleIntegration(RoleIntegrationType.OVERRIDEROLES);
                    component.addRole(ParticipationProductUri);
                }
            }

            a = a + 1;
        }
    }

    public void out(String filename) throws SBOLConversionException, IOException {
        document.write(filename);
        System.out.println("Done converting SBOLdoc to " + filename);
    }

}
