package org.Retrosynthesis;
import org.C5;
import org.Retrosynthesis.models.*;
import org.sbolstandard.core2.*;

import java.io.IOException;
import java.net.URI;
import java.util.*;

public class Pathway2SBOLdoc {

    private List<Chemical> allChems;
    private Set<String> natives;
    private URI BIOPAX_BiochemicalRxn_uri;
    private URI InteractionBiochemicalRxnuri;
    private URI BIOPAX_SmallMolecule_uri;
    private URI ParticipationReactantUri;
    private URI ParticipationProductUri;
    private URI RoleBiochemicalRxnuri;
    private URI ParticipationSideSubstrateUri;
    private URI ParticipationSideProductUri;
    private SBOLDocument document;
    private URI SequenceSmilesURI;
    private Set<Chemical> visitedChem;
    private HashMap<String, ComponentDefinition> visitedComponents;

    public void initiate() throws Exception {
        String uriPrefix = "http://sbols.org/dummy_uri";
        document = new SBOLDocument ();
        document.setDefaultURIprefix (uriPrefix);
        document.setComplete(true);
        document.setCreateDefaults (true);
        visitedChem = new HashSet<>();
        visitedComponents = new HashMap<>();

        BIOPAX_BiochemicalRxn_uri = URI.create("http://www.biopax.org/release/biopax-level3.owl#BiochemicalReaction");
        RoleBiochemicalRxnuri = URI.create("http://identifiers.org/so/SO:0000176");
        InteractionBiochemicalRxnuri = URI.create("http://identifiers.org/biomodels.sbo/SBO:0000176");
        BIOPAX_SmallMolecule_uri = URI.create("http://www.biopax.org/release/biopax-level3.owl#SmallMolecule");
        ParticipationReactantUri = URI.create("http://identifiers.org/biomodels.sbo/SBO:0000010");
        ParticipationProductUri = URI.create("http://identifiers.org/biomodels.sbo/SBO:0000011");
        ParticipationSideSubstrateUri = URI.create("http://identifiers.org/biomodels.sbo/SBO:0000604");
        ParticipationSideProductUri = URI.create("http://identifiers.org/biomodels.sbo/SBO:0000603");
        SequenceSmilesURI = URI.create("http://www.opensmiles.org/opensmiles.html");

        ChemExtractor ce = new ChemExtractor();
        ce.initiate();
        RxnExtractor re = new RxnExtractor();
        re.initiate();
        MetaboliteExtractor me = new MetaboliteExtractor();
        me.initiate();
        Synthesizer s = new Synthesizer();
        s.initiate();

        String chempath = new C5().getClass().getResource("data" + "/" + "good_chems.txt").getFile();
        String rxnpath = new C5().getClass().getResource("data" + "/" + "good_reactions.txt").getFile();
        List<String> metPaths = new ArrayList<>();
        metPaths.add(new C5().getClass().getResource("data" + "/" + "minimal_metabolites.txt").getFile());
        metPaths.add(new C5().getClass().getResource("data" + "/" + "universal_metabolites.txt").getFile());
        allChems = ce.run(chempath);
        natives = me.run(metPaths);
    }

    private void run(Pathway pathway) throws SBOLValidationException, SBOLConversionException, IOException {
        List<Reaction> pathwayIn = pathway.getReactions();
        List<Chemical> metabolites = new ArrayList<>();
        ModuleDefinition moduleDef = document.createModuleDefinition("pathway" + pathway.toString());

        int step = 0;

        for (Reaction r: pathwayIn){
            ComponentDefinition genericEnzyme= document.createComponentDefinition("genericStep" + String.valueOf(step) + "enzyme", BIOPAX_BiochemicalRxn_uri);
            genericEnzyme.addRole(RoleBiochemicalRxnuri);
            Interaction Reaction = moduleDef.createInteraction(r.getecnum(),InteractionBiochemicalRxnuri);
            for (Chemical c : r.getSubstrates()){
                if (natives.contains(c.getName())&& !visitedChem.contains(c)){
                    ComponentDefinition sideReactant = visitedComponents.get(c.getName());
                    FunctionalComponent sideReactantFC = moduleDef.createFunctionalComponent(c.toString(), AccessType.PUBLIC, sideReactant.getPersistentIdentity(),DirectionType.NONE);
                    Participation participation = Reaction.createParticipation(c.getName(),sideReactantFC.getPersistentIdentity(),ParticipationSideSubstrateUri);
                    continue;
                }
                if (!visitedChem.contains(c) && !natives.contains(c.getName())){
                    ComponentDefinition genericReactant= document.createComponentDefinition(c.getName(), BIOPAX_SmallMolecule_uri);
                    visitedComponents.put(c.getName(), genericReactant);
                    FunctionalComponent reactant = moduleDef.createFunctionalComponent(c.toString(), AccessType.PUBLIC, genericReactant.getPersistentIdentity(),DirectionType.NONE);
                    Participation participation = Reaction.createParticipation(c.getName(),reactant.getPersistentIdentity(),ParticipationReactantUri);
                    visitedChem.add(c);
                    Sequence chemSmile =  document.createSequence(c.getName(),c.getSmiles(),SequenceSmilesURI);
                    continue;
                }
                if (visitedChem.contains(c) && !natives.contains(c.getName())){
                    ComponentDefinition genericReactant = visitedComponents.get(c.getName());
                    FunctionalComponent reactant = moduleDef.createFunctionalComponent(c.toString(), AccessType.PUBLIC, genericReactant.getPersistentIdentity(),DirectionType.NONE);
                    Participation participation = Reaction.createParticipation(c.getName(),reactant.getPersistentIdentity(),ParticipationReactantUri);
                    continue;
                }
            }
            for (Chemical c : r.getProducts()){
                if (!visitedChem.contains(c) && !natives.contains(c.getName())){
                    ComponentDefinition genericProduct= document.createComponentDefinition(c.getName(), BIOPAX_SmallMolecule_uri);
                    visitedComponents.put(c.getName(), genericProduct);
                    FunctionalComponent product = moduleDef.createFunctionalComponent(c.toString(), AccessType.PUBLIC, genericProduct.getPersistentIdentity(),DirectionType.NONE);
                    Participation participation = Reaction.createParticipation(c.getName(),product.getPersistentIdentity(),ParticipationReactantUri);
                    visitedChem.add(c);
                    Sequence chemSmile =  document.createSequence(c.getName(),c.getSmiles(),SequenceSmilesURI);
                    continue;
                }
                if (visitedChem.contains(c) && !natives.contains(c.getName())){
                    ComponentDefinition genericProduct= visitedComponents.get(c.getName());
                    FunctionalComponent product= moduleDef.createFunctionalComponent(c.toString(), AccessType.PUBLIC, genericProduct.getPersistentIdentity(),DirectionType.NONE);
                    Participation participation = Reaction.createParticipation(c.getName(),product.getPersistentIdentity(),ParticipationProductUri);
                    continue;
                }
            }
            step = step + 1;
        }
    }

    public void out(String filename) throws SBOLConversionException, IOException {
        document.write(filename);
    }

}
