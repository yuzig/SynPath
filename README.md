# SBOLmetPathDesign
An integrated platform for the prediction and evaluation of novel biosynthetic pathways for the production of target metabolites (SBOL compatible). The program is comprised of two parts: pathway discovery and pathway evaluation. Pathways are discovered via a retrosynthesis algorithm from the Metacyc database. Pathways are evaluated by inserting them into user-provide genome-scale models and analying them using COBRA methods, such as flux balance analysis and flux variability analysis. The user defines a target metabolite, maxium number of pathway steps, any precursor and a chassis model. Pathways will be evaluated based on several parameters: theoreticla yield (aerobic vs anaerobic), energy use (aerobic vs anaerobic), cofactor balance (aerobic vs anaerobic), and thermodynamic feasibility. The users also have the option to have the results saved as SBOL-compatible .xml files.

## Software requirements
The Docker image will run on any OS that Docker supports (e.g. macOS, Linux, Windows).
### Docker
Docker is the preferred environment, as it creates reproduceable runs in a tested runtime environment. Docker also avoids the installation headaches and potential pitfalls of directly installing applications into your system. Install docker here: https://docs.docker.com/get-docker/

## Running MetPathDesign
1. First clone the package: `git clone https://github.com/yuzig/SBOLmetPathDesign.git`
2. Navigate to the docker directory within the project by `cd SBOLmetPathDesign/docker`
3. Build and create a virtual environment within docker by using `make build`, and this should take several minutes
4. Enter the docker container by using this command `make dev`
5. Navigate to the directory where you saved this package to, and enter the SBOLmetPathDesign/cobrapyconverter/GenomeScaleModels
6. run `python main.py' to start running an example using 1,3-PDO as target within the E.cli model iML1515. The results will be printed on the terminal and the SBOL-format files saved to SBOLmetPathDesign/cobrapyconverter/GenomeScaleModels/result

## Instructions for use
To run the software with other inputs, open SBOLmetPathDesign/cobrapyconverter/GenomeScaleModels/args.yml and change the inputs manually and then run `python main.py`.

To run with a customized model, add the model file to SBOLSBOLmetPathDesign/cobrapyconverter/GenomeScaleModels/bigg_models and change the model input in args.yml.
