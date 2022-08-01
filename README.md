# Welcome to SynPath!
An integrated platform for the prediction and evaluation of novel biosynthetic pathways for the production of target metabolites (SBOL compatible). The program is comprised of two parts: pathway discovery and pathway evaluation. Pathways are discovered via a retrosynthesis algorithm from the Metacyc database. Pathways are evaluated by inserting them into user-provide genome-scale models and analying them using COBRA methods, such as flux balance analysis and flux variability analysis. The user defines a target metabolite, maxium number of pathway steps, any precursor and a chassis model. Pathways will be evaluated based on several parameters: theoreticla yield (aerobic vs anaerobic), energy use (aerobic vs anaerobic), cofactor balance (aerobic vs anaerobic), and thermodynamic feasibility. The users also have the option to have the results saved as SBOL-compatible .xml files.

## Software requirements
The Docker image will run on any OS that Docker supports (e.g. macOS, Linux, Windows).
### Docker
Docker is the preferred environment, as it creates reproduceable runs in a tested runtime environment. Docker also avoids the installation headaches and potential pitfalls of directly installing applications into your system. Install docker here: https://docs.docker.com/get-docker/

## Running SynPath
There are several ways of running SynPath.
1. Navigate to your terminal, and pull the project from docker registry using the command `docker pull carol111/sbolmetpathdesign`, followed by `docker run -dp 5000:5000 carol111/sbolmetpathdesign`. After the program starts running, open any browser and navigate to the link `http://127.0.0.1:5000` to see the interface.
2. First clone the package: `git clone https://github.com/yuzig/SBOLmetPathDesign.git`. Navigate to the docker directory within the project by `cd SBOLmetPathDesign/docker`. Build and create a virtual environment within docker by using `make build`, and this should take several minutes. Enter the docker container by using this command `make dev`. Navigate to the directory where you saved this package to, and enter the `SBOLmetPathDesign/cobrapyconverter/GenomeScaleModels`. run `python main.py` to start running. Open `http://127.0.0.1:5000` in a broswer to see the interface.

## Screenshot of the workflow
Figure 1. Overview of the SynPath workflow
![overview](/demo_screenshots/overview.jpg)
Figure 2.  screenshot of home page
![demo1](/demo_screenshots/home_page.jpg)
Figure 3. screenshot of result page, with Lycopene, CPD1F-114, as in put, FARNESYL-PP as precursor, maximum step set to 7, and with the E.coli model iML1515
![demo2](/demo_screenshots/result_page.jpg )
