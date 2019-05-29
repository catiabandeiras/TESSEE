## TESSEE

TESSEE - Tool for Early Stem cellS Economic Evaluation

Stem Cell Engineering Research Group, Instituto Superior Tecnico, Universidade de Lisboa, Portugal

Check out the Read the Docs page for more information on the project:


Abstract:

Stem cell therapies are promising for diverse clinical indications. However, there are manufacturing and reimbursement challenges that must be addressed towards widespread adoption. This thesis presents TESSEE, a new tool for early health technology assessment (eHTA), supported on bioprocess and/or health economics models. TESSEE is developed specifically for stem cell therapies, incorporating biological, process, clinical, and economic uncertainty. Unlike other eHTA, TESSEE is an open source tool, freely available and customizable to several case studies. In order to develop and demonstrate the different features of TESSEE, at the moment, three industrially and clinically relevant case studies were published in academic journals.

Introduction:

Stem cell based therapies may be a breakthrough for several unmet medical needs. Their efficacy
has already been proven for graft vs host disease, osteoarthritis, acute myocardial infarction and diabetic retinopathy, and clinical trials on several other prospective indications on the field of neurological diseases, diabetes and autoimmune diseases are also being explored. The global market for cell based therapies currently generates annual profits of more than $ 1 billion, with an estimated revenue of $ 20 billion in 2025. In particular, stem cells have regenerative and immunomodulatory potential to address a diverse number of unmet medical needs. Over 5400 clinical trials related to stem cells as an intervention have been reported until now, with the majority of the trials being related with adult stem cells, like the hematopoietic (1763 trials) and mesenchymal (811). Pluripotent stem cells, like the induced pluripotent (50 trials) and embryonic stem cells (34 trials) are in an earlier stage of development; despite the large interest, regulatory approval for these therapies has been difficult. However, currently, there are six approved products in specific countries and three reimbursed products, with the price of one course of therapy rising up to dozens or even hundreds of thousands of dollars.

The widespread application of stem cell based therapies would benefit from reducing reimbursement price, while maintaining product profitability. Moreover, the set reimbursement price must cover the research and development and clinical trials costs, but also the manufacturing costs of such therapies, which are still extremely high when compared with conventional pharma or biotherapeutic products.
These large costs are due to largely manual product handling and manipulation, product and process variability, impractical scaling-up of production, use of xenogeneic materials, high culture media costs, and high costs of quality control. Commonly used small scale planar expansion platforms, with cells cultivated in 2D surfaces, such as T-flasks, are not enough to meet market demands and ensuring maintenance of the therapeutic potential of the product. Apart from difficult scaling-up, they do not allow control and monitoring of culture parameters, lead to development of concentration gradients and require a lot of incubation space and manual operation. Therefore, other manufacturing methods need to be adopted in order to provide more cost competitive therapies and with higher possibilities of being lucrative upon the thresholds for reimbursement that the payers from several countries impose on the therapies.

Current process design is guided by the envisioned demand and compliance with regulatory requirements. However, design of manufacturing processes streamlining for cost efficiency, while preparing a new therapy for approval and reimbursement, is often neglected. After initial regulatory approval, further process manufacturing changes are usually administrative and their validation is cost prohibitive. In order to have a more thorough and less time consuming risk assessment of changes in process design, computational modeling of bioprocessing and bioeconomics is of great value to consider the impact of those changes on the process costs and quality of the final products given biological and technological parameters informed by past experiments. Computational decision support tools can contribute to faster, safer and less expensive production of therapies. Namely, through design of logical processes and optimization of several manufacturing parameters to achieve the lowest cost of goods (CoG) for a given demand of doses and lots of the therapy, as well as providing recommendations on which unit operations have higher impact on the process costs and need to be further optimized. These tools can also allow to select, for a given demand of therapeutic doses, the production configuration that ensures manufacture profitability within and maximum reimbursement price thresholds accessible to relevant payers.

The area of computational modeling for stem cell manufacturing is a recent one and there are a few academic contributions in the field, either using commercial flowsheeting software, like Superpro Designer or on custom-made code. The published models are focused on either the simulation of bioprocessing of allogeneic mesenchymal stem cell based therapies or the simulation of manufacturing of induced pluripotent stem cell derived differentiated cells, such as cardiomyocytes, neurons and pancreatic progenitors. In terms of manufacturing challenges, the limitations of 2D culture systems were addressed in terms of cost of goods of expansion and the inability to meet high dose demands for doses containing high numbers of cells. The process strategies for overcoming these limitations, include the automation for processing of large multi stack systems as opposed to regular T-flasks, and the use of stirred tank bioreactors with microcarriers to increase expansion area. Process modeling allows to evaluate cost effectiveness of stem cell production in suspension over planar technologies, as well as the selection of the best combinations between upstream and downstream technologies. While most of these studies focused on deterministic parameters and employed sensitivity analysis to determine the impact of key model parameters on final costs, stochasticity was also considered with appropriate statistical distributions through the Monte Carlo method. These tools can be employed for a multitude of manufacturing problems and ultimately propose technological and pricing changes in materials employed in cell manufacturing.

Pre-Installation Requirements

TESSEE is a cross platform tool that requires a Python 3 implementation and a text editor to run. The Anaconda implementation of Python 3 is recommended. For a Text Editor, Sublime Text 3 is also recommended.

Steps:
    1 - Install Anaconda Python 3: https://www.anaconda.com/distribution/
        NOTE: In case you are running Anaconda on a Windows operating system, it is required to change the Path variables to add Anaconda: https://medium.com/@GalarnykMichael/install-python-on-windows-anaconda-c63c7c3d1444
    2 - Install Sublime Text 3: https://www.sublimetext.com/3
    3 - Open the Command Window/Terminal in your operating system
    4 - Insert the following commands on your Command Window/Terminal. Each command should be followed by Enter
        $conda update conda
        $conda update numpy
        $conda update pandas
        $conda update scipy
        $pip install simpy

   
   
Installation

There are two ways to install/download the source code of TESSEE to your local machine:

1 - Using GitHub (requires pre-installation of Git in your machine):

    1 - Open the Command Window/Terminal in your operating system
    
    2 - Change the directory to download the TESSEE source code on your command window:
    
        $cd <Insert path to desired location>
        
    3 - Clone the Git repository 

