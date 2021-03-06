## TESSEE

<b>TESSEE - Tool for Early Stem cellS Economic Evaluation</b>

Stem Cell Engineering Research Group, Instituto Superior Tecnico, Universidade de Lisboa, Portugal

## Abstract:

Stem cell therapies are promising for diverse clinical indications. However, there are manufacturing and reimbursement challenges that must be addressed towards widespread adoption. This thesis presents TESSEE, a new tool for early health technology assessment (eHTA), supported on bioprocess and/or health economics models. TESSEE is developed specifically for stem cell therapies, incorporating biological, process, clinical, and economic uncertainty. Unlike other eHTA, TESSEE is an open source tool, freely available and customizable to several case studies. In order to develop and demonstrate the different features of TESSEE, at the moment, three industrially and clinically relevant case studies were published in academic journals.

## Pre-Installation Requirements

TESSEE is a cross platform tool that requires a Python 3 implementation and a text editor to run. The Anaconda implementation of Python 3 is recommended. For a Text Editor, Sublime Text 3 is also recommended.

Steps:

1 - Install [Anaconda Python 3](https://www.anaconda.com/distribution/) using all the default options.
       
 NOTE: In case you are running Anaconda on a Windows operating system, it is required to change the Path variables to add Anaconda manually:
 
 1 - Find the Anaconda installation path. Open the Anaconda Prompt and type:
       
       where conda
       
 2 - Copy the path to the files (ends in Anaconda3)
 
 3 - Open the command window/Terminal and type the following command. Click Enter. Where it reads $Anaconda Path$, replace by the path copied in the previous step
 
       setx path "%path%;$Anaconda Path$;$Anaconda Path$/Scripts"
       
       
 4 - Close and reopen the Command Window/Terminal
    
2 - Install [Sublime Text 3](https://www.sublimetext.com/3)

3 - Open the Command Window/Terminal in your operating system. In order to confirm if Python is well installed, run:

        python
        
in the terminal command window.

4 - Insert the following commands on your Command Window/Terminal. Each command should be followed by Enter


        conda update conda
        conda update numpy
        conda update pandas
        conda update scipy
        python -m pip install simpy


## Installation

There are two ways to install/download the source code of TESSEE to your local machine:

<b>Using GitHub (requires [pre-installation of Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git) in your machine):</b>

1 - Open the Command Window/Terminal in your operating system
    
2 - Change the directory to download the TESSEE source code on your command window:
    

    $cd <Insert path to desired location>

        
3 - Clone the Git repository 


    $ git clone https://github.com/catiabandeiras/TESSEE.git


<b>Direct download as a ZIP folder</b>

1 - Click on the [root GitHub repository of TESSEE](https://github.com/catiabandeiras/TESSEE)

2 - Click on the green "Clone or download" button

3 - Select "Download ZIP"

4 - Unzip the folder and copy to your desired location on your local machine.


## Demo Cases

For first feedback, please refer to the folder "demo1", where a simple demonstration case for production of doses of 10 million cells, using automatically selected planar technologies, is presented. In this case, only bioprocess modeling is provided. This version showcases automatically generated CSV files and plots for the total costs per donor/batch, total costs per dose, the costs division per stage, and the costs division per resource.

<b>Run the demo on the Command Line Interface (CLI)</b>

1 - Open the Command Window/Terminal in your operating system

2 - Change to the directory where you downloaded the "demo1" folder to:

    
    $cd <Insert path to desired location>
   

3 - Run the main file of the source code as follows:

    
    $ python -m main
    

4 - When the file stops, open the spreadsheets with results on the "csv_files" subdirectory and the plots on the "plots" subdirectory.

## Credits

Credits in the development of this tool are given to the PhD supervisors of [Catia Bandeiras](http://scerg.tecnico.ulisboa.pt/cbandeiras.html): [Prof. Dr. Frederico Ferreira (iBB-IST)](http://scerg.tecnico.ulisboa.pt/fcferreira.html) and [Prof. Dr. Stan Neil Finkelstein (IDSS-MIT and DCI-BIDMC-HMS)](https://idss.mit.edu/staff/stan-finkelstein/). 

Additionally, members of the [Stem Cell Engineering Research Group](http://scerg.tecnico.ulisboa.pt/index.html) provided important feedback on the development of this tool due to their lab expertise in Stem Cell Manufacturing. Furthermore, members of the Cell Therapy Core Facilities at [Dana-Farber Cancer Institute](https://www.dana-farber.org/), [Harvard Stem Cell Institute](https://hsci.harvard.edu/), and [Case Western Reserve University](https://case.edu/) also provided feedback, along with several industry members.

Invaluable feedback on cost-effectiveness models on type 1 diabetes and cystic fibrosis was provided by clinicians at the [Joslin Diabetes Center](https://www.joslin.org/), [Boston Children's Hospital](http://www.childrenshospital.org/), and [Indiana University School of Medicine](https://medicine.iu.edu/).

## License

GPL © [Stem Cell Engineering Research Group](http://scerg.tecnico.ulisboa.pt)
