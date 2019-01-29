# OPT4e
Optimal Prediction for Type4 secretion system effector proteins

(Important please make sure to download first and second part of OPT4e sofware and put them all in one folder.) 

To run this software you need to have python and Emboss installed:

Installation:

-For python it is highly recommended to install Jupyter Notebook and run the OPT4e notebook. For this, you can use this link and https://www.anaconda.com/download/ and install anaconda for python 3 which automatically adds Jupyter Notebook to your operating system. 
Although you can just have python installed on your system and run use the OPT4e.py script. 

-Go to "TO INSTALL" folder and install EMBOSS using "mEMBOSS-6.5.0.0-setup" which is very straightforward. Although, you need to have java installed on your system. In case you do not have it, you can use "JavaSetup8u191" from the same folder to install java.

Run:
-Before running, you need to provide your input protein sequence(s). The input file should has ".fa" format and should be named "input_sequences". It is recommended to copy and paste your fasta format sequences in the "input_sequences.fa" file which is provided in the folder. 
-To run the software, simply run OPT4e notebook or OPT4e.py script and the GUI will open up giving you some options to choose based on what you want to do. 
--optional but recommended: to make sure that you have all necessary packages in your system, you can first run "packages notebook" or "packages.py" script
--please note that after running the OPT4e notebook, a toggle button will show up t the end of the code that you can use to show/hide the code when using the software.

-When the GUI pops up, you have 4 buttons to choose. 
1-"Predict effector candidates for input sequences": this button is used to predict of you input protein sequence(s) is a candidate T4SS effector. 
Please not that this step might take long to get completed based on the number of input sequences. Please do not close the GUI window. When the task is done, a new window will open and let you know that it is done. 
The predicted candidate effectors will be saved in 4 files (in the same folder of the software package), based on their likelihood of being an effector: predicted_effectors_mostlikely, predicted_effectors_probable, predicted_effectors_lesslikely, predicted_effectors_total. 
The number of predicted candidates are reported at the beginning of each file. 

2-"Add new validated effectors to training set": in case you have some experimentally verified effectors and want to add them to the software to make its future predictions better, you can use this button. 

3-"Add new validated non-effectors to training set": in case you have some known non-effectors and want to add them to the software to make its future predictions better, you can use this button. 

4-"Reset training set to original version": In case you realize that the changes you made to the training set, using the second or third button, were not correct and want to undo it, you can use this button. 

Thank you for using OPT4e.
Zhila Esna Ashari
z.esnaashariesfahan@wsu.edu
