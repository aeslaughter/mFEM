mFEM Toolbox



MAE4700/5700 Homework Assignments

Please read these intructions carefully if you are adding or editting the
homework assignments.

LaTeX File Structure
All of the LaTeX source code is located in the assignments/latex directory,
with each assignment located in a folder within this directory.

To add a new assignment follow these steps:
(1) Create a new folder for storing the LaTeX source (e.g., HW12), it is a
good idea to use the number in the name because the building of the
assignment for release to the students is done by referencing the 
assignment number (see (2) below)
(2) Add the folder to the list in assignments.tex (e.g. \assignment{HW12}), 
the assignment that is added is refered to later by the assignment number
which is determined by the order of the folders listed in the 
assignments.tex file.
(3) Within the new folder create a tex file with the same name as the
folder (e.g., ../latex/HW12/HW12.tex). This file contains the text for the 
homework problem and it should contain a \chapter command, see the other
problems for example.
(4) Questions within each assigment are added using the question command.
For example, \question{HW12-a} automatically adds the file located in 
../latex/HW12/H12-a/prob.tex. If a figure for this problem is desired a 
fig.pdf file should be included in the same location. The solution may also 
be added with a soln.tex file. If the names of the files for each question 
follow this convention they will be added automatically.

buildHW.m
In the main directory, the buildHW.m function is used to update the MATLAB
source code, build the LaTeX of the homework and solution and package the
assignment in a zip file for distribution.
    


