# Needleman-Wunsch / Smith Waterman algorithms
## Overview
This program provides the implementation of two of the most important 
algorithms in bioinformatics. Both are dynamic programming based
algorithms confronting the sequence alignment problem. That is, their 
goal is to find the best agreement between the sequences based on 
the statistics of probability of changes in the sequences, provided  
in the form of substitution matrix and gap penalty value.
The older Needleman-Wunsch algorithm is used to find the optimal global
alignment. That is, the best score across the entirety of both sequences, 
whereas the Smith-Waterman algorithm finds the optimal local alignment, 
parts of both sequences that have the most in common. 
## Usage and Input Format
It can be seen from this brief introduction that a lot of data is required
to make the procedure run properly, which may make passing the arguments
a bit impractical. 
That is the main reason why the basic run of the program will expect its
arguments mostly in file form. To be precise, it expects four or five arguments,
where the first is the path to the first sequence, second is the path to the
second, third is path to the file with the substitution matrix and finally
one or two numerical arguments are for the gap penalty value. 
The reason why there are two options where it comes to the gap penalty is that
while the function of the algorithms is usually demonstrated with a simple 
linear gap, practical use favors an affine gap penalty, because experiments 
show that it is far more likely to have few larger insertions or deletion, 
then to have many smaller. 
Consequently, the first number will be read as a linear gap value  and 
the optional second one as an affine one. 
Consider that the value is meant to be subtracted, therefore it is logical 
for it to be positive, but since technically there is no such rule it is 
up to the user. 
Both of the files with the sequences are expected to be in the simple and common
FASTA format, where lines starting with > are considered comments, while the 
rest forms the sequence itself. 
The matrix file is expected to provide first row of symbols separated by 
arbitrary number of whitespaces. All of the symbols are expected to appear in 
the first column, thus creating a table in which numerical values can be found. 
The 0,0 position of the table is not useful and so may be left blank. 
If this description doesn't provide sufficient clarity, please consult 
the test resource files with some examples. 
## Operation Outline
After passing the anticipated amount of arguments, the user is then prompted 
to choose between the two algorithms and is consequently presented with their 
results in the form of possibly multiple alignment options. The "-" symbol 
commonly used to denote a gap is replaced by "_", but the meaning stays 
the same. 
If no arguments are provided, the program will instead commence an interactive 
session where the user is expected to input the necessary information manually, 
albeit in a simplified form. Before that, there is an option to see the help 
string, which is presented on every call of the program with other than expected 
argument count, so it appears also in the case of zero arguments. 
If the user then chooses to run the interactive mode itself, the program will 
first ask for first and then second sequence, counting for them everything 
provided on the standard input. 
It would be very impractical to type the whole matrix so the user is instead
asked to successively provide values for match and mismatch between any two
symbols. It might be beneficial to notice that the mismatch value should
most likely be negative. Similar simplification goes for the gap penalty value, 
where only the linear option with one number remains. Ultimately, the user is 
asked to pick the desired algorithm and is consequently provided with 
the result in the same format as in the regular run. 
The package also provides four additional methods usable to compute 
the alignments for further uses elsewhere; they are named smithWaterman 
and needlemanWunsch respectively, both with two overloads for the two 
gap penalty options. 
