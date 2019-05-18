[![Language Demo Video](https://github.com/rafo23/biolpp/blob/master/img/Sequence%2001.00_00_00_01.Still002.jpg)](http://www.youtube.com/watch?v=1ArjKt102Ao "Video")

## Introduction

Computational biology is a broad and expansive field of study centered mostly on the storage, processing and analysis of
biological data pertinent to many fields such as genetics, molecular biology and biophysics, among others. It uses applied
mathematics and statistics to aid in the analysis of large amounts of data, often requiring a way to output the results onto
graphs and other displays to aid in the visualization of the data and results of the analyses. The field is intricately linked
to computer science as it makes heavy use of computers to acquire, store and process all the data, as well as designing
algorithms to process it and analyze it.

The main motivation comes from the fact that many of the tasks performed by computational biologists are repetitive and
require several steps to achieve the final result. The goal of this language is to streamline the workflow and condense many of
these tasks into one or two simple instructions which will automatically output the results as desired. Another reason behind
creating such a language is that it allows biologists with little to no programming knowledge or experience to easily analyze
their own data without having to rely on computational biologists to do it for them. This way they can benefit immensely without
having to spend months trying to learn how to program or learn a complex programming language. Using simple and straightforward
functions, anyone can learn how to use this programming language in just one afternoon.


## Installation

In order to work with Biol++ the following elements must be met:
- [Python 3.0](https://www.python.org/downloads/) or above - (latest - [3.7.3](https://www.python.org/downloads/release/python-373/) as of writing this document)
- A Python virtual environment and/or a development environment (i.e. [Pycharm](https://www.jetbrains.com/pycharm/))
- Install [Pip](https://pip.pypa.io/en/stable/installing/)
- Install [Biopython](https://biopython.org/wiki/Documentation)
- Install Pylab using pip 

## Language Features
- Simple and effective syntax - Easy to learn and implement
- DNA/RNA translation
- Transcription/Complementation of DNA/RNA strands
- Basic recurrence relationship functions
- Strands Analysis (GC-Content, Hamming distances in mutations, motif intervals)
- Simple Mendelian Inheritance Probabilities
- FASTA format file analysis supported
- Protein translation and inferring
- Opening Reading Frames (ORFs) analysis
- Phylogenetic Tree Generation
- Protein inference & weight
- Punnet diagrams generation
- Others

## Sample Program
```
      BIOL++ >>> print('IOFiles/file.txt')
      BIOL++ >>> a = seq('GATGGAACTTGACTACGTAAATT')
      BIOL++ >>> a
      BIOL++ >>> b = comp(a)
      BIOL++ >>> b = transc(a)
      BIOL++ >>> b = transl(a, dna)
      BIOL++ >>> b = read('Seq 1', 'IOFiles/file.txt')
      BIOL++ >>> write(a, 'IOFiles/myfile')
      BIOL++ >>> punnett('Aa Bb', 'Cc Dd')
```
## Language Development 
### Translator Architecture 

The language was built using the Python programming language as well as the Python library PLY, which combines both a lexer and a parser. The language runs directly in the terminal where commands can be entered and executed immediately, with no need to write out the code first and then compile and run. This makes the language more versatile and easier to use and test. The following figure shows the basic structure of the language.

![Image](https://github.com/rafo23/biolpp/blob/master/img/PL%20Sketch.png)

      Figure 1 Translator Architecture for BioL++ 
      
As seen in Figure 1, the user inputs the BioL++ code into the terminal where it is passed on to the lexer. This generates the different tokens which are then passed on to the parser which matches the correct grammar rule. If no errors are found in the syntax, the underlying intermediate code in Python is executed to produce the requested output, which is then displayed to the user in the terminal. 

### Module Interfacing 

The entire project was developed with modularity in mind. All aspects are separated based on their functionality. The main program is responsible for instantiating the parser and the loop to retrieve user input from the terminal window. The input is sent to the parser for analysis. Inside the parser, the lexer module is referenced and the token list defined there is obtained. Based on the syntax, grammar rules and the specific function that was input, the parser either throws a syntax error to the user, or on a syntactically correct statement calls the respective algorithm in the algorithms module to process the data and produce the desired output, which is then displayed in the terminal or saved to a file on the computer.  

![Image](https://github.com/rafo23/biolpp/blob/master/img/PL%20Diag.png)

      Figure 2 Module interfacing 
      
The algorithms module contains our own implementations as well as external third-party libraries. The main library is [Biopython](https://biopython.org/wiki/Documentation), which is used in combination with [pylab](https://scipy.github.io/old-wiki/pages/PyLab) to read and draw the phylogenetic trees. This library allows us to draw an ASCII representation of the tree directly in the terminal as well as an actual image which can be modified and saved to the machine using the [matplotlib](https://matplotlib.org/) viewer. 

### Software Development Environment 

The entire project development was done with the help of the following programs/applications: JetBrains PyCharm, PIP/virtualenv and GitHub Desktop. PyCharm is a free integrated development environment (IDE) tailored for Python development and contains a plethora of tools and aids to help with coding in Python. It provides tools such as syntax highlighting and code completion, as well as project management and debugging tools. PIP is Python’s own package manager that allows one to download and install packages and libraries seamlessly to one’s machine through the terminal window. Virtualenv allows one to create a separate Python environment for development which only contains the necessary packages and the required versions to properly run the application. The final component is GitHub Desktop which greatly simplifies teamwork and version control. It allows groups of developers to push their work and retrieve updates done by others and manage the work progress. 

### Test Methodology 

A separate test program was written to make sure all the algorithms worked as they were supposed to before actually getting into the language itself. The language on its own can be easily tested by running the commands directly on the terminal and comparing the output with the expected result obtained manually or from a trusted source. The test program used to verify the algorithms separately before being incorporated into the project can be found in the [tester.py](https://github.com/rafo23/biolpp/blob/master/tester.py) file

The lexer was tested by hardcoding different expressions and verifying the console output to make sure the lexer correctly tokenized the expression. Below is the sample test program which is located in the lexer module itself. 

```
Tester

lexer.input("   bseq abc = btree .read()")

while True:
     tok = lexer.token()
     if not tok:
         break
     print(tok)

```
Once these are tested and proven to work, the next step is to test the entire language in the terminal since there is no way to test the parser separately and at this point the entire structure of the language is already done. As was mentioned previously, the test methodology consisted of running the commands in the terminal for which we already knew the result. We compared the output and made sure it matched. Trying to break the language is the best way to test it and make sure it is robust. Throwing in syntax errors and mismatched data types helps to confirm it is working in order.

## Reference Manual

The [reference manual]() contains the list of all the functions the language contains, as well as a detailed guide on how to use these.
