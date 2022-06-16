# PairwiseSequenceAlignment
Implementation of local and global pairwise sequence alignment algorithms.

# Requirements
- Python(version >= 3)

# Usage
- python3 pairwise_sequence_alignment.py &#45;&#45;input &lt;path to input text file containing sequences&gt; &#45;&#45;alignment &lt;local or global&gt; &#45;&#45;scoring&#45;matrix &lt;path to scoring matrix&gt; &#45;&#45;gap&#45;opening&#45;penalty &lt;a negative number&gt; &#45;&#45;gap&#45;extension&#45;penalty &lt;a negative number&gt; &#45;&#45;output &lt;path to output file&gt;
- &#45;&#45;output is optional.
- The input text file must contain only two lines, each of which contains a sequence. An example for expected input file is sequences.txt.
- As long as it is in a similar form with the one in BLOSUM62.txt file, any scoring matrix (e.g., PAM110, BLOSUM50, etc.) can be given as input and used.
- The tool can also be used for DNA sequences provided that an appropriate scoring matrix is given as input.

# Sample Run
- Command:
    > python3 pairwise_sequence_alignment.py &#45;&#45;input sequences.txt &#45;&#45;alignment local &#45;&#45;scoring&#45;matrix BLOSUM62.txt &#45;&#45;gap&#45;opening&#45;penalty &#45;5 &#45;&#45;gap&#45;extension&#45;penalty &#45;2
- Output:


    > ![Screenshot](https://raw.githubusercontent.com/ender-s/PairwiseSequenceAlignment/main/ss.png)
