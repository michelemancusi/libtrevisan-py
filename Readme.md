# Trevisan-extractor
TREVISANEXTRACTOR.py is the python code for the extractor, I thank Robert Campbell (campbell@math.umbc.edu) http://www.math.umbc.edu/~campbell/ from whom I took the necessary functions to be able to do the calculations in finite fields.

source.txt and seed.txt are two sample files because the code requests in input a source file and a seed file. Obviously these two files will have to be replaced with your own source and seed files (always in the .txt extension containing zeros and ones).

When TREVISANEXTRACTOR.py is executed, it requires input data: the name of the source file, the name of the seed file, the value of the min-entropy per bit of the source and finally the desired error per bit (the smaller the error the more uniform the output string will be, but the longer the execution time will be).

P.S. Unlike the librevisan repository (written in C code), in this python version, only the SWD (standard weak design) is implemented and not also the BWD (block weak design), see Theoretical_guidelines_for_Trevisan_extractor.pdf for more details. 
