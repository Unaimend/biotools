# biotools

## Tetra
Calculates different variations of tetra-nucleotide frequencies of different input data

### Functions 
`countTetraNucleotide(sequence: string):CountTable[string]`  
Counts 4-mer amount of a given string. If a 4-mer does not not occur in the string, it will not occur in Table, to get a Table in which not occuring 4-mers have a value of 0 in the table use `countTetraNucleotideAll`  
Parameter 1: String of which to count the mount of 4-mers.  
Return: Table in which keys are the 4-mers and the corresponding values are the amounts of it.  


`countTetraNucleotideAll(sequence: string):CountTable[string]`  
Counts 4-mer amount of a given string. If a 4-mer does not not occur in the string, it will still occur in the table with a value of zero.  
Parameter 1: String of which to count the mount of 4-mers.  
Return: Table in which keys are the 4-mers and the corresponding values are the amounts of it. Not occuring 4-mers have a value of zero.  
