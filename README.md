## Optimal-Translate
For each sequence in a list of sequences, optimal_translate will translate in the reading frame with the least number of stop codons.

_Includes python script and geneious plugin extension._


**Usage**: 
```shell 
python optimal_translate.py [input.fasta] [output.fasta] 
```
Output.fasta will contain each sequence transated in the frame which produces the least number of stop codons. The frame will be appended to the sequence name. In cases where there are multiple optimum frames, all frames will be written to the output file. 

**Geneious Plugin Installer:** Install the function as a geneious plugin on windows. The plugin will appear as 'Optimal Translate' in the 'Tools' menu on Geneious. When running for the first time running you will need to specify the location of the .exe file. It will be present in your geneious plugins folder which is in C:\ProgramData\Geneious\plugins on Windows 7 or C:\Users\AppData\Local\Geneious on Windows 10. Requires Geneious 11.0.2.
