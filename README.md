# BridgeEvaluator
This package designs bridgeRNAs targeting a given locus and and provides metrics to evaluate their efficiencies and specificities

*__Author__: Jaymin Patel: jayman1466@gmail.com*

## Installation
### Streamlined Installation
You can install this package via pip:
```
pip install BridgeEvaluator
```
You must also install the [Vienna RNA Suite](https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/install.html) if you want to score the predicted folding of the designed bRNAs. If you don't have this installed, set *score_structure=False* in the *design_bridges()* command. 

### Manual Installation
Alternatively, for manual installation, you can place the files from the __"src/BridgeEvaluator/"__ directory directly into your working directory. If you use this manual installation, make sure you have the following dependencies installed with a __python version >=3.9__:


[biopython >= 1.85](https://pypi.org/project/biopython/)

[Levenshtein >= 0.27.1](https://pypi.org/project/Levenshtein/)

[viennarna >= 2.7.0](https://pypi.org/project/ViennaRNA/)

[pandas >= 2.2.0](https://pypi.org/project/pandas/)


## Usage
The simplest usage within a python script is as follows:
```python
from bridge_evaluator import design_bridges

#Specify the locus that will be scanned to design bRNAs
target_locus = "ATGAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCCGTGGAGAGGGTGAAGGTGATGCTACAAACGGAAAACTCACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCGTGGCCAACACTTGTCACTACTCTGACCTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCACATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGACCTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAGGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGACACAAACTCGAGTACAACTTTAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACGTTGAAGATGGTTCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGTCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAAtaa"

#Name of this locus. The results will be outputted with this filename
target_name = "sfGFP"

#Genbank file of the full genome of the target. This will be used to evaluate possible off targets 
genbank_file = "MG1655.gb"

design_bridges(target_locus, target_name, genbank_file)
```

### Required Parameters
__target_locus__: Nucleotide sequence of the locus you want to target. This script will scan this locus to identify all permissive 14mer target sites and provide attributes to score their predicted efficiency and specificity.

__target_name__: Name of the locus you are targeting. This will be used to name the designed bridgeRNAs and the output file of the script. bridgeRNA names are given the syntax: *bridge_IS621_T_{target_name}_{index}_D_{donor_name}*. The {target_name} and {donor_name} are pulled from the arguments of this function.


__genbank_file__: Genbank file of the recipient genome. This will be used to identify potential off target sites for each designed bridgeRNA.


### Optional Parameters
__donor_seq__: Sequence of the Donor Sequence (14mer) being used. Default is the native IS621 donor __"ACAGTATCTTGTAT"__

__donor_name__: Name of the Donor. This will be used to name the designed bridgeRNAs. Default is __"1"__

__cores__: Core sequences that can be used, provided as a list. Note, the core of the provided donor sequence will be modified to match the core of the target sequence for each designed bridgeRNA. Default is __['CT']__

__kmer__: The first X bp of the target sequence that will be used to identify perfect and imperfect offtargets in the recipient genome. Default is __11__

__avoid_restriction__: Restriction sites (and other sequences) to avoid in the designed bridgeRNAs, provided as a list. Reverse complements must be provided manually. Default is __[]__

__check_imperfect__: Should this script look for imperfect offtarget sequences in the recipient genome? This increases computational time dramatically. All imperfect offtargets with a Levenshtein distance <= 2 are tabulated. Indels are given a Levenshtein score of 2 and mismatches are given a Levenshtein score of 1. Default is 
__True__

__score_structure__: Should this script evaluate the predicted secondary structure of the designed bridgeRNAs? Predicted MFE structures of designed bridgeRNAs are compared to the reference native IS621 secondary structure using Vienna RNA's RNAforester forest-based structural aligner. A similarity score from 0 to 1 is provided. Default is __True__

__feature_type__: For internal metadata, you can include a feature type (eg CDS, ncRNA) for the locus you are targeting. Default is __""__

__primer_seqs__: In addition to full bridgeRNA sequences, this script converts these sequences into DNA fragments for synthesis which can be cloned into the bsaI sites in *Patel et al.* IS110 vectors. If you'd like to include PCR priming sites to amplify these fragments, you can specify them as a dictionary of lists, segmented by core sequence, as follows: <span style="color: red;">{'CT': ["for_primer_seq1", "rev_primer_seq1"], 'GT': ["for_primer_seq2", "rev_primer_seq1"]}</span>. Default is __{"CT": ["",""], "GT": ["",""], "AT": ["",""], "TT": ["",""]}__   

### Example python script utilizing some optional parameters:
```python
from bridge_evaluator import design_bridges

#Specify the locus that will be scanned to design bRNAs
target_locus = "ATGAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCCGTGGAGAGGGTGAAGGTGATGCTACAAACGGAAAACTCACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCGTGGCCAACACTTGTCACTACTCTGACCTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCACATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGACCTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAGGGTATTGATTTTAAAGAAGATGGAAACATTCTTGGACACAAACTCGAGTACAACTTTAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACGTTGAAGATGGTTCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGTCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAAtaa"

#Name of this locus. The results will be outputted with this filename
target_name = "sfGFP"

#Genbank file of the full genome of the target. This will be used to evaluate possible off targets 
genbank_file = "MG1655.gb"

#Donor A is the native Donor of IS621
donor_seq = 'ACAGTATCTTGTAT' #donor A

#Cores to use
cores = ["CT","GT"]

#avoid the following sequences - manually include the reverse complement
avoid_restriction = ["GGTCTC","GAGACC","GCTCTTC","GAAGAGC"]

design_bridges(target_locus, target_name, genbank_file, donor_seq = donor_seq, cores = cores, avoid_restriction = avoid_restriction)
```

## Output

The results are exported in the present directory as a csv file named __{target_name}.csv__. Each row represents a potential bridgeRNA that can target the provided locus

### Output Columns: 

__target_gene__: Name of locus being targeted

__feature_type__: Feature type of locus being targeted

__donor_seq__: Donor Sequence (14mer)

__target_seq__: Target Sequence (14mer)

__index__: Position of the Target Sequence relative to target locus

__strand__: Orientation of the Target Sequence (+ or -) relative to the target locus

__core__: Core sequence being used

__perfect_match_targets__: Number of perfect matches of the Target Sequence present in the provided recipient genome. Note, only the first X bp specified by the kmer attribute (default = 11) is used for matching.

__levenshtein_distance_1_targets__: Number of matches of the Target Sequence present in the provided recipient genome with a Levenshtein distance of 1. Note, only the first X bp specified by the kmer attribute (default = 11) is used for matching. Indels are scored as a distance of 2. SNPs are scored as a distance of 1.

__levenshtein_distance_2_targets__: Number of matches of the Target Sequence present in the provided recipient genome with a Levenshtein distance of 2. Note, only the first X bp specified by the kmer attribute (default = 11) is used for matching. Indels are scored as a distance of 2. SNPs are scored as a distance of 1.

__bridge_sequence__: Full sequence of the bridgeRNA

__p6p7_warning__: Warning if this bridgeRNA violates preferred handshake rules.

__RNA_structural_similarity__: The predicted MFE secondary structure of designed bridgeRNA is compared to the reference native IS621 secondary structure using Vienna RNA's RNAforester forest-based structural aligner. A similarity score from 0 to 1 is provided. This can be used to assess whether this bridgeRNA is likely to misfold. 

__eblock_seq__: The bridgeRNA converted into a DNA fragment for synthesis, which can be cloned into the bsaI sites in *Patel et al.* IS110 vectors.

## References 
[Bridge RNAs direct programmable recombination of target and donor DNA](https://www.nature.com/articles/s41586-024-07552-4)<br>
[Arc Institure Bridge RNA Design Tool](https://bridge.hsulab.arcinstitute.org/)