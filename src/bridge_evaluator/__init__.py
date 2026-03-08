from Bio.Seq import Seq
from Bio import SeqIO
from rapidfuzz import process, distance
import pandas as pd
import numpy as np
import heapq
import subprocess, textwrap
import RNA
from collections import Counter
from classes import SCAFFOLD_NAME_TO_CLASS

#author: Jaymin Patel jayman1466@gmail.com
VERSION = "1.0.0"
__all__ = ["design_bridges"]

#wrapper to design the IS621 bridgeRNA using Matt's script https://github.com/hsulab-arc/BridgeRNADesigner
def design_bridge_rna(target, donor, scaffold):
    target = target.upper()
    donor = donor.upper()
    brna_scaffold = SCAFFOLD_NAME_TO_CLASS[scaffold]

    # Hard checks, will raise errors if not met
    brna_scaffold.check_target_length(target)
    brna_scaffold.check_donor_length(donor)
    brna_scaffold.check_target_is_dna(target)
    brna_scaffold.check_donor_is_dna(donor)

    # Warning if not met
    brna_scaffold.check_core_mismatch(target, donor)
    p6p7_warning = brna_scaffold.check_p6p7_match(target, donor)

    bridge_design = brna_scaffold()
    bridge_design.update_target(target)
    bridge_design.update_donor(donor)
    bridge_design.update_hsg()
    TBL_seq, DBL_seq = bridge_design.generate_split_loops()

    return bridge_design, p6p7_warning, TBL_seq, DBL_seq



#function to design an eblock/DNA fragment for cloning via golden gate in bsaI sites
def eblock_design(input_seq,left_primer="",right_primer=""):
    
    left_remove = 'AGTGCAGAGAAAATCGGCCAGTTTTCTCTGCCTGCAGTCCGCATGCCGT'
    right_remove = 'TGGTTTCACT'
    left_gg = 'GAGAGggtctcTCCGT' #left bsaI site - keep lowercase to avoid the restriction screening
    right_gg = 'TGGTAgagaccGAGAG' #right bsaI site
    stuffer = "GACATTGTCCCTGATTTCTCCACTACTAATAGCACACACGGGGCAATACCAGCACAAGCTAGTCTCGCGGGAACGCTCGTCAGCATACGAAAGAGCTTAAGGCACGCCAATTCGCACTGTCAGGGTCACTTGGGTGTTTTGCACTACCGT" #150bp stuffer to get to 300bp fragment size

    
    #append additional PCR primers to ends for amplification of eBlock. This expects both primers to be in 5' to 3' orientation 
    left_OH = left_primer + left_gg
    right_OH = right_gg + str(Seq(right_primer).reverse_complement()) + stuffer

    mod_seq = input_seq.replace(left_remove,left_OH)
    mod_seq = mod_seq.replace(right_remove,right_OH)
    return mod_seq

#check for restriction sites in the eblock/DNA fragment
def check_restriction(avoid_restriction, eblock):
    avoid_restriction = [item.upper() for seq in avoid_restriction for item in (seq, str(Seq(seq).reverse_complement()))] #add in the reverse complements 
    avoid_restriction = set(avoid_restriction) #remove any duplicate entries

    total_restriction_sites = 0
    for restriction_site in avoid_restriction:
        total_restriction_sites = total_restriction_sites + eblock.count(restriction_site)

    return total_restriction_sites

#quick function to get all indexes of substring within a text
def find_all_loop(sub, text):
    indexes = []
    start = 0
    while True:
        start = text.find(sub, start)
        if start == -1: return indexes
        indexes.append(start)
        start += 1

#function to convert genome into a counter of kmers for efficient string lookup
def kmer_counter(genome_seq,kmer_length,core):
    core=core.upper()
    genome_seq_rc = str(Seq(genome_seq).reverse_complement())

    #only find the kmers that contain the relevant core sequence
    plus_indexes = find_all_loop(core, genome_seq)
    minus_indexes = find_all_loop(core, genome_seq_rc)

    #coordinates to pull
    left = 7
    right = kmer_length-left 


    #create list of all possible targets, based on presence of the core
    all_possible_targets = []
    for index in plus_indexes:
        if index>left and index<len(genome_seq)-right:
            all_possible_targets.append(genome_seq[index-left:index+right])
    for index in minus_indexes:
        if index>left and index<len(genome_seq_rc)-right:
            all_possible_targets.append(genome_seq_rc[index-left:index+right])

    #convert to a counter
    kmers = Counter(all_possible_targets)

    #free up memory
    del genome_seq_rc
    del plus_indexes
    del minus_indexes
    del all_possible_targets

    return kmers

#search for offtargets to the current target by Levenshtein distance
def off_target_analysis(targets:list, genome_seq:str,core:str, kmer_length:int, kmers:Counter, distances = [0,1,2]):
    
    #trim targets to length kmer_length
    targets = [target[:kmer_length] for target in targets]

    choices = list(kmers.keys()) #list of unique k-mers
    occurences = list(kmers.values()) #list of occurence of each kmer
    
    #calculate Levenshtein distance of all targets across all kmers 
    result_matrix = np.array(process.cdist(targets, choices, scorer=distance.Levenshtein.distance, score_cutoff=2, workers=-1))
    
    results = {} #dict to hold the results

    for d in distances:

        #Pull out all comparisons that have a given Levenshtein distance. Apply the occurance values for each kmer
        occurence_matrix = (result_matrix == d) * np.array(occurences)
        
        #sum the results
        occurence_sum = np.sum(occurence_matrix, axis=1)
       
        #add to results
        results[d] = occurence_sum

    return results


#function to compute the similarity in the predicted RNA folding of the designed bRNA to the reference bRNA using the RNAforester algorithm from Vienna RNA Suite
def rnaforester_score(seq1:str, seq2:str, db2:str):
    #seq*, db* are strings (sequence and dot-bracket notation of secondary structure).

    #predict the MFE folding of the designed bRNA
    seq1 = seq1
    db1, mfe = RNA.fold(seq1)

    # prepare FASTA-like input that RNAforester accepts
    payload = textwrap.dedent(f"""
    >x
    {seq1}
    {db1}
    >y
    {seq2}
    {db2}
    """).lstrip()

    # build CLI
    args = ["RNAforester"]
    args.append("-r")          #relative scoring from 0-1: sr(a,b) = 2*s(a,b)/(s(a,a)+s(b,b))
    args.append("--score")   # print only the optimal score

    # run
    proc = subprocess.run(args, input=payload.encode(), capture_output=True, check=True)
    out = proc.stdout.decode().strip()

    # when --score is used, stdout is usually just a number
    try:
        score = float(out.splitlines()[-1].split()[0])
    except Exception:
        score = None

    return out, score


#function to iterate through possible targets in the locus and design bridges
def iterate_bridge_design(target_locus, locus_name, **kwargs):

    #unpack the variables,
    genbank_files = kwargs.get("genbank_files","") #list of genbank files, against which to screen off targets
    donor_seq = kwargs.get("donor_seq","ACAGTATCTTGTAT") #default is the native donor of IS621
    donor_name = kwargs.get("donor_name","1") #name of donor 
    cores = kwargs.get("cores",['CT']) #cores to screen
    scaffold = kwargs.get("scaffold","IS621_WT") #bRNA scaffold to create. Can be IS621_WT, IS621_enhanced, ISCro4_WT, ISCro4_enhanced
    include_left = kwargs.get("include_left",7) #bases to include left of the core in the target
    include_right = kwargs.get("include_right",5) #bases to include right of the core in the target (include_left+include_right must equal 14)
    kmer_length = kwargs.get("kmer_length",11) #only consider the first Xbp for off target analysis. For IS621, 11 is relevant
    primer_seqs = kwargs.get("primer_seqs",{"CT":["",""],"GT":["",""],"AT":["",""],"TT":["",""]}) #if making a bRNA library this allows you to embed amplification primers for subpools 
    avoid_restriction = kwargs.get("avoid_restriction",["GGTCTC"]) #avoid these restriction sites 
    feature_type = kwargs.get("feature_type","") #define what type of feature this locus is
    check_offtargets = kwargs.get("check_offtargets", True) #check for off targets across sequences, provided as genbank files
    score_structure = kwargs.get("score_structure", True) #scores how well the predicted secondary structure od the designed bRNA matches the native
    distances = kwargs.get("distances",[0,1,2]) #Levenshtein distances to calculate for off targets 


    #open and read genbank files against which to score off targets if check_offtargets = True
    if check_offtargets == True:
        genome_seq = ''
        for genbank_file in genbank_files:
            for record in SeqIO.parse(genbank_file, "genbank"):
                genome_seq = genome_seq + str(record.seq).strip().upper() + 'AAAAAAAAAAAA' #in the future find some other way to separate the contigs

        genome_seq_rc = str(Seq(genome_seq).reverse_complement())

    #convert target locus to upper
    target_seq = target_locus.upper()
    donor_seq = donor_seq.upper()

    #start a dataframe to hold the contents
    target_df = pd.DataFrame()

    #find all potential target seqeunces 

    #create a counter for naming the bRNAs
    i=0
    for core in cores:

        #create an array to store target data
        target_array = []

        core = core.upper()
        #adjust the donor seq so that the cores match between donor and target
        donor_seq = donor_seq[0:7] + core + donor_seq[9:]

        core_rc = str(Seq(core).reverse_complement())

        #pull out primer seqs for eblock. In 5' to 3' orientation
        left_primer = primer_seqs[core][0]
        right_primer = primer_seqs[core][1]

        #get indexes for all potential + and - strand target sequences based on core sequence
        target_indexes = [[m,"+"] for m in find_all_loop(core, target_seq)]
        if core != core_rc:
            target_indexes.extend([[m,"-"] for m in find_all_loop(core_rc, target_seq)])

        #iterate through all potential target sequences
        for target_index_array in target_indexes:
            target_index = target_index_array[0] #pull out the index for the target

            #make sure index is not out too close to the edges
            if target_index > 7 and target_index < len(target_seq) - 10:
                i = i+1

                #create dict to hold data for this target site
                this_target_data = {}
                this_target_data['target_gene'] = locus_name
                if feature_type != "":
                    this_target_data['feature_type'] = feature_type
                this_target_data['donor_seq'] = donor_seq

                #create name for bridge
                bridge_name = f"bridge_IS621_T_{locus_name}_{i}_D_{donor_name}"
                this_target_data['bridge_name'] = bridge_name

                #check if this target is on the + or - strand and then pull out the target sequence 
                if target_index_array[1] == '+':
                    this_target_seq = target_seq[target_index-include_left:target_index] + core + target_seq[target_index + len(core): target_index + len(core) + include_right]
                    this_target_data['target_seq'] = this_target_seq #sequence of target
                    this_target_data['index'] = target_index + 1 #index is middle of core
                    this_target_data['strand'] = '+' #strand of target
                else:
                    this_target_seq = str(Seq(target_seq[target_index-include_right:target_index] + core_rc + target_seq[target_index + len(core): target_index + len(core) + include_left]).reverse_complement())
                    this_target_data['target_seq'] = this_target_seq #sequence of target
                    this_target_data['index'] = target_index + 1 #index is middle of core 
                    this_target_data['strand'] = '-' #strand of target

                this_target_data['core'] = core

                #design the bridge
                bridge_design, p6p7_warning, TBL_seq, DBL_seq = design_bridge_rna(this_target_seq, donor_seq, scaffold)
                this_target_data['full_bRNA_seq'] = bridge_design.bridge_sequence
                this_target_data['split_TBL_seq'] = TBL_seq
                this_target_data['split_DBL_seq'] = DBL_seq
                this_target_data['p6p7_warning'] = p6p7_warning


                #calculate RNA structural similarity of the designed bRNA to the reference Is621 bRNA
                if score_structure == True:
                    out, score = rnaforester_score(seq1 = bridge_design.bridge_sequence, seq2 = bridge_design.NATIVESQ, db2 = bridge_design.DOT_STRC)
                    this_target_data['RNA_structural_similarity'] = score

                #design the eblock
                eblock = eblock_design(bridge_design.bridge_sequence, left_primer=left_primer, right_primer=right_primer)
                this_target_data['cloning_fragment_seq_Patel_et_al'] = eblock

                #find 14mer perfect matches in the genome
                if check_offtargets == True:
                    perfect_matches = len(find_all_loop(this_target_seq, genome_seq)) + len(find_all_loop(this_target_seq, genome_seq_rc))
                    this_target_data[f'fulllength_matches_LevDist_0'] = perfect_matches

                #check to make sure there aren't disallowed restriction sites in the eBlock - if there are, then don't include this target
                total_restriction_sites = check_restriction(avoid_restriction, eblock)
                
                if total_restriction_sites == 0:
                    #append to array
                    target_array.append(this_target_data)

        #convert to a dataframe
        this_target_df = pd.DataFrame(target_array)

        #calculate off targets
        if check_offtargets == True:
            #pull all target sequences identified
            targets = this_target_df['target_seq'].tolist()
            #calculate all core-containing sequences of length kmer in the genome
            kmers = kmer_counter(genome_seq = genome_seq,kmer_length = kmer_length,core = core)
            off_targets = off_target_analysis(targets=targets, genome_seq=genome_seq, core=core, kmer_length=kmer_length, kmers = kmers,distances = distances)
            
            #add to dataframe
            for d in distances:
                this_target_df[f'{kmer_length}mer_matches_LevDist_{d}'] = off_targets[d]

        #append to existing dataframe
        target_df = pd.concat([target_df, this_target_df])
                
    #export a csv
    target_df.to_csv(f"{locus_name}.csv")

    #return a dataframe of the designed bRNAs
    return target_df

# function for designing a single bRNA
def single_bridge_design(target_seq, **kwargs):

    #unpack the variables,
    genbank_files = kwargs.get("genbank_files","") #list of genbank files, against which to screen off targets
    donor_seq = kwargs.get("donor_seq","ACAGTATCTTGTAT") #default is the native donor of IS621
    scaffold = kwargs.get("scaffold","IS621_WT") #bRNA scaffold to create. Can be IS621_WT, IS621_enhanced, ISCro4_WT, ISCro4_enhanced
    kmer_length = kwargs.get("kmer_length",11) #only consider the first Xbp for off target analysis. For IS621, 11 is relevant
    avoid_restriction = kwargs.get("avoid_restriction",[]) #avoid these restriction sites 
    check_offtargets = kwargs.get("check_offtargets", False) #check for off targets across sequences, provided as genbank files
    score_structure = kwargs.get("score_structure", True) #scores how well the predicted secondary structure od the designed bRNA matches the native
    distances = kwargs.get("distances",[0,1,2]) #Levenshtein distances to calculate for off targets 

    #identify core
    core = target_seq[7:9]

    #open and read genbank files against which to score off targets if check_offtargets = True
    if check_offtargets == True:
        genome_seq = ''
        for genbank_file in genbank_files:
            for record in SeqIO.parse(genbank_file, "genbank"):
                genome_seq = genome_seq + str(record.seq).strip().upper() + 'AAAAAAAAAAAA' #in the future find some other way to separate the contigs

        genome_seq_rc = str(Seq(genome_seq).reverse_complement())

    target_seq = target_seq.upper()
    donor_seq = donor_seq.upper()

    #start a dict to hold the data
    this_target_data = {}
    this_target_data['target_seq'] = target_seq #sequence of target
    this_target_data['donor_seq'] = donor_seq #sequence of donor
    this_target_data['core'] = core #core of target

    #design the bridge
    bridge_design, p6p7_warning, TBL_seq, DBL_seq = design_bridge_rna(target_seq, donor_seq, scaffold)
    this_target_data['full_bRNA_seq'] = bridge_design.bridge_sequence
    this_target_data['split_TBL_seq'] = TBL_seq
    this_target_data['split_DBL_seq'] = DBL_seq
    this_target_data['p6p7_warning'] = p6p7_warning

    #calculate RNA structural similarity of the designed bRNA to the reference Is621 bRNA
    if score_structure == True:
        out, score = rnaforester_score(seq1 = bridge_design.bridge_sequence, seq2 = bridge_design.NATIVESQ, db2 = bridge_design.DOT_STRC)
        this_target_data['RNA_structural_similarity'] = score

    #design the eblock
    eblock = eblock_design(bridge_design.bridge_sequence, left_primer="", right_primer="")
    this_target_data['cloning_fragment_seq_Patel_et_al'] = eblock

    #check to make sure there aren't disallowed restriction sites in the eBlock - if there are, then don't include this target
    total_restriction_sites = check_restriction(avoid_restriction, eblock)
    if total_restriction_sites > 0:
        this_target_data['restriction_site_warnings'] = f'{total_restriction_sites} disallowed restriction sites present'

    #find offtargets in the provided genbank files
    if check_offtargets == True:
        perfect_matches = len(find_all_loop(target_seq, genome_seq)) + len(find_all_loop(target_seq, genome_seq_rc))
        this_target_data[f'fulllength_matches_LevDist_0'] = perfect_matches

        targets = [target_seq]
        #calculate all core-containing sequences of length kmer in the genome
        kmers = kmer_counter(genome_seq = genome_seq,kmer_length = kmer_length,core = core)
        off_targets = off_target_analysis(targets=targets, genome_seq=genome_seq, core=core, kmer_length=kmer_length, kmers = kmers, distances = distances)
        
        #add to dataframe
        for d in distances:
            this_target_data[f'{kmer_length}mer_matches_LevDist_{d}'] = off_targets[d][0]

    #return a dict of the designed bRNA    
    return this_target_data

#function to assess off targets of a given 14mer
def off_target_assess(target_seq, genbank_files, **kwargs):

    #unpack the variables,
    kmer_length = kwargs.get("kmer_length",11) #only consider the first Xbp for off target analysis. For IS621, 11 is relevant
    distances = kwargs.get("distances",[0,1,2]) #Levenshtein distances to calculate for off targets    

    #open and read genbank files against which to score off targets if check_offtargets = True
    genome_seq = ''
    for genbank_file in genbank_files:
        for record in SeqIO.parse(genbank_file, "genbank"):
            genome_seq = genome_seq + str(record.seq).strip().upper() + 'AAAAAAAAAAAA' #in the future find some other way to separate the contigs

    genome_seq_rc = str(Seq(genome_seq).reverse_complement())

    target_seq = target_seq.upper()

    #determine core sequence
    core = target_seq[7:9]

    #dict to hold the data
    this_target_data={}
    this_target_data['target_seq'] = target_seq
    this_target_data['core'] = core #core of target

    #calculate the offtargets
    perfect_matches = len(find_all_loop(target_seq, genome_seq)) + len(find_all_loop(target_seq, genome_seq_rc))
    this_target_data[f'fulllength_matches_LevDist_0'] = perfect_matches

    targets = [target_seq]
    #calculate all core-containing sequences of length kmer in the genome
    kmers = kmer_counter(genome_seq = genome_seq,kmer_length = kmer_length,core = core)
    off_targets = off_target_analysis(targets=targets, genome_seq=genome_seq, core=core, kmer_length=kmer_length, kmers=kmers, distances = distances)
    
    #add to dataframe
    for d in distances:
        this_target_data[f'{kmer_length}mer_matches_LevDist_{d}'] = off_targets[d][0]
    
    #return a dict of the offtarget information
    return this_target_data

# function for finding unique 14mers in a genome
def find_unique_targets(genbank_files, **kwargs):

    #unpack the variables,
    cores = kwargs.get("cores",['CT']) #cores to screen
    kmer_length = kwargs.get("kmer_length",11) #only consider the first Xbp for off target analysis. For IS621, 11 is relevant
    distances = kwargs.get("distances",[0,1,2]) #Levenshtein distances to calculate for off targets
    cutoff = kwargs.get("cutoff",10) #Number of unique sequences to find per core 

    #open and read genbank files against which to score off targets if check_offtargets = True
    genome_seq = ''
    for genbank_file in genbank_files:
        for record in SeqIO.parse(genbank_file, "genbank"):
            genome_seq = genome_seq + str(record.seq).strip().upper() + 'AAAAAAAAAAAA' #in the future find some other way to separate the contigs

    target_df = pd.DataFrame()

    for core in cores:
        #first find the most unique 14mers with core sequence. Pull out 10x the cutoff so we can then sort by mismatches 
        kmers_14mers = kmer_counter(genome_seq = genome_seq,kmer_length = 14,core = core)
        least_common_14mers = heapq.nsmallest(cutoff*10, kmers_14mers.items(), key=lambda x: x[1])
        del kmers_14mers

        this_target_df = pd.DataFrame(least_common_14mers, columns=['target_seq', 'fulllength_matches_LevDist_0'])

        targets= this_target_df["target_seq"].tolist()

        #then calculate off targets of the shorter kmer
        kmers = kmer_counter(genome_seq = genome_seq,kmer_length = kmer_length,core = core)
        off_targets = off_target_analysis(targets=targets, genome_seq=genome_seq, core=core, kmer_length=kmer_length, kmers=kmers, distances = distances)

        #add to dataframe
        for d in distances:
            this_target_df[f'{kmer_length}mer_matches_LevDist_{d}'] = off_targets[d]

        #append to existing dataframe
        target_df = pd.concat([target_df, this_target_df])

    if target_df.empty:
        return target_df

    #sort the list by columns
    columns = ['fulllength_matches_LevDist_0'] + [f'{kmer_length}mer_matches_LevDist_{d}' for d in distances]
    target_df_sorted = target_df.sort_values(by=columns, ignore_index=True)

    #return a dataframe of the identified sites
    return target_df_sorted.head(cutoff)

#for testing 

#genbank_files = ["/Users/jayman1466/My Drive (jayman1466@gmail.com)/Biome/New Project Organization/IS110/Code/Insertion Mapper/genomes/MG1655 split rRNA.gb"]
#print(single_bridge_design("TTTCACCCTGGAGG", genbank_files = genbank_files, check_offtargets = True))
#iterate_bridge_design("GTTTTCACCCTGGAGGATTTTGTGGGTGATTGGCGGCAGACCGCCGGTTACAACCTGGACCAGGTACTGGAACAGGGAGGTGTCTCCTCCCTGTTCCAGAACTTGGGGGTGTCGGTGACGCCGATACAGCGCATCGTGCTGTCGGGTGAAAACGGCCTGAAGATTGACATTCATGTGATCATTCCATATGAAGGCCTATCTGGCGACCAGATGGGCCAGATTGAGAAGATCTTCAAGGTCGTTTACCCGGTGGATGATCACCATTTTAAGGTGATCCTCCACTATGGTACCCTGGTCATTGACGGGGTGACCCCGAATATGATCGATTACTTCGGGCGACCCTATGAAGGGATCGCAGTATTCGATGGGAAGAAGATAACCGTTACCGGAACGCTCTGGAATGGTAACAAAATCATTGATGAGCGGCTGATCAACCCGGACGGATCACTGCTGTTCCGGGTGACCATTAACGGGGTCACTGGATGGCGGCTATGTGAACGCATCCTTGCCTAA", "nanoluc", genbank_files = ["/Users/jayman1466/My Drive (jayman1466@gmail.com)/Biome/New Project Organization/IS110/Code/Insertion Mapper/genomes/MG1655 split rRNA.gb"])
#print(off_target_assess("TTTCACCCTGGAGG", genbank_files))
#print(find_unique_targets(genbank_files))