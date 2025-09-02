from Bio.Seq import Seq
from Bio import SeqIO
import Levenshtein
import re
import pandas as pd
import subprocess, tempfile, textwrap, shlex, pathlib
import RNA
from .classes import WTBridgeRNA177nt

#author: Jaymin Patel jayman1466@gmail.com
VERSION = "1.0.0"
__all__ = ["design_bridges"]

#function to design the IS621 bridgeRNA using Matt's script https://github.com/hsulab-arc/BridgeRNADesigner
def design_bridge_rna(target, donor):

    target = target.upper()
    donor = donor.upper()

    WTBridgeRNA177nt.check_target_length(target)
    WTBridgeRNA177nt.check_donor_length(donor)
    WTBridgeRNA177nt.check_core_match(target, donor)
    WTBridgeRNA177nt.check_target_is_dna(target)
    WTBridgeRNA177nt.check_donor_is_dna(donor)

    brna = WTBridgeRNA177nt()
    brna.update_target(target)
    brna.update_donor(donor)
    brna.update_hsg()

    return brna


#function to design an eblock for cloning via golden gate in bsaI sites
def eblock_design(input_seq,left_primer="",right_primer=""):
    #print(input_seq)
    left_remove = 'AGTGCAGAGAAAATCGGCCAGTTTTCTCTGCCTGCAGTCCGCATGCCGT'
    right_remove = 'TGGTTTCACT'
    left_gg = 'GAGAGggtctcTCCGT' #left bsaI site
    right_gg = 'TGGTAgagaccGAGAG' #right bsaI site
    stuffer = "GACATTGTCCCTGATTTCTCCACTACTAATAGCACACACGGGGCAATACCAGCACAAGCTAGTCTCGCGGGAACGCTCGTCAGCATACGAAAGAGCTTAAGGCACGCCAATTCGCACTGTCAGGGTCACTTGGGTGTTTTGCACTACCGT" #150bp stuffer to get to 300bp fragment size

    
    #append additional PCR primers to ends for amplification of eBlock. This expects both primers to be in 5' to 3' orientation 
    left_OH = left_primer + left_gg
    right_OH = right_gg + str(Seq(right_primer).reverse_complement()) + stuffer

    mod_seq = input_seq.replace(left_remove,left_OH)
    mod_seq = mod_seq.replace(right_remove,right_OH)
    return mod_seq


#function to quantify number of perfect match targets in the genome. This only evaluates the first N bp defined by the kmer argument (default: 11bp)
def perfect_match_offtargets(this_target_seq,genome_seq,kmer):
    offtarget_query = this_target_seq[0:kmer]            
    offtarget_query_rc = str(Seq(offtarget_query).reverse_complement())

    offtarget_count = genome_seq.count(offtarget_query)
    offtarget_rc_count = genome_seq.count(offtarget_query_rc)

    offtargets = offtarget_count + offtarget_rc_count

    return offtargets


#function to quantify number of imperfect off targets in the genome (but still retain a perfect match at the core). This is based on levenshtein distance. With this implementation, indels are given a distance score of 2, and substitutions a score of 1. Scores <= 2 are reported. This only evaluates the first N bp defined by the kmer argument (default: 11bp).
def imperfect_match_offtargets(this_target_seq, genome_seq, kmer, core, core_rc, include_left, include_right):
    
    lev_1_count=0 #running total of offtargets with Levenshtein distance of 1
    lev_2_count=0 #running total of offtargets with Levenshtein distance of 2
    
    #identify all index values of potential off targets based on core matches
    offtarget_indexes = [[m.start(),"+"] for m in re.finditer(core, genome_seq)]
    if core != core_rc:
        offtarget_indexes.extend([[m.start(),"-"] for m in re.finditer(core_rc, genome_seq)])
    
    #iterate through all potential off target sequences
    for offtarget_index_array in offtarget_indexes:
        offtarget_index = offtarget_index_array[0] #pull out the index for the target

        #make sure index is not out too close to the edges
        if offtarget_index > 7 and offtarget_index < len(genome_seq) - 10:

            #check if this potential off target is on the + or - strand and then pull out the target sequence 
            if offtarget_index_array[1] == '+':
                this_offtarget_seq = genome_seq[offtarget_index-include_left:offtarget_index] + core + genome_seq[offtarget_index + len(core): offtarget_index + len(core) + include_right]
            else:
                this_offtarget_seq = str(Seq(genome_seq[offtarget_index-include_right:offtarget_index] + core_rc + genome_seq[offtarget_index + len(core): offtarget_index + len(core) + include_left]).reverse_complement())
            
            levenshtein_distance = Levenshtein.distance(this_target_seq[0:kmer], this_offtarget_seq[0:kmer], score_cutoff=2)
            
            #add to running totals
            if levenshtein_distance == 1:
                lev_1_count += 1
            if levenshtein_distance == 2:
                lev_2_count += 1
    
    return [lev_1_count,lev_2_count]

#function to compute the similarity in the predicted RNA folding of the designed bRNA to the reference bRNA using the RNAforester algorithm from Vienna RNA Suite
def rnaforester_score(seq1):
    #seq*, db* are strings (sequence and dot-bracket).

    #reference folding of the native IS621 bRNA
    seq2 = "AGTGCAGAGAAAATCGGCCAGTTTTCTCTGCCTGCAGTCCGCATGCCGTATCGGGCCTTGGGTTCTAACCTGTTGCGTAGATTTATGCAGCGGACTGCCTTTCTCCCAAAGTGATAAACCGGACAGTATCATGGACCGGTTTTCCCGGTAATCCGTATTTACAAGGCTGGTTTCACT"
    db2 =  "...(((((((((((......))))))))))).((((((((((.(((............(((((....)))))..............))))).)))))))).........((((....(((.............((((((.....)))))...)...............)))..))))"

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
def iterate_bridge_design(target_locus, target_name, genbank_file, **kwargs):

    #unpack the variables,
    donor_seq = kwargs.get("donor_seq","ACAGTATCTTGTAT") #default is the native donor of IS621
    donor_name = kwargs.get("donor_name","1")
    cores = kwargs.get("cores",['CT'])
    include_left = kwargs.get("include_left",7)
    include_right = kwargs.get("include_right",5)
    kmer = kwargs.get("kmer",11)
    primer_seqs = kwargs.get("primer_seqs",{"CT":["",""],"GT":["",""],"AT":["",""],"TT":["",""]})
    avoid_restriction = kwargs.get("avoid_restriction",[])
    feature_type = kwargs.get("feature_type","")
    check_imperfect = kwargs.get("feature_type", True)
    score_structure = kwargs.get("score_structure", True)

    #open and read genbank file for the reference genome
    genome_seq = ''
    for record in SeqIO.parse(genbank_file, "genbank"):
        genome_seq = genome_seq + record.seq + 'AAAAAAAAAAAA' #in the future find some other way to separate the contigs

    #convert everything to upper
    target_seq = target_locus.upper()
    genome_seq = str(genome_seq.upper())

    #create an array to store target data
    target_array = []

    #find all potential target seqeunces 
    i=0
    for core in cores:

        #adjust the donor seq so that the cores match between donor and target
        donor_seq = donor_seq[0:7] + core + donor_seq[9:]

        core_rc = str(Seq(core).reverse_complement())

        #pull out primer seqs for eblock. In 5' to 3' orientation
        left_primer = primer_seqs[core][0]
        right_primer = primer_seqs[core][1]

        #get indexes for all potential + and - strand target sequences based on core sequence
        target_indexes = [[m.start(),"+"] for m in re.finditer(core, target_seq)]
        if core != core_rc:
            target_indexes.extend([[m.start(),"-"] for m in re.finditer(core_rc, target_seq)])

        #iterate through all potential target sequences
        for target_index_array in target_indexes:
            target_index = target_index_array[0] #pull out the index for the target

            #make sure index is not out too close to the edges
            if target_index > 7 and target_index < len(target_seq) - 10:
                i = i+1

                #create dict to hold data for this target site
                this_target_data = {}
                this_target_data['target_gene'] = target_name
                this_target_data['feature_type'] = feature_type
                this_target_data['donor_seq'] = donor_seq

                #create name for bridge
                bridge_name = f"bridge_IS621_T_{target_name}_{i}_D_{donor_name}"
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


                #Quantify number of perfect targets in the genome.
                this_target_data['perfect_match_targets'] = perfect_match_offtargets(this_target_seq,genome_seq,kmer)

                #Quantify number of imperfect off targets in the genome.
                if check_imperfect == True:
                    imperfect_matches = imperfect_match_offtargets(this_target_seq, genome_seq, kmer, core, core_rc, include_left, include_right)
                    this_target_data['levenshtein_distance_1_targets'] = imperfect_matches[0]
                    this_target_data['levenshtein_distance_2_targets'] = imperfect_matches[1]

                #TO DO: WOBBLE MISMATCHES

                #design the bridge
                bridge_design = design_bridge_rna(this_target_seq, donor_seq)
                this_target_data['bridge_sequence'] = bridge_design.bridge_sequence
                if bridge_design.p6p7_match == True:
                    this_target_data['p6p7_warning'] = "Target P6-P7 and Donor P6-P7 match, efficiency is unclear. DBL HSG forced to be anti-complementary."
                else:
                    this_target_data['p6p7_warning'] = ""

                #calculate RNA structural similarity of the designed bRNA to the reference Is621 bRNA
                if score_structure == True:
                    out, score = rnaforester_score(bridge_design.bridge_sequence)
                    this_target_data['RNA_structural_similarity'] = score

                #design the eblock
                eblock = eblock_design(bridge_design.bridge_sequence, left_primer=left_primer, right_primer=right_primer)
                this_target_data['eblock_seq'] = eblock

                #check to make sure there aren't disallowed restriction sites in the eBlock
                total_restriction_sites = 0
                for restriction_site in avoid_restriction:
                    total_restriction_sites = total_restriction_sites + eblock.count(restriction_site)
                
                if total_restriction_sites == 0:
                    #append to array
                    target_array.append(this_target_data)
                
    #return the array
    return target_array

def design_bridges(target_locus, target_name, genbank_file, **kwargs):

    #run the pipeline
    target_array = iterate_bridge_design(target_locus=target_locus, target_name=target_name, genbank_file=genbank_file, **kwargs)

    #convert to pandas dataframe
    target_df = pd.DataFrame(target_array)

    #export
    target_df.to_csv(f"{target_name}.csv")
