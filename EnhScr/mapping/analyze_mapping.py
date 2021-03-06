#import pandas as pd
#import seaborn as sns
#import matplotlib.pyplot as plt
import pysam
import numpy as np
from Bio import SeqIO
import gzip
from collections import Counter
import joblib


class AlignmentRead:
    def __init__(self, rec):
        if not isinstance(rec, str):
            raise TypeError
        rec = rec.split('\t')
        
        self._rec = rec
        self.readid = rec[0]
        self.flag = "{0:b}".format(int(rec[1])).zfill(12)
        self.refid = rec[2]
        self.start = rec[3]
        self.mapq = rec[4]
        self.cigar = rec[5]
        
        
        self.seq = rec[9]
        
    def isR1(self):
        return self.flag[5]=='1'
    def isR2(self):
        return self.flag[4]=='1'

    
class AlignmentPair:
    def __init__(self, pair):
        if isinstance(pair[0], AlignmentRead):
            read1 = pair[0]
        elif isinstance(pair[0], str):
            read1 = AlignmentRead(pair[0])
        if isinstance(pair[1], AlignmentRead):
            read2 = pair[1]
        elif isinstance(pair[1], str):
            read2 = AlignmentRead(pair[1])
        
        assert read1.isR1()!=read2.isR1() and read1.isR2()!=read2.isR2(), f"{read1.readid}\n{read2.readid}"
        assert read1.refid ==read2.refid
        
        if read1.isR2():
            read1, read2 = read2, read1
        
        self.R1 = read1
        self.R2 = read2
        
class Alignments:
    def __init__(self, alignment_pair_list):
        self.alignments = alignment_pair_list
        self.indexing = dict()
        
        for i,ap in enumerate(self.alignments):
            self.indexing[ap.R1.readid] = i
                
    def __getitem__(self, key):
        if isinstance(key, str):
            return self.alignments[self.indexing[key]]
        elif isinstance(key, int):
            return self.alignments[key]
        else:
            raise TypeError
    def __len__(self):
        return len(self.alignments)
        
        
def parse_sam(fn, sam_filter_args=('-q', '30', '-f', '0x2')):
    recs = pysam.view(*sam_filter_args, fn).strip('\n')
    recs = recs.split('\n')
    recs = [AlignmentRead(rec) for rec in recs]
    
    tmp = dict()
    for rec in recs:
        if rec.readid not in tmp:
            tmp[rec.readid] = []
        tmp[rec.readid].append(rec)
    recs = tmp
    
    alignment_pairs = []
    for pair in recs.values():
        if len(pair)!=2:
            continue
        alignment_pairs.append(AlignmentPair(pair))
#     print(len(recs))
#     print(recs[0])
#     print(recs[-1])
#     assert len(recs)%2==0
    
#     alignment_pairs = [ AlignmentPair(recs[i:i+2]) for i in range(0,len(recs),2) ]
    
    return Alignments(alignment_pairs)
        

def match_reads_and_refs(alignment_pairs):
    read_to_ref = dict()
    ref_to_read = dict()

    for ap in alignment_pairs:
        refid, readid  = ap.R2.refid, ap.R2.readid
        if refid not in ref_to_read:
            ref_to_read[refid] = []
        ref_to_read[refid].append(readid)
        read_to_ref[readid] = refid
    return read_to_ref, ref_to_read

def match_reads_and_barcodes(r2fn, barcode_length, use_reads=None, barcode_trim_length=0):
    read_to_barcode = dict()
    barcode_to_read = dict()
    with gzip.open(r2fn, "rt") as f:
        for record in SeqIO.parse(f, "fastq"):
            readid = record.id
            barcode = str(record.seq[barcode_trim_length:barcode_length])
            
            if use_reads is not None:
                if readid not in use_reads:
                    continue
            
            read_to_barcode[readid] = barcode
            if barcode not in barcode_to_read: 
                barcode_to_read[barcode] = []
            barcode_to_read[barcode].append(readid)
    return read_to_barcode, barcode_to_read
            
def match_refs_and_barcodes(reads_to_refs, refs_to_reads, reads_to_barcodes, barcodes_to_reads):
    barcodes_to_refs = dict()
    for bc in barcodes_to_reads:
        refs = list( reads_to_refs[rd] for rd in barcodes_to_reads[bc] )
        barcodes_to_refs[bc] = refs
        
    refs_to_barcodes = dict()
    for ref in refs_to_reads:
        barcodes = list( reads_to_barcodes[rd] for rd in refs_to_reads[ref] )
        refs_to_barcodes[ref] = barcodes
       
    return barcodes_to_refs, refs_to_barcodes
        
            
        
        

class AlignmentResults:
    def __init__(self, samfn, r2fn, barcode_length=34, barcode_trim_length=0, 
                 sam_filter_args='-q 30 -f 0x2'):
        sam_filter_args = sam_filter_args.split()
        self.alignments = parse_sam(samfn,sam_filter_args=sam_filter_args)
        reads_to_refs, refs_to_reads = match_reads_and_refs(self.alignments)
        reads_to_barcodes, barcodes_to_reads = match_reads_and_barcodes(r2fn, barcode_length, 
                                                                        set(reads_to_refs.keys()), 
                                                                        barcode_trim_length)
        barcodes_to_refs, refs_to_barcodes = match_refs_and_barcodes(reads_to_refs, refs_to_reads, 
                                                                     reads_to_barcodes, barcodes_to_reads)
        
        self.reads_to_refs = reads_to_refs 
        self.refs_to_reads = refs_to_reads
        self.reads_to_barcodes = reads_to_barcodes
        self.barcodes_to_reads = barcodes_to_reads
        self.barcodes_to_refs = barcodes_to_refs
        self.refs_to_barcodes = refs_to_barcodes
        
    #def TODO    


def analyze_mapping(alignment_path, R2_path, output_prefix, barcode_length, barcode_trim_length, sam_filter_args):
    alignment_results = AlignmentResults(alignment_path, R2_path, barcode_length, barcode_trim_length, sam_filter_args)
    joblib.dump(alignment_results, f'{output_prefix}.joblib.gz')
