#!/usr/bin/python
# Author: SJ Riesenfeld
# Feb 24, 2018
# Cite-seq analysis

from Bio import SeqIO
import sys
import argparse

## Updated 2018 Jan 23: Should work now for arbitrary barcode length. Tested at least for 6-bp barcodes with max_dist=2 and 15-bp barcodes with max_dist==2.
### Potential corrections/updates needed:
## Now that the BC is 15 bp, it might make sense to really implement
## or borrow a true edit distance, which should take into account
## insertions and deletions. This does not currently do that; it just
## allows a sliding window.

# for rec in BC_ref.keys():
#    print (BC_ref[rec])

def simpledist (a, b):
    return sum(map(lambda (x, y): 0 if x == y else 1, zip(a, b)))

def distance_at_most(k, a, b):
    return simpledist (a, b) <= k
    
def do_align(reads_fastq, out_sam, max_dist, start, num, ref_fasta):
    BC_ref=SeqIO.to_dict(SeqIO.parse(ref_fasta, "fasta"))    
    keysv=sorted(BC_ref.keys())
    n_recs=0
    n_proc=0
    n_aln=0
    if (len(set([len(BC_ref[x].seq) for x in keysv]))>1):
        print("Using length of barcode "+keysv[1]+" for default barcode length")
    BC_len=len(BC_ref[keysv[1]].seq)
    print("Antibody barcode length is " + str(BC_len) + " bp")
    with open(reads_fastq, "rU") as fin:
        with open(out_sam, "w") as fout:
            fout.write("@HD\tVN:1\n")
            for x in keysv:
                fout.write("@SQ\tSN:"+x+"\tLN:"+str(BC_len)+"\n")
            fout.write("@PG\tID:align_reads_BCs.py\tPN:align_reads_BCS.py\n")
            for record in SeqIO.parse(fin, "fastq"):
                n_recs=n_recs+1
                #print("n_recs: "+str(n_recs))
                if (n_recs < start):
                    continue
                if ((num!=-1) and (n_recs >= (start+num))):
                    break
                if (n_recs == start):
                    print("ID of first sequence processed: " + record.id)
                if ((num!=-1) and (n_recs == (start+num-1))):
                    print("ID of last sequence processed: " + record.id)
                n_proc=n_proc+1
                checkhits=[record.seq.find(BC_ref[x].seq) for x in keysv]            
                hits=[i for i, e in enumerate(checkhits) if e != -1]
                cigarmatch=""
                if len(hits)>0:
                    cigarmatch=str(BC_len)+"M"                
                else:
                    #print("No exact match")
                    checkhits=[-1]*len(keysv)
                    for j,x in enumerate(keysv):
                        # print str(j) + ",  "+x
                        for i in range(len(record.seq)-BC_len+1):
                            subst=str(record.seq[i:(i+BC_len)])
                            #print(subst)
                            #print(BC_ref[x].seq)
                            d_tf=distance_at_most(max_dist, subst, BC_ref[x].seq)
                            #print(d_tf)
                            if d_tf==True:
                                checkhits[j]=i
                                break
                    hits=[i for i, e in enumerate(checkhits) if e != -1]
                if len(hits)==0: continue
                n_aln=n_aln+1
                p=checkhits[hits[0]]
                match_key=keysv[hits[0]]
                matched=record.seq[p:(p+BC_len)]
                if (cigarmatch==""):
                    mismatchind=[i for i in range(BC_len) if matched[i]!=BC_ref[match_key].seq[i]]
                    #mismatchind=[0,10] # test
                    if len(mismatchind)>max_dist:
                        print("Unexpected number of mismatches: "+len(mismatchind))
                        exit()
                    for i,j in enumerate(mismatchind):
                        #print str(i)+", "+str(j)
                        if (i==0): # initialize
                            x_run_start=j # set start index of current run of consecutive mismatches
                            j_prev=-1
                            n_leadm=j # n_leadm could be 0 if string starts with mismatches
                        elif ((j-j_prev)>1): # j_prev should be >0; nonconsecutive mismatches
                            if (n_leadm>0): 
                                cigarmatch=cigarmatch+str(n_leadm)+"M" # record number of leading matches
                                n_leadm=0 # done recording leading M's; reset to 0
                            cigarmatch=cigarmatch+str(j_prev-x_run_start+1)+"X" # record length of previous run of mismatches
                            cigarmatch=cigarmatch+str(j-j_prev-1)+"M" # record number of matches up to current mismatch                            
                            x_run_start=j # reset start index of current run of mismatches
                        if (i==(len(mismatchind)-1)): # on final mismatch, still need to record it
                            if (n_leadm>0): # in case leading M's were not yet recorded
                                cigarmatch=cigarmatch+str(n_leadm)+"M" # record number of leading matches
                            cigarmatch=cigarmatch+str(j-x_run_start+1)+"X" # record length of current run of mismatches
                            if (j<(BC_len-1)): # end with remaining matches
                                cigarmatch=cigarmatch+str(BC_len-j-1)+"M"
                        j_prev=j                        
                        #print(cigarmatch)
                s=len(record.seq)-BC_len-p # extra characters
                cigar=""
                if (p>0):
                    cigar=str(p) + "S"
                cigar=cigar+cigarmatch
                if (s>0):
                    cigar=cigar + str(s) + "S"
                # print(cigar)
                fout.write(record.id + "\t0\t"+ match_key + "\t1\t255\t" + cigar + "\t*\t0\t0\t" + str(record.seq) +"\t*\n")
                if (n_recs % 1e5)==0:                        
                    print(str(n_proc) + " sequences processed and " + str(n_aln) + " aligned")
    if (num==-1):
        print("ID of last sequence processed: " + record.id)
    print(str(n_proc) + " total sequences processed from seq index " + str(start) + " through " + str(start+n_proc-1) + "; " + str(n_aln) + " aligned")
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-ref_fasta', help='fasta file of reference bar codes', default="/ahg/regevdata/projects/CiteSeq/Reference/AB_BC.ADT.fa")
    parser.add_argument('-reads_fastq','-i', help='input reads fastq file',
                        default="fastq_ADT/unaligned_mc_tagged_polyA_filtered.fastq")
    parser.add_argument('-out_dir','-o', help='output directory for SAM files', default="./")
    parser.add_argument('-seq_num_start', '-s', type=int, help='1-based index of sequence to start processing', default=1)
    parser.add_argument('-do_num', '-n', type=int, help='number of sequences to align (-1 to process all)', default=-1)
    parser.add_argument('-max_dist', '-d', type=int, help='maximum allowed simple distance between read and barcode', default=1)
    args = parser.parse_args()
    out_sam=args.out_dir+"/sjrAligned.d_" + str(args.max_dist) + ".s_" + str(args.seq_num_start) + ".n_" + str(args.do_num) + ".sam"
    print("Aligning reads from " +args.reads_fastq + " to sequences in " + args.ref_fasta)
    print("SAM output will be written to file: " +out_sam)
    print("arguments parsed: ")
    print(args)
    do_align(args.reads_fastq, out_sam, args.max_dist, args.seq_num_start, args.do_num, args.ref_fasta)


## test from command line:
# python scripts/align_reads_to_BCs.py -s 10 -n 2
# ref_fasta="/ahg/regevdata/projects/CiteSeq/Reference/AB_BC.ADT.fa"
# reads_fastq="testdir/unaligned_mc_tagged_polyA_filtered.B16_S3_AB.head100.fastq"
# out_dir="testdir"
# seq_num_start=1
# do_num=25
# max_dist=2
# out_sam=out_dir+"/sjrAligned.d_" +  str(max_dist) + ".s_" + str(seq_num_start) + ".n_" + str(do_num) + ".sam"
