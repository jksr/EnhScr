from pathlib import Path
import pandas as pd
import tempfile
import subprocess
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import shlex
from shutil import which



def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""
    return which(name) is not None


def mmseqs2_from_fasta_files(qfn, tfn, outfn=None, tmpdir=None, ):
    assert is_tool('mmseqs'), 'no mmseqs found! recomend "mmseqs2\t7.4e23d\th21aa3a5_1    bioconda/label/cf201901"'
    
    if tmpdir is not None:
        tmpdir = Path(tmpdir)
        tmpdir.mkdir(exist_ok=True, parents=True)
    else:
        tmpdir = Path(tempfile.gettempdir())/next(tempfile._get_candidate_names())
        tmpdir.mkdir()


    mmseqs_tmp = f'{tmpdir}/tmpmmseqs'
    if outfn is None:
        ofn = f'{tmpdir}/tmpout.tsv'
    else:
        ofn = outfn
    
    assert not Path(ofn).exists(), f'output file {ofn} already exists!'
    
    mmseqs_format = '"query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"'
#     cmd = f'mmseqs easy-search {qfn} {tfn} {ofn} {mmseqs_tmp} --format-output {mmseqs_format}'
    cmd = f'mmseqs easy-search {qfn} {tfn} {ofn} {mmseqs_tmp}'

    pr = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # blocks until process finishes
    out, err = pr.communicate()
    
    # check the return code
    if pr.returncode != 0:
        raise RuntimeError('Something wrong in running mmseqs2!')

    outdf = pd.read_csv(ofn, sep='\t', names = ['query','target','pident','alnlen','mismatch','gapopen',
                                                'qstart','qend','tstart','tend','evalue','bits'])
    return outdf



def mmseqs2_from_lists(qfas, qids, tfas, tids, tmpdir=None):
    assert is_tool('mmseqs'), 'no mmseqs found! recomend "mmseqs2\t7.4e23d\th21aa3a5_1    bioconda/label/cf201901"'
    
    qfasta_seqs = ( SeqRecord(Seq(qfa),qid) for qid,qfa in zip(qids,qfas) ) 
    tfasta_seqs = ( SeqRecord(Seq(tfa),tid) for tid,tfa in zip(tids,tfas) ) 
    
    if tmpdir is not None:
        tmpdir = Path(tmpdir)
        tmpdir.mkdir(exist_ok=True, parents=True)
    else:
        tmpdir = Path(tempfile.gettempdir())/next(tempfile._get_candidate_names())
        tmpdir.mkdir()
        
        
    qfn = f'{tmpdir}/tmpquery.fasta'
    tfn = f'{tmpdir}/tmptarget.fasta'
    SeqIO.write(qfasta_seqs, qfn, "fasta")
    SeqIO.write(tfasta_seqs, tfn, "fasta")
    
    return mmseqs2_from_fasta_files(qfn, tfn, None, tmpdir)
