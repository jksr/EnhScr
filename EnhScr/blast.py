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



def run_mmseqs2(qseqs, qids, tseqs, tids, tmpdir=None):
    assert is_tool('mmseqs'), 'no mmseqs found! recomend "mmseqs2\t7.4e23d\th21aa3a5_1    bioconda/label/cf201901"'
    
    qfasta_seqs = ( SeqRecord(Seq(qfa),qid) for qid,qfa in zip(qids,qfas) ) 
    tfasta_seqs = ( SeqRecord(Seq(tfa),tid) for tid,tfa in zip(tids,tfas) ) 
    
    if tmpdir is not None:
        tmpdir = Path(tmpdir)
        tmpdir.mkdir(exist_ok=True, parents=True)
    else:
        import tempfile
        tmpdir = Path(next(tempfile._get_candidate_names()))
        tmpdir.mkdir()
        
        
    qfn = f'{tmpdir}/tmpquery.fasta'
    tfn = f'{tmpdir}/tmptarget.fasta'
    SeqIO.write(qfasta_seqs, qfn, "fasta")
    SeqIO.write(tfasta_seqs, tfn, "fasta")

    mmseqs_tmp = f'{tmpdir}/tmpmmseqs'
    ofn = f'{tmpdir}/tmpout.tsv'
    
    assert not ofn.exists(), f'there is previous results in {tmpdir}, please specify another temp folder or omit it'
    
    mmseqs_format = '"query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits"'
    cmd = f'mmseqs easy-search {qfn} {tfn} {ofn} {mmseqs_tmp} --format-output {mmseqs_format}'
    cmd = f'mmseqs easy-search {qfn} {tfn} {ofn} {mmseqs_tmp}'
#     print(cmd)


    pr = subprocess.Popen(shlex.split(cmd), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    # blocks until process finishes
    out, err = pr.communicate() 
    # check the return code
    if pr.returncode != 0:
        pass

    outdf = pd.read_csv(ofn, sep='\t', names = ['query','target','pident','alnlen','mismatch','gapopen',
                                                'qstart','qend','tstart','tend','evalue','bits'])
#     print(ofn)
    return outdf
