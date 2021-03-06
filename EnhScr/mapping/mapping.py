import subprocess
from pathlib import Path
import tempfile
import shlex
import shutil


def _r2_premapping(r2fn, n_barcode_bps=34, n_shared_bps=72, trim_barcode=True, trim_shared=False, tmpdirn=None):
    n_bp_trimming = 0
    if trim_barcode:
        n_bp_trimming += n_barcode_bps
    if trim_shared:
        n_bp_trimming += n_shared_bps

    if tmpdirn is not None:
        tmpdir = Path(tmpdirn)
        tmpdir.mkdir(exist_ok=True)
    else:
        tmpdir = Path( next(tempfile._get_candidate_names()) )
        tmpdir.mkdir()

    tmpfn = tmpdir/f'{next(tempfile._get_candidate_names())}.fastq.gz'
    outfn = tmpfn.with_suffix('').with_suffix('.cutadapt.txt')

    cmd = f'cutadapt -u {n_bp_trimming} -o {tmpfn} {r2fn}'

    with open(outfn,'w') as f:
        subprocess.run(shlex.split(cmd), stdout=f, stderr=f)

    return tmpfn, outfn

def _map(r1fn, r2fn, prefix, refdir, e2e=True, maxins=1000, minins=0, bowtie2_params=''):
    opts = ''
    if e2e:
        opts += ' --end-to-end'
        postfix = 'e2e'
    else:
        opts += ' --local'
        postfix = 'local'
    
    samfn = f'{prefix}.{postfix}.sam'
    reportfn = f'{prefix}.{postfix}.report.txt'
    opts += f' -X {maxins} -I {minins} -x {refdir} -1 {r1fn} -2 {r2fn} -S {samfn}'


    cmd = f'bowtie2 {opts} {bowtie2_params}'
    
    with open(reportfn,'w') as f:
        subprocess.run(shlex.split(cmd), stdout=f, stderr=f)
    
    return samfn, reportfn



def mapping(R1_path, R2_path, output_prefix, reference_folder, 
        max_insertion,
        min_insertion,
        barcode_length,
        shared_length,
        trim_shared,
        local_mapping,
        bowtie2_params,):


    prefixdir = Path(output_prefix).parent
    prefixdir.mkdir(exist_ok=True)
    tmpdirn = next(tempfile._get_candidate_names())

    tmpr2fn,_ = _r2_premapping(R2_path, n_barcode_bps=barcode_length, n_shared_bps=barcode_length, trim_barcode=True, trim_shared=trim_shared, tmpdirn=tmpdirn)

    samfn, reportfn = _map(R1_path, tmpr2fn, output_prefix, reference_folder, e2e=not local_mapping,
                           maxins=max_insertion, minins=min_insertion, bowtie2_params=bowtie2_params)

    shutil.rmtree(tmpdirn)

    return samfn
