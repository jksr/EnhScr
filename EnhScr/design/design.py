import pandas as pd
import pyfaidx
import networkx as nx
import numpy as np
from ..tools.blast import mmseqs2_from_lists
from Bio.Seq import Seq

class DesignParameters:
    def __init__(self, fa_path, oligo_length, n_oligos_per_target, 
                 anneal_length, annealing_filters, cushion_length):
        
        self.faidx = pyfaidx.Faidx(fa_path)
        self.olg_len = oligo_length
        self.n_olg_per_tgt = n_oligos_per_target        
        self.anl_len = anneal_length
        self.anl_filters = annealing_filters
        self.csh_len = cushion_length
        
        self.tot_len = oligo_length * n_oligos_per_target - anneal_length * (n_oligos_per_target-1)
     
    
    
class OligoDesign:
    def __init__(self, candidate_id, target_designer, anneal_starts):
        self.id = f'{target_designer.id}.{candidate_id}'
        self._tgt_designer = target_designer
        self._anneal_starts = anneal_starts # these positions are corresponding to the scan region in TargetDesigner
        self.oligos = []
        self.oligos_on_pos_strand = []
        self.annealings = []
        self.seq = ''
        self.tgt_offset = None
        
        self._assemble_seq()
        self._reverse_complement_even_number_oligos()
        
    def _assemble_seq(self):
        anl_len = self._tgt_designer.params.anl_len
        olg_len = self._tgt_designer.params.olg_len
        olg_n = self._tgt_designer.params.n_olg_per_tgt
        adt_lft =self._tgt_designer.adt_lft
        adt_rgt =self._tgt_designer.adt_rgt
        
        
        seq_start = self._anneal_starts[0] + anl_len - olg_len
        assert seq_start>=0,f'{self._anneal_starts[0]} + {anl_len} - {olg_len} = {seq_start}'
        seq_end = self._anneal_starts[-1] + olg_len
        assert seq_end>0
        
        self.seq = self._tgt_designer.scn_seq[seq_start:seq_end]
        
        self.seq = adt_lft + self.seq[len(adt_lft):]
        self.seq = self.seq[:-len(adt_rgt)] + adt_rgt
        
        for i in range(olg_n):
            start = i*(olg_len-anl_len)
            self.oligos.append( self.seq[start:start+olg_len] )
            self.oligos_on_pos_strand.append( self.seq[start:start+olg_len] )
            
            if i<olg_n-1:
                self.annealings.append( self.seq[start+olg_len-anl_len:start+olg_len] )

        self.tgt_offset = self.seq.index(self._tgt_designer.tgt_seq)
        
    def _reverse_complement_even_number_oligos(self):
        for i in range(1,len(self.oligos),2):
            self.oligos[i] = str(Seq(self.oligos[i]).reverse_complement())



class TargetDesigner:
    def __init__(self, id, chrom, start, end, design_params, adaptor_left, adaptor_right ):
        self.id = id
        self.tgt_chrom = chrom
        self.tgt_start = start
        self.tgt_end = end
        self.tgt_len = end-start
        
        self.params = design_params

        self.adt_lft = adaptor_left
        self.adt_rgt = adaptor_right
        self.adt_lft_len = len(adaptor_left)
        self.adt_rgt_len = len(adaptor_right)
        assert self.params.csh_len >= max(self.adt_lft_len, self.adt_rgt_len)
        
        self.tgt_seq = self.params.faidx.fetch(chrom, start, end).seq

        slop = self.params.tot_len - self.tgt_len - self.params.csh_len
        self.scn_seq = self.params.faidx.fetch(chrom, start-slop, end+slop).seq
        
        self.anneal_df = None
        self.candidates = []
    
#     def __repr__(self):
#         return 'haha'
    
    def __str__(self):
        return self.__repr__()
    
    def _filter_annealing_position(self):
        offsets = []
        cands_flt_rlts=[]
        seqs = []
        
        offset_start = self.params.olg_len - self.params.anl_len
        offset_end = len(self.scn_seq) - self.params.olg_len
        for offset in range(offset_start, offset_end):
            seq = self.scn_seq[offset:offset+self.params.anl_len]
            flt = [ flt.check(seq) for flt in self.params.anl_filters ]
            
            offsets.append(offset)
            cands_flt_rlts.append(flt)
            seqs.append(seq)
            
        anneal_df = pd.DataFrame(dict(offset=offsets,
                                       target_id=self.id,
                                       anneal_seq=seqs,
                                       ))
        anneal_df = anneal_df.join(pd.DataFrame(data = cands_flt_rlts,
                                                  columns = [f'{flt.name} :: pass' for flt in self.params.anl_filters],))
        anneal_df['Pass'] = np.all(cands_flt_rlts, axis=1)
        anneal_df.index = anneal_df['target_id']+'.'+anneal_df['offset'].astype(str)
        
        self.anneal_df = anneal_df
        return self.anneal_df
    
    
    
    def _filter_candidate_sets(self):
        anl_sep_len = self.params.olg_len - self.params.anl_len
        anl_base_pos = np.array([ i*anl_sep_len for i in range(self.params.n_olg_per_tgt-1) ])
        
        self.candidates = []
        
        for i in range(len(self.anneal_df)):
            anl_cand_pos = i + anl_base_pos
            if max(anl_cand_pos)>=len(self.anneal_df):
                break
            
            cand_df = self.anneal_df.iloc[anl_cand_pos]
            if cand_df['Pass'].all():
                self.candidates.append(
                    OligoDesign( len(self.candidates), self, cand_df['offset'].tolist()) )
        

    def design(self):
        self._filter_annealing_position()
        self._filter_candidate_sets()
        
    def get_design_candidate(self, cand_id):
        if isinstance(cand_id, int):
            pass
        elif isinstance(cand_id, str):
            cand_id = int(cand_id.rsplit('.',1)[1])
        else:
            raise NotImplementedError
        return self.candidates[cand_id]
        
        

class PoolDesigner:
    def __init__(self, target_df, adatpor_df, design_params, verbose=0):
        
        self.tgt_df = target_df
        self.adt_df = adatpor_df.set_index('name')
        self.params = design_params
        self.verbose = verbose
        
        self.tgt_designers = []
        self.cand_to_tgt_designer = dict()
        self.conflict_graph = None

        
    def design(self, override=False):
        if override or len(self.tgt_designers)==0:
            self._individual_target_design()
        if override or self.conflict_graph is None:
            self._detect_design_conflicts()
    
    
    def _individual_target_design(self):
        self.tgt_designers = []
        for idx,rec in self.tgt_df.iterrows():
            tgt_designer = TargetDesigner(rec['name'], rec['chr'], rec['start'], rec['end'], self.params,
                                          self.adt_df.loc[rec['name'], 'left'], self.adt_df.loc[rec['name'], 'right']
                                         )
            tgt_designer.design()
            
            self.tgt_designers.append(tgt_designer)
            for cand in tgt_designer.candidates:
                self.cand_to_tgt_designer[cand.id] = tgt_designer
                
            
    def _detect_design_conflicts(self):
        anneal_df = []
        target_df = []

        for tdr in self.tgt_designers:
            for cand in tdr.candidates:
                target_df.append([cand.id, cand.seq])
                for i,anl_seq in enumerate(cand.annealings):
                    anneal_df.append([cand.id, f'{cand.id}.{i}', anl_seq])
                ## TODO add adaptor to anneal_df
        
        
        anneal_df = pd.DataFrame(anneal_df, columns=['target_id','anneal_id','anneal_seq'])
        target_df = pd.DataFrame(target_df, columns=['target_id','final_seq'])    

        anlvanl = mmseqs2_from_lists(anneal_df['anneal_seq'], anneal_df['anneal_id'],
                             anneal_df['anneal_seq'], anneal_df['anneal_id'])
        anlvtgt = mmseqs2_from_lists(anneal_df['anneal_seq'], anneal_df['anneal_id'],
                             target_df['final_seq'], target_df['target_id'])
        
        anlvanl = anlvanl[['query','target']]
        anlvanl.loc[:,'query'] = anlvanl['query'].str.rsplit('.',1).str[0]
        anlvanl.loc[:,'target'] = anlvanl['target'].str.rsplit('.',1).str[0]
        
        anlvtgt = anlvtgt[['query','target']]
        anlvtgt.loc[:,'query'] = anlvtgt['query'].str.rsplit('.',1).str[0]
        
        
        self.design_conflicts = pd.concat([anlvanl,anlvtgt]).drop_duplicates()
        self.design_conflicts = self.design_conflicts[self.design_conflicts['query']!=self.design_conflicts['target']]
        
        
        self.conflict_graph = nx.Graph()
        self.conflict_graph.add_nodes_from(self.cand_to_tgt_designer.keys())
        self.conflict_graph.add_edges_from(self.design_conflicts.values)
        

    def construct_oligo_pools(self, n_pools, selector, exclude_pooled_designs = True, exclude_pooled_targets = True):
        graph = self.conflict_graph.copy()
        
        pools = []
        for _ in range(n_pools):
            currpool = []
            sub_graphs = list(nx.connected_component_subgraphs(graph))
            for sg in sub_graphs:
                candset = selector.select(sg.nodes())
                currpool.append(candset)
            pools.append(currpool)

            if exclude_pooled_designs:
                graph.remove_nodes_from(currpool)

            if exclude_pooled_targets:
                pooled_targets = set( self.cand_to_tgt_designer[cand_id] for cand_id in currpool )
                to_remove = []
                for tgt in pooled_targets:
                    to_remove += [cand.id for cand in tgt.candidates]
                graph.remove_nodes_from(to_remove)

        return pools

    
    def get_design_candidate(self, cand_id):
        return self.cand_to_tgt_designer[cand_id].get_design_candidate(cand_id)
        