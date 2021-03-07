import numpy as np

class DesignSelector:
    def select(self, *args, **kwargs) -> str:
        pass
    
    
class RandomDesignSelector(DesignSelector):
    def __init__(self, *args, **kwargs):
        super(RandomDesignSelector, self).__init__()
    
    def select(self, candidate_list):
        return np.random.choice(list(candidate_list))

    
class MostCenteredDesignSelector(DesignSelector):
    def __init__(self, pool_designer, *args, **kwargs):
        super(MostCenteredDesignSelector, self).__init__()
        self.pool_designer = pool_designer
    
    def select(self, candidate_list):
        candidate_list = np.array(list(candidate_list))
        
        centerness = []
        for candid in candidate_list:
            cand = self.pool_designer.get_design_candidate(candid)
            c = len(cand.seq)/2 - ( cand.tgt_offset + len(cand._tgt_designer.tgt_seq)/2 )
            centerness.append( abs(c) )
        
        return candidate_list[np.argsort(centerness)][0]
    