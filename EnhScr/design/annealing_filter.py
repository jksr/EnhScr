"""Spreadsheet Column Printer

This script allows the user to print to the console all columns in the
spreadsheet. It is assumed that the first row of the spreadsheet is the
location of the columns.

This tool accepts comma separated value files (.csv) as well as excel
(.xls, .xlsx) files.

This script requires that `pandas` be installed within the Python
environment you are running this script in.

This file can also be imported as a module and contains the following
functions:

    * get_spreadsheet_cols - returns the column headers of the file
    * main - the main function of the script
"""


import primer3


class HardFilter:
    """
    Base class for all filters with only pass or fail results

    ...

    Attributes
    ----------
    name : str
        name
    min_Tm : float
        min_Tm
    max_Tm : float
        max_TM

    Methods
    -------
    check(seq)
        check if seq passes this filter
    """
    def check(self) -> bool:
        """Check if seq passes this filter

        Parameters
        ----------
        seq : str
            The file location of the spreadsheet

        Returns
        -------
        bool
            Filtering result.
        """
        pass
    
class SoftFilter:
    def check(self) -> float:
        pass

    
class TmFilter(HardFilter):
    def __init__(self, min_Tm, max_Tm):
        super(TmFilter, self).__init__()
        self.min_Tm = min_Tm
        self.max_Tm = max_Tm
        self.name = 'Tm'
    
    def check(self, seq):
        aseq = str(seq).upper()
        Tm = primer3.calcTm(aseq)
        if Tm < self.min_Tm or Tm > self.max_Tm:
            rtn = False
        else:
            rtn = True
        return rtn


class GCContentFilter(HardFilter):
    def __init__(self, max_gc):
        super(GCContentFilter, self).__init__()
        self.max_gc = max_gc
        self.name = 'GC content'
        
    def check(self, seq):
        aseq = str(seq).upper()
        if isinstance(self.max_gc, float):
            max_gc = int(len(aseq)*self.max_gc)
        elif isinstance(self.max_gc, int):
            max_gc = self.max_gc
        else:
            raise TypeError
            
        rtn = False if aseq.count('CG')>max_gc else True
        return rtn

class PolyXFilter(HardFilter):
    def __init__(self, max_poly_X, X=None):
        super(PolyXFilter, self).__init__()
        self.max_poly_X = max_poly_X
        self.X = X.upper() if X is not None else X
        self.name = 'Poly X' if self.X is None else f'Poly {self.X}'
        
    def check(self, seq):
        aseq = str(seq).upper()
        if isinstance(self.max_poly_X, float):
            max_poly_X = int(len(aseq)*self.max_poly_X)
        elif isinstance(self.max_poly_X, int):
            max_poly_X = self.max_poly_X
        else:
            raise TypeError
        
        max_poly_X += 1
        if self.X is None:
            check_list = ['A'*max_poly_X, 'T'*max_poly_X, 'C'*max_poly_X, 'G'*max_poly_X]
        else:
            check_list = [ x*max_poly_X for x in self.X ]
            
        rtn = sum( [aseq.count(poly) for poly in check_list] )<1

        return rtn 

    
class HairpinFilter(HardFilter):
    def __init__(self, max_tm):
        super(HairpinFilter, self).__init__()
        self.name = 'Hairpin'
        self.max_tm = max_tm
        
    def check(self, seq):
        aseq = str(seq).upper()
        p3rlt = primer3.calcHairpin(aseq, output_structure=True)
        
        rtn = p3rlt.tm < self.max_tm
        return rtn

    
class HomodimerFilter(HardFilter):
    def __init__(self, max_tm):
        super(HomodimerFilter, self).__init__()
        self.name = 'Homodimer'
        self.max_tm = max_tm
        
    def check(self, seq):
        aseq = str(seq).upper()
        p3rlt = primer3.calcHomodimer(aseq, output_structure=True)
        
        rtn = p3rlt.tm < self.max_tm
        return rtn
    