import numpy as np

def allcomb(*varargin):
    NC = len(varargin)

    if NC == 0:
        return np.zeros((0, 0))

    if isinstance(varargin[-1], str) and (varargin[-1].lower() == 'matlab' or varargin[-1].lower() == 'john'):
        NC -= 1
        ii = list(range(NC))
    else:
        ii = list(range(NC - 1, -1, -1))

    if any(len(v) == 0 for v in [varargin[i] for i in ii]):
        import warnings
        warnings.warn('ALLCOMB:EmptyInput Empty inputs result in an empty output.')
        return np.zeros((0, NC))

    if NC > 1:
        isCellInput = [isinstance(v, (list, tuple)) for v in varargin[:NC]]
        if any(isCellInput):
            if not all(isCellInput):
                raise ValueError('ALLCOMB:InvalidCellInput For cell input, all arguments should be cell arrays.')

            ix = [list(range(len(c))) for c in varargin[:NC]]

            grids = np.meshgrid(*[ix[i] for i in ii], indexing='ij')
            # reorder grids to original order
            grids_ordered = [None] * NC
            for idx, val in zip(ii, grids):
                grids_ordered[idx] = val

            A = []
            for k in range(NC):
                vals = [varargin[k][i] for i in grids_ordered[k].flatten()]
                A.append(vals)
            A = list(zip(*A))
            return A
        else:
            grids = np.meshgrid(*[varargin[i] for i in ii], indexing='ij')
            grids_ordered = [None] * NC
            for idx, val in zip(ii, grids):
                grids_ordered[idx] = val
            A = np.stack(grids_ordered, axis=-1).reshape(-1, NC)
            return A
    elif NC == 1:
        return np.array(varargin[0]).reshape(-1, 1)
    else:
        return np.zeros((0, 0))