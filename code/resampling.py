import dask
import numpy as np
import xarray as xr


def _resample_iterations_idx(
    init, iterations, dim="member", replace=True, chunk=True, dim_max=None
):
    """Resample over ``dim`` by index ``iterations`` times.

    .. note::
        This is a much faster way to bootstrap than resampling each iteration
        individually and applying the function to it. However, this will create a
        DataArray with dimension ``iteration`` of size ``iterations``. It is probably
        best to do this out-of-memory with ``dask`` if you are doing a large number
        of iterations or using spatial output (i.e., not time series data).

    Args:
        init (xr.DataArray, xr.Dataset): Initialized prediction ensemble.
        iterations (int): Number of bootstrapping iterations.
        dim (str): Dimension name to bootstrap over. Defaults to ``'member'``.
        replace (bool): Bootstrapping with or without replacement. Defaults to ``True``.
        chunk: (bool): Auto-chunk along chunking_dims to get optimal blocksize
        dim_max (int): Number of indices from `dim` to return. Not implemented.

    Returns:
        xr.DataArray, xr.Dataset: Bootstrapped data with additional dim ```iteration```

    """
    #print(init)
    if dask.is_dask_collection(init):
        init = init.chunk({dim: -1})
        init = init.copy(deep=True)
        
    if dim_max is not None and dim_max <= init[dim].size:
        # select only dim_max items
        select_dim_items = dim_max
        new_dim = init[dim].isel({dim: slice(None, dim_max)})
    else:
        select_dim_items = init[dim].size
        new_dim = init[dim]
    #print(init)


    if dask.is_dask_collection(init):
        if chunk:
            chunking_dims = [d for d in init.dims if d not in CLIMPRED_DIMS]
            init = _chunk_before_resample_iterations_idx(
                init, iterations, chunking_dims
            )

    # resample with or without replacement
    if replace:
        idx = np.random.randint(0, init[dim].size, (iterations, select_dim_items))
    elif not replace:
        # create 2d np.arange()
        idx = np.linspace(
            (np.arange(init[dim].size)),
            (np.arange(init[dim].size)),
            iterations,
            dtype="int",
        )
        # shuffle each line
        for ndx in np.arange(iterations):
            np.random.shuffle(idx[ndx])
    idx_da = xr.DataArray(
        idx,
        dims=("iteration", dim),
        coords=({"iteration": range(iterations), dim: new_dim}),
    )
    transpose_kwargs = (
        {"transpose_coords": False} if isinstance(init, xr.DataArray) else {}
    )
    #print(init.transpose(dim, ..., **transpose_kwargs))
    #print(idx_da)
    return xr.apply_ufunc(
        select_bootstrap_indices_ufunc,
        init.transpose(dim, ..., **transpose_kwargs),
        idx_da.rename({'time':'sample'}),
        input_core_dims=[['time'],['iteration',]],
        vectorize=True, 
        output_core_dims=[['iteration']], 
        dask="parallelized",
        output_dtypes=[float],
    ).mean('sample')

def select_bootstrap_indices_ufunc(x, idx):
        """Selects multi-level indices ``idx`` from xarray object ``x`` for all
        iterations."""
        # `apply_ufunc` sometimes adds a singleton dimension on the end, so we squeeze
        # it out here. This leverages multi-level indexing from numpy, so we can
        # select a different set of, e.g., ensemble members for each iteration and
        # construct one large DataArray with ``iterations`` as a dimension.
        return np.moveaxis(x.squeeze()[idx.squeeze().transpose()], 0, -1)