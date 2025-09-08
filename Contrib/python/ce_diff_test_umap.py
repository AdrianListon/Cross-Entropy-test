#!/usr/bin/python3
"""
Python implementation of cross-entropy based statistical test for t-SNE and
UMAP representation comparisons from Roca CP et al. Cell Rep Methods (2023).

Author(s):
    Michael Yao @michael-s-yao
    Allison Chae @allisonjchae

Citation(s):
    [1] Roca CP, Burton OT, Neumann J, Tareen S, Whyte CE, Gergelits V, Veiga
        RV, Humblet-Baron S, Liston A. A cross entropy test allows quantitative
        statistical comparison of t-SNE and UMAP representations. Cell Rep
        Methods 3(1): 100390. (2023). doi: 10.1016/j.crmeth.2022.100390

Licensed under the MIT License. Copyright Main Line Health 2025.
"""
import numpy as np
import torch
from scipy.optimize import brentq
from scipy.stats import ks_2samp
from scipy.stats._stats_py import KstestResult
from tqdm import tqdm
from typing import Optional, Union


def umap_sigma(
    distance_matrix: torch.Tensor,
    verbose: bool = False
) -> torch.Tensor:
    """
    Computes the sigma parameters to fix the perplexity of a Gaussian
    probability distribution. Implements the method from the following URL:
    https://github.com/AdrianListon/Cross-Entropy-test/blob/main/
    CSV-based%20flow%20analysis/00_src/ce_diff_test_umap.r
    Input:
        distance_matrix: a matrix of the pairwise distances of shape NxN.
        verbose: whether to display a verbose progress bar. Default False.
    Returns:
        A vector of length N of the sigma values to use.
    """
    assert distance_matrix.ndim == 2
    assert distance_matrix.size(dim=0) == distance_matrix.size(dim=1)
    n = distance_matrix.size(dim=0)

    sigmas = torch.zeros(
        n, dtype=distance_matrix.dtype, device=distance_matrix.device
    )
    for i, row in enumerate(tqdm(distance_matrix, disable=(not verbose))):
        row_pos = torch.sort(row[row > 0]).values
        if len(row_pos) == 0:
            sigmas[i] = 1.0
            continue
        ss_lower = row_pos[0]

        row_finite = torch.sort(row[torch.isfinite(row)]).values
        if len(row_finite) == 0:
            sigmas[i] = 1.0
            continue
        ss_upper = row_finite[-1]

        def umap_sigma_error(ss: Union[torch.Tensor, float]) -> torch.Tensor:
            p = torch.exp(-1.0 * torch.clamp(row - row.min(), min=0.0) / ss)
            return torch.sum(p) - np.log2(len(p))

        while umap_sigma_error(ss_upper) < 0:
            ss_lower = ss_upper
            ss_upper *= 2.0

        while umap_sigma_error(ss_lower) > 0:
            ss_upper = ss_lower
            ss_lower /= 2.0

        try:
            sigmas[i] = brentq(
                umap_sigma_error,
                ss_lower.item(),
                ss_upper.item(),
                xtol=np.power(torch.finfo(distance_matrix.dtype).eps, 0.25)
            )
        except ValueError:
            sigmas[i] = 0.5 * (ss_lower + ss_upper)

    return torch.clamp(sigmas, min=torch.finfo(distance_matrix.dtype).eps)


def pairwise_distance_matrix(X: torch.Tensor, metric: str) -> torch.Tensor:
    """
    Computes the pairwise distances between points in a dataset.
    Input:
        X: a dataset of points of shape Nxd.
        metric: the distance metric to use. One of [`euclidean`, `cosine`].
    Returns:
        A dataset of pairwise distances of shape NxN.
    """
    if metric.lower() == "euclidean":
        X_norm = torch.square(X).sum(dim=1, keepdim=True)
        D_sq = X_norm + X_norm.T - 2 * X @ X.T
        return torch.sqrt(torch.clamp(D_sq, min=0.0))
    elif metric.lower() == "cosine":
        X_norm = torch.nn.functional.normalize(X, p=2, dim=1)
        return 1.0 - (X_norm @ X_norm.T)
    raise NotImplementedError


@torch.no_grad()
def compute_ce(
    x: torch.Tensor,
    umap_x: torch.Tensor,
    metric: str,
    subsample: float = 1.0,
    random_state: Optional[int] = 2025,
    verbose: bool = False
) -> torch.Tensor:
    """
    Computes the vector of cross entropies as defined in Roca CP et al. Cell
    Rep Methods (2023).
    Input:
        x: a dataset of points of shape Nxd.
        umap_x: a dataset of compressed points of shape Nxd', where d' is the
            chosen UMAP dimension.
        metric: the distance metric to use. One of [`euclidean`, `cosine`].
        subsample: a random subset of the dataset of points to use. By default,
            no subsampling is performed.
        random_state: an optional random seed. Default 2025.
        verbose: whether to display a verbose progress bar.
    Returns:
        The vector of cross entropies of shape N.
    """
    if not np.isclose(subsample, 1.0):
        rng = np.random.default_rng(seed=random_state)
        idxs = rng.choice(len(x), int(subsample * len(x)), replace=False)
        x, umap_x = x[idxs], umap_x[idxs]

    dist_matrix = pairwise_distance_matrix(x, metric=metric)
    sigma = umap_sigma(dist_matrix, verbose=verbose)
    p = torch.exp(
        (-1.0 / sigma.unsqueeze(dim=-1)) * torch.clamp(
            dist_matrix - dist_matrix.min(dim=1).values.unsqueeze(dim=-1),
            min=0.0
        )
    )
    p /= torch.sum(p, dim=-1, keepdim=True) - torch.diag(p).unsqueeze(dim=-1)
    p = (p + p.T) / (2.0 * len(p))

    assert np.isclose((torch.sum(p) - torch.diag(p).sum()).item(), 1.0)
    assert torch.all(torch.isclose(p, p.T))

    umap_dist_matrix = pairwise_distance_matrix(umap_x, metric=metric)
    q = 1.0 / (1.0 + torch.square(umap_dist_matrix))
    q /= torch.sum(q, dim=-1, keepdim=True) - torch.diag(q).unsqueeze(dim=-1)
    q /= float(len(q))

    assert np.isclose((torch.sum(q) - torch.diag(q).sum()).item(), 1.0)

    eps = torch.finfo(q.dtype).eps
    q = torch.clamp(q, min=eps, max=(1.0 - eps))
    return (p * torch.log(q) + ((1.0 - p) * torch.log(1.0 - q))).sum(dim=1)


def ce_diff_test_umap(
    x: torch.Tensor,
    umap_x: torch.Tensor,
    y: torch.Tensor,
    umap_y: torch.Tensor,
    subsample: Union[int, float] = 10000,
    metric: str = "euclidean",
    random_state: Optional[int] = 2025,
    verbose: bool = False
) -> KstestResult:
    """
    Computes the statistical test from Roca CP et al. Cell Rep Methods (2023).
    Input:
        x: a dataset of points of shape Nxd.
        umap_x: a dataset of compressed points of shape Nxd'.
        y: a dataset of points of shape Mxd.
        umap_y: a dataset of compressed points of shape Mxd'.
        subsample: a random subset of the dataset of points to use. By default,
            no subsampling is performed.
        metric: the distance metric to use. One of [`euclidean`, `cosine`].
        random_state: an optional random seed. Default 2025.
        verbose: whether to display a verbose progress bar. Default False.
    Returns:
       The statistical test results.
    """
    ce_x = compute_ce(
        x,
        umap_x,
        metric,
        subsample=(
            subsample
            if isinstance(subsample, float)
            else float(subsample) / len(x)
        ),
        random_state=random_state,
        verbose=verbose
    )
    ce_y = compute_ce(
        y,
        umap_y,
        metric,
        subsample=(
            subsample
            if isinstance(subsample, float)
            else float(subsample) / len(y)
        ),
        random_state=random_state,
        verbose=verbose
    )
    return ks_2samp(
        ce_x.detach().cpu().numpy(),
        ce_y.detach().cpu().numpy(),
        alternative="two-sided",
        method="auto"
    )
