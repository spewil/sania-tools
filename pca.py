from sklearn.decomposition import PCA
from sequence import sc


def compute_pca(stack):
    pca = PCA()
    pca.fit(stack)
    components = pca.components_
    projections = components[:2] @ stack.T
    return projections, components


bin_projections, bin_components = compute_pca(sc.bin_stack)
num_projections, num_components = compute_pca(sc.num_stack)
mut_projections, mut_components = compute_pca(sc.mutation_stack)
