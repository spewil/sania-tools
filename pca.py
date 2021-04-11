from sklearn.decomposition import PCA
from sequence import sc


def compute_pca(stack):
    pca = PCA()
    pca.fit(stack)
    components = pca.components_
    projections = components[:2] @ stack.T
    variance = pca.explained_variance_ratio_
    return projections, components, variance


bin_projections, bin_components, bin_variance = compute_pca(sc.bin_stack)
num_projections, num_components, num_variance = compute_pca(sc.num_stack)
mut_projections, mut_components, mut_variance = compute_pca(sc.mutation_stack)
