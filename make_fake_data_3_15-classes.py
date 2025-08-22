
# simulate a strip: sort of
# invert X to be 255 minus the intensity DONE
# use pandas to output labeled RGB data DONE
# Github push DONE
# HCA DONE
# ASCII text for headers (KIM: ask Mark) DONE
# 2. PCA with 3D ellipses: DONE
# allowed for multiple colors DONE
# 3. NEXT 15 classes strip (get RGB values for 15 classes from Josselyn)
# Goal of 15 classes: adopt code from 3 classes, make educated guesses at binding affinities, see if it recapitulates PCA data


# need NP colors (3 of them, which one is Innova?)
# need binding affinities

# cooperativity?



# look in literature for stuff about actual noise of antibody-antigen binding
# We need a name for this thing. Currently it's just called "the simulator"



# Where are we at F 8/22/2025: 
# Pasted in make_fake_data_2.py
# changed how colors for plotting (random hexadecimal) so can extend to 15 classes
# obtained direct test data for 15 classes as stand-in for Ab binding profiles, made educated guesses for missing numbers
# started changing immunoprobe colors and binding affinities, midway through. Code currently broken
# Things to fix: importing the Ab binding profiles, and getting the row names (i.e., the antigens) Line 89
# Still to do: have the RGB values of the 3 immunoprobes from JMC in excel file, need to define the immunoprobe colors







import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy.stats import chi2
import random

################################################################################
################################################################################
 
# ██╗███╗   ██╗██╗████████╗██╗ █████╗ ██╗     ██╗███████╗███████╗
# ██║████╗  ██║██║╚══██╔══╝██║██╔══██╗██║     ██║╚══███╔╝██╔════╝
# ██║██╔██╗ ██║██║   ██║   ██║███████║██║     ██║  ███╔╝ █████╗  
# ██║██║╚██╗██║██║   ██║   ██║██╔══██║██║     ██║ ███╔╝  ██╔══╝  
# ██║██║ ╚████║██║   ██║   ██║██║  ██║███████╗██║███████╗███████╗
# ╚═╝╚═╝  ╚═══╝╚═╝   ╚═╝   ╚═╝╚═╝  ╚═╝╚══════╝╚═╝╚══════╝╚══════╝

################################################################################
################################################################################

# =================================
# Make table of nanoparticle colors
# =================================

color_channels = ("R", "G", "B",)
NP_names = ("Innova", "Red", "Blue")
# NP_colors = [
# #  Innova     Red    Blue
#     [122.15, 50.39, ],  # R
#     [119.85, 67.67, ],  # G
#     [110.61, 53.24, ],] # B


# NP_colors = pd.DataFrame(
#     NP_colors, 
#     index = color_channels, 
#     columns = NP_names,)

################################################################################

# ==========================================
# Make table of idealized binding affinities
# ==========================================
# antibody binding profiles are 5 columns X 15 rows

xls = pd.ExcelFile('15-class_Ab_binding_profiles.xlsx')
df1 = pd.read_excel(xls, 'Ab_profiles_08222025')
print("hello")
print((df1.columns.tolist()[1:]))
print((df1.index.tolist()[1:]))
print("goodbye")


antigens = ("Control", "Alpha", "BA.5", "BA.1",)
antibodies = (list(df1.head(5))[1:])


mean_affinities = [
#   D001  R007
    [0.0, 0.0],  # Control
    [0.6, 0.7],  # Alpha
    [0.0, 0.4],  # BA.5
    [0.0, 1.0],] # BA.1

mean_affinities = pd.DataFrame(
    mean_affinities, 
    index = antigens, 
    columns = antibodies,)

################################################################################

# ==========================================================================
# Make tables of antibody pair interferences, initialize with 1.0 as default
# ==========================================================================

# 1.0 = no interference

interferences = {
    antigen: pd.DataFrame(
        data = 1.0, # default
        index = antibodies,
        columns = antibodies)
    for antigen in antigens}

# Declare specific pairs that do interfere
# Example: declare a symmetric interference of 0.5 for antibodies A and B:
# interferences['some antigen'].at['antibody A', 'antibody B'] = 0.5
# interferences['some antigen'].at['antibody B', 'antibody A'] = 0.5

################################################################################

# ==================================
# Define test lines and immunoprobes
# ==================================

test_lines = ("D001", "R007",)
immunoprobe_NPs = {"D001": "Red", "R007": "Blue",}
immunoprobes = immunoprobe_NPs.keys()

# sanity checks
for test_line in test_lines:
    assert test_line in antibodies
for probe in immunoprobes:
    assert probe in antibodies
    NP_name = immunoprobe_NPs[probe]
    assert NP_name in NP_names

################################################################################




# ====================================================
# Define distribution of noise about each cluster mean 
# ====================================================

# define target variance for each antigen
variances = {
    "Control": 0.001, 
    "Alpha": 0.001,
    "BA.5": 0.001,
    "BA.1": 0.001,}

# TODO write comment
noise_scaling_spot = 0.5
noise_scaling_sandwich = 0.2

# Get number of dimensions in the antibody space
ABs_in_use = tuple(set(test_lines) | set(immunoprobes))
dimensions = len(ABs_in_use)

# Spherical (i.i.d.) noise: use diagonal covariance matrices (scaled identity)
cov_matrices = {
    antigen: np.identity(dimensions) * variances[antigen]
    for antigen in variances}

################################################################################


# ==============================
# Narrow down to antigens in use
# ==============================


antigens_in_use = tuple(list(antigens))

antigens_in_use_indices = {}
for antigen in antigens_in_use:
    antigens_in_use_indices[antigen]=antigens_in_use.index(antigen)




################################################################################
################################################################################

# ███████╗ █████╗ ██╗  ██╗███████╗    ██████╗  █████╗ ████████╗ █████╗ 
# ██╔════╝██╔══██╗██║ ██╔╝██╔════╝    ██╔══██╗██╔══██╗╚══██╔══╝██╔══██╗
# █████╗  ███████║█████╔╝ █████╗      ██║  ██║███████║   ██║   ███████║
# ██╔══╝  ██╔══██║██╔═██╗ ██╔══╝      ██║  ██║██╔══██║   ██║   ██╔══██║
# ██║     ██║  ██║██║  ██╗███████╗    ██████╔╝██║  ██║   ██║   ██║  ██║
# ╚═╝     ╚═╝  ╚═╝╚═╝  ╚═╝╚══════╝    ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚═╝  ╚═╝
                                                                     
################################################################################
################################################################################

# ==========================
# Math to generate fake data
# ==========================

def add_noise_to_affins(
        antigen: str, 
        affins: pd.DataFrame,
        scaling: float = 1.0,
        ) -> pd.DataFrame:
    noisy_affins = affins.copy()
    deviations = np.random.multivariate_normal(
        mean = np.zeros(cov_matrices[antigen].shape[0]), 
        cov = cov_matrices[antigen])
    i = 0
    for antibody in ABs_in_use:
        noisy_affins.loc[antigen, antibody] += (deviations[i] * scaling)
        i += 1
    return noisy_affins

def sandwich_color(
        antigen: str, 
        test_line: str, 
        immunoprobe: str, 
        affins: pd.DataFrame,
        add_noise: bool,
        ) -> pd.Series:
    '''Generate color channel data for a single sandwich'''
    if add_noise: affins = add_noise_to_affins(antigen, affins, noise_scaling_sandwich)
    NP_color = NP_colors[immunoprobe_NPs[immunoprobe]].copy()
    immunoprobe_affinity = affins.at[antigen, immunoprobe]
    test_line_affinity = affins.at[antigen, test_line]
    sandwich_affinity = test_line_affinity * immunoprobe_affinity
    interference = interferences[antigen].at[test_line, immunoprobe]
    return NP_color * sandwich_affinity * interference

def spot_color(
        antigen: str, 
        test_line: str, 
        affins: pd.DataFrame, 
        add_noise: bool,
        ) -> pd.Series:
    '''Generate combined color channel data for all sandwiches at a spot'''
    color = pd.Series([255.0]*len(color_channels), index = color_channels, dtype = float,)
    if add_noise: affins = add_noise_to_affins(antigen, affins, noise_scaling_spot)
    for probe in immunoprobes: 
        color -= sandwich_color(antigen, test_line, probe, affins, add_noise)
    return color

def strip_colors(antigen: str, add_noise: bool) -> pd.DataFrame:
    '''Generate a full-strip dataset of all the spot colors'''
    affins = mean_affinities
    if add_noise: affins = add_noise_to_affins(antigen, affins)
    result = pd.DataFrame(index = test_lines, columns = color_channels, dtype=float,)
    for test_line in test_lines: 
        result.loc[test_line] = spot_color(antigen, test_line, affins, add_noise)
    return result

def generate_strips(strips_per_antigen_in_use=1,add_noise = False) -> tuple[np.ndarray, list]:
    ''' Generate strips for each antigen in use'''
    all_data = []
    all_labels = []


    for antigen in antigens_in_use:
        for _ in range(strips_per_antigen_in_use):
            strip = strip_colors(antigen, add_noise=True)
            all_data.append(strip)
            all_labels.append(antigen)

    # Convert each strip to a 6D RGB vector: [D001_R, D001_G, D001_B, R007_R, R007_G, R007_B]
    X = []
    for strip in all_data:
        D001_rgb = strip.loc["D001", ["R", "G", "B"]].values
        R007_rgb = strip.loc["R007", ["R", "G", "B"]].values
        X.append(np.concatenate([D001_rgb, R007_rgb]))
    return (np.array(X), all_labels)

# 'strips' is the raw colorimetric data for all generated strips (all antigens in use)
# 'strip_labels' corresponds to strips, provides the antigen name for each entry in strips
strips,strip_labels = generate_strips(strips_per_antigen_in_use=100, add_noise=True)


################################################################################
################################################################################

#  █████╗ ███╗   ██╗ █████╗ ██╗  ██╗   ██╗███████╗██╗███████╗
# ██╔══██╗████╗  ██║██╔══██╗██║  ╚██╗ ██╔╝██╔════╝██║██╔════╝
# ███████║██╔██╗ ██║███████║██║   ╚████╔╝ ███████╗██║███████╗
# ██╔══██║██║╚██╗██║██╔══██║██║    ╚██╔╝  ╚════██║██║╚════██║
# ██║  ██║██║ ╚████║██║  ██║███████╗██║   ███████║██║███████║
# ╚═╝  ╚═╝╚═╝  ╚═══╝╚═╝  ╚═╝╚══════╝╚═╝   ╚══════╝╚═╝╚══════╝

################################################################################
################################################################################

# Perform PCA to reduce from 6D to 3D
pca = PCA(n_components=3)
X_pca = pca.fit_transform(strips)


# print data as test
# print(X_pca)

# cluster_points: break up all_data into subsets by antigen
cluster_points = []
for antigen in antigens_in_use:
    mask = [label == antigen for label in strip_labels]
    cluster_points.append(X_pca[mask])


# cluster_means: the mean of the points in each cluster
cluster_means=[]
for cluster in cluster_points:
    mean_of_points = np.mean(cluster, axis=0)
    cluster_means.append(mean_of_points)


# cluster_covs: the covariance matrices of the points in each cluster
cluster_covs=[]

for cluster in cluster_points:
    cov_of_points = np.cov(cluster.T)
    cluster_covs.append(cov_of_points)



################################################################################
################################################################################

# ██████╗ ██╗      ██████╗ ████████╗    ██████╗  █████╗ ████████╗ █████╗ 
# ██╔══██╗██║     ██╔═══██╗╚══██╔══╝    ██╔══██╗██╔══██╗╚══██╔══╝██╔══██╗
# ██████╔╝██║     ██║   ██║   ██║       ██║  ██║███████║   ██║   ███████║
# ██╔═══╝ ██║     ██║   ██║   ██║       ██║  ██║██╔══██║   ██║   ██╔══██║
# ██║     ███████╗╚██████╔╝   ██║       ██████╔╝██║  ██║   ██║   ██║  ██║
# ╚═╝     ╚══════╝ ╚═════╝    ╚═╝       ╚═════╝ ╚═╝  ╚═╝   ╚═╝   ╚═╝  ╚═╝
                                                                       
################################################################################
################################################################################

# Ellipsoid code from https://github.com/CircusMonkey/covariance-ellipsoid/blob/master/ellipsoid.py 
# August 2025

def get_cov_ellipsoid(cov, mu=np.zeros((3)), nstd=3):
    """
    Return the 3d points representing the covariance matrix
    cov centred at mu and scaled by the factor nstd.

    Plot on your favourite 3d axis. 
    Example 1:  ax.plot_wireframe(X,Y,Z,alpha=0.1)
    Example 2:  ax.plot_surface(X,Y,Z,alpha=0.1)
    """
    assert cov.shape==(3,3)

    # Find and sort eigenvalues to correspond to the covariance matrix
    eigvals, eigvecs = np.linalg.eigh(cov)
    idx = np.sum(cov,axis=0).argsort()
    eigvals_temp = eigvals[idx]
    idx = eigvals_temp.argsort()
    eigvals = eigvals[idx]
    eigvecs = eigvecs[:,idx]

    # Set of all spherical angles to draw our ellipsoid
    n_points = 100
    theta = np.linspace(0, 2*np.pi, n_points)
    phi = np.linspace(0, np.pi, n_points)

    # Width, height and depth of ellipsoid
    rx, ry, rz = nstd * np.sqrt(eigvals)

    # Get the xyz points for plotting
    # Cartesian coordinates that correspond to the spherical angles:
    X = rx * np.outer(np.cos(theta), np.sin(phi))
    Y = ry * np.outer(np.sin(theta), np.sin(phi))
    Z = rz * np.outer(np.ones_like(theta), np.cos(phi))

    # Rotate ellipsoid for off axis alignment
    old_shape = X.shape
    # Flatten to vectorise rotation
    X,Y,Z = X.flatten(), Y.flatten(), Z.flatten()
    X,Y,Z = np.matmul(eigvecs, np.array([X,Y,Z]))
    X,Y,Z = X.reshape(old_shape), Y.reshape(old_shape), Z.reshape(old_shape)
   
    # Add in offsets for the mean
    X = X + mu[0]
    Y = Y + mu[1]
    Z = Z + mu[2]
    
    return X,Y,Z





 # Setup the plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
plt.xlabel('PC1')
plt.ylabel('PC2')

# plot data points and ellipses


nstd = 2 # number of standard deviations of ellipsoid; determines the ellipsoid volume

for antigen in antigens_in_use:
    i = antigens_in_use_indices[antigen]
    points = cluster_points[i]
    mean = cluster_means[i]
    cov = cluster_covs[i]
    r = lambda: random.randint(0,255)
    plot_color = str('#%02X%02X%02X' % (r(),r(),r()))
    X1,Y1,Z1 = get_cov_ellipsoid(cov, mean, nstd)
    ax.plot_wireframe(X1,Y1,Z1, color = plot_color, alpha=0.1)
    ax.scatter(points[:,0],points[:,1],points[:,2],c =plot_color)




################################################################################
################################################################################

# ██████╗ ██████╗ ██╗███╗   ██╗████████╗     ██████╗ ██████╗ ██╗      ██████╗ ██████╗ ███████╗
# ██╔══██╗██╔══██╗██║████╗  ██║╚══██╔══╝    ██╔════╝██╔═══██╗██║     ██╔═══██╗██╔══██╗██╔════╝
# ██████╔╝██████╔╝██║██╔██╗ ██║   ██║       ██║     ██║   ██║██║     ██║   ██║██████╔╝███████╗
# ██╔═══╝ ██╔══██╗██║██║╚██╗██║   ██║       ██║     ██║   ██║██║     ██║   ██║██╔══██╗╚════██║
# ██║     ██║  ██║██║██║ ╚████║   ██║       ╚██████╗╚██████╔╝███████╗╚██████╔╝██║  ██║███████║
# ╚═╝     ╚═╝  ╚═╝╚═╝╚═╝  ╚═══╝   ╚═╝        ╚═════╝ ╚═════╝ ╚══════╝ ╚═════╝ ╚═╝  ╚═╝╚══════╝

################################################################################
################################################################################

# ChatGPT code

# ========================================
# Generate a simulated image of a strip 
# ========================================

strips, strip_labels=generate_strips()
strips = list(strips)


column_names = []
for test_line in test_lines:
    column_names.append(f"R_{test_line}")
    column_names.append(f"G_{test_line}")
    column_names.append(f"B_{test_line}")


#column_names= ("R_D001", "G_D001", "B_D001", "R_R007" ,"G_R007","B_R007" )
#row_names_classes = ("Control", "alpha", "BA.5", "BA.1")


strips = pd.DataFrame(
    strips, 
    index = antigens_in_use, 
    columns = column_names,)

print(strips)




################################################################################
################################################################################

# ██████╗ ██╗      ██████╗ ████████╗    ██╗  ██╗ ██████╗ █████╗ 
# ██╔══██╗██║     ██╔═══██╗╚══██╔══╝    ██║  ██║██╔════╝██╔══██╗
# ██████╔╝██║     ██║   ██║   ██║       ███████║██║     ███████║
# ██╔═══╝ ██║     ██║   ██║   ██║       ██╔══██║██║     ██╔══██║
# ██║     ███████╗╚██████╔╝   ██║       ██║  ██║╚██████╗██║  ██║
# ╚═╝     ╚══════╝ ╚═════╝    ╚═╝       ╚═╝  ╚═╝ ╚═════╝╚═╝  ╚═╝

################################################################################
################################################################################


################################################################################
# New HCA code block: from GEMINI
################################################################################

import scipy.cluster.hierarchy as sch

# Perform hierarchical clustering on the same 6D data (X)
# We use the 'ward' method, which is good for minimizing the variance
# within each cluster, and 'euclidean' metric is the default and a good choice.
Z = sch.linkage(strips, method='ward')

# Create a figure for the dendrogram
plt.figure(figsize=(10, 7))
plt.title("Hierarchical Cluster Analysis Dendrogram")
plt.xlabel("Sample Index")
plt.ylabel("Distance")

# Create the dendrogram
sch.dendrogram(
    Z,
    labels=strip_labels,
    leaf_rotation=90.,  # Rotate the leaf labels for better readability
    leaf_font_size=10.,
)

# You can add a horizontal line to "cut" the dendrogram and
# visually determine the clusters. A value of 't' for the fcluster function.
# plt.axhline(y=10, c='k', linestyle='--') # Example cut-off line

plt.tight_layout()
plt.show()

# To get the cluster assignments from the dendrogram, you can 'cut' it
# at a specific distance threshold. For example, a threshold of 10.
# You can adjust this value based on what you see in the dendrogram.
# from scipy.cluster.hierarchy import fcluster
# cluster_labels_hca = fcluster(Z, t=10, criterion='distance')
# print("HCA Cluster Labels:", cluster_labels_hca)

# You could also get a specific number of clusters.
# cluster_labels_hca_3 = fcluster(Z, t=3, criterion='maxclust')
# print("HCA Cluster Labels (3 clusters):", cluster_labels_hca_3)