# simulate a strip: sort of
# invert X to be 255 minus the intensity DONE
# use pandas to output labeled RGB data
# HCA
# Github push
# PCA with 3D ellipses

# Next: 15 classes strip



# cooperativity?



# look in literature for stuff about actual noise of antibody-antigen binding
# We need a name for this thing. Currently it's just called "the simulator"



import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
from scipy.stats import chi2

################################################################################

# =================================
# Make table of nanoparticle colors
# =================================

color_channels = ("R", "G", "B",)
NP_names = ("Blue", "Red",)
NP_colors = [
#   Blue     Red
    [122.15, 50.39],  # R
    [119.85, 67.67],  # G
    [110.61, 53.24],] # B


NP_colors = pd.DataFrame(
    NP_colors, 
    index = color_channels, 
    columns = NP_names,)

################################################################################

# ==========================================
# Make table of idealized binding affinities
# ==========================================

antigens = ("Control", "Alpha", "BA.5", "BA.1",)
antibodies = ("D001", "R007")
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

################################################################################

# ===========================================================
# Generate the cluster mean or "ideal value" for each antigen
# ===========================================================

cluster_means = {
    antigen: strip_colors(antigen, False) 
    for antigen in antigens}

################################################################################

# =============================================================
# Generate a fake data cluster around the mean for each antigen
# =============================================================

# data_for_antigen




# PCA
# try 3d and 2d




    
print(strip_colors("BA.1", True))
print(strip_colors("BA.5", True))

# Confirm that, when there is no interference, our fake data math creates
# spots that are all "the same color," i.e., have colinear RGB vectors
print(81.0173/69.4434)
print(87.1479/74.6982)
print(76.5597/65.6226)






################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

# ChatGPT code

# ========================================
# ChatGPT code: Generate and visualize data
# ========================================



# Store flattened data and labels
all_data = []
all_labels = []

n_samples = 100

for antigen in antigens:
    for _ in range(n_samples):
        strip = strip_colors(antigen, add_noise=True)
        all_data.append(strip)
        all_labels.append(antigen)

# Convert each strip to a 6D RGB vector: [D001_R, D001_G, D001_B, R007_R, R007_G, R007_B]
X = []
for strip in all_data:
    D001_rgb = strip.loc["D001", ["R", "G", "B"]].values
    R007_rgb = strip.loc["R007", ["R", "G", "B"]].values
    X.append(np.concatenate([D001_rgb, R007_rgb]))
X = np.array(X)

# Perform PCA to reduce from 6D to 2D
pca = PCA(n_components=2)
X_pca = pca.fit_transform(X)

# Map antigen labels to display colors
label_colors = {
    "Control": "gray",
    "Alpha": "red",
    "BA.5": "blue",
    "BA.1": "green",
}

# Ellipse parameters
conf_level = 0.95
chi2_val = chi2.ppf(conf_level, df=2)  # Mahalanobis radius for 95% CI

def plot_confidence_ellipse(data, ax, edgecolor='black'):
    """
    Plot a 95% confidence ellipse for a 2D dataset on the given axis.
    """
    if len(data) < 2:
        return  # Not enough data for ellipse
    mean = np.mean(data, axis=0)
    cov = np.cov(data.T)

    # Eigen-decomposition
    eigvals, eigvecs = np.linalg.eigh(cov)
    order = eigvals.argsort()[::-1]
    eigvals, eigvecs = eigvals[order], eigvecs[:, order]

    # Width and height are 2 * sqrt of chi2_val * eigenvalues
    width, height = 2 * np.sqrt(chi2_val * eigvals)
    angle = np.degrees(np.arctan2(*eigvecs[:, 0][::-1]))

    ellipse = Ellipse(
        xy=mean,
        width=width,
        height=height,
        angle=angle,
        edgecolor=edgecolor,
        facecolor='none',
        lw=2,
    )
    ax.add_patch(ellipse)

# Plotting
fig, ax = plt.subplots(figsize=(8, 6))

for antigen in antigens:
    indices = [i for i, label in enumerate(all_labels) if label == antigen]
    cluster_data = X_pca[indices]
    ax.scatter(
        cluster_data[:, 0],
        cluster_data[:, 1],
        color=label_colors.get(antigen, "black"),
        alpha=0.6,
        s=40,
        label=antigen,
    )
    plot_confidence_ellipse(cluster_data, ax, edgecolor=label_colors[antigen])

# Final touches
ax.set_xlabel("PCA Component 1")
ax.set_ylabel("PCA Component 2")
ax.set_title("PCA of Simulated RGB Spot Colors\nwith 95% Confidence Ellipses")
ax.legend()
ax.grid(True)
plt.tight_layout()
plt.show()



################################################################################
################################################################################
################################################################################
################################################################################
################################################################################



# ChatGPT code

# ========================================
# Generate a simulated image of a strip 
# ========================================



# Store flattened data and labels
all_data = []
all_labels = []

n_samples = 1

for antigen in antigens:
    for _ in range(n_samples):
        strip = strip_colors(antigen, add_noise=False)
        all_data.append(strip)
        all_labels.append(antigen)

# Convert each strip to a 6D RGB vector: [D001_R, D001_G, D001_B, R007_R, R007_G, R007_B]
X = []
for strip in all_data:
    D001_rgb = strip.loc["D001", ["R", "G", "B"]].values
    R007_rgb = strip.loc["R007", ["R", "G", "B"]].values
    X.append(np.concatenate([D001_rgb, R007_rgb]))
X = np.array(X)

print(X)
