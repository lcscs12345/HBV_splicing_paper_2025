#!/usr/bin/env Rscript

# 1. Get the script and repo directories

# Get all command-line arguments, including the script file itself
cmd_args <- commandArgs(trailingOnly = FALSE)

# Find the argument that specifies the script file
# This argument usually starts with "--file="
script_arg_idx <- grep("--file=", cmd_args)

if (length(script_arg_idx) > 0) {
  # Extract the path after "--file="
  script_path <- sub("--file=", "", cmd_args[script_arg_idx[1]])
} else {
  # Fallback for scenarios where --file= is not present (e.g., sourcing in R console)
  script_path <- Sys.getenv("R_SCRIPT_SOURCE") # RStudio specific env var
  if (script_path == "" && interactive()) {
    message("Cannot reliably determine script path in this interactive session.")
    script_path <- getwd() # Fallback to current working directory, might not be accurate unless running from the script's directory 
  } else if (script_path == "") {
    stop("Could not determine script path in non-interactive mode. Please run with Rscript.")
  }
}

# Normalize the path (resolves '..', '.', '~', etc.) and get the directory of the script
script_dir <- dirname(normalizePath(script_path))
repo_dir <- normalizePath(file.path(script_dir, ".."))


# args <- commandArgs(trailingOnly = TRUE)

# # Check if arguments were provided
# if (length(args) == 0) {
#   stop("No arguments provided. Usage: Rscript track.r <repo_dir>", call. = FALSE)
# }

# repo_dir <- args[1]
intermediate_dir <- file.path(repo_dir, "results", "intermediate_files")

# -- Ensure the directory exists before trying to set it as WD --
if (!dir.exists(intermediate_dir)) {
  message("Creating directory: ", intermediate_dir)
  dir.create(intermediate_dir, recursive = TRUE)
}

# --- Set the working directory ---
setwd(intermediate_dir)
message("Current working directory set to: ", getwd())


# 2. Install required packages if not already installed

# --- Configuration ---
# List of Bioconductor packages
required_bioc_packages <- c(
  "GenomicFeatures",
  "Gviz",
  "trackViewer"
)

# Define the exact R major.minor version and Bioconductor *release number* to be used
enforced_r_major <- "4"
enforced_r_minor <- "5"

target_bioc_version <- "3.22"

# --- Function to check and install packages ---
install_bioc_packages <- function(packages, enforced_r_major, enforced_r_minor, target_bioc_version) {

  # R version check
  current_r_major <- R.version$major
  current_r_minor <- sub("^(\\d+).*$", "\\1", R.version$minor) # Get only the minor part

  if (current_r_major != enforced_r_major || current_r_minor != enforced_r_minor) {
    stop(
      "\nFATAL ERROR: Your current R version (", R.version$major, ".", R.version$minor, ") is not compatible.\n",
      "This script REQUIRES R version ", enforced_r_major, ".", enforced_r_minor, ".x (e.g., R ", enforced_r_major, ".", enforced_r_minor, ".0).\n",
      "This R version is necessary for Bioconductor release ", target_bioc_version, ".\n",
      "Please ensure R ", enforced_r_major, ".", enforced_r_minor, " is installed."
    )
  } else {
    message("R version (", R.version$major, ".", R.version$minor, ") is compatible with required R ", enforced_r_major, ".", enforced_r_minor, ".x for Bioconductor ", target_bioc_version, ".")
  }

  # Install BiocManager if not already installed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    message("BiocManager not found. Installing BiocManager...")
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      stop("Failed to install BiocManager. Please check your internet connection and R permissions.")
    }
  }

  # Check and install required Bioconductor packages
  message("\nChecking for required Bioconductor packages for Bioconductor version '", target_bioc_version, "'...")

  installed_pkgs <- installed.packages()[, "Package"]
  packages_to_install <- character(0)

  for (pkg in packages) {
    if (!(pkg %in% installed_pkgs)) {
      packages_to_install <- c(packages_to_install, pkg)
    } else {
      message("Package '", pkg, "' is already installed.")
    }
  }

  if (length(packages_to_install) > 0) {
    message("\nAttempting to install the following Bioconductor packages (version = '", target_bioc_version, "'):")
    print(packages_to_install)

    BiocManager::install(packages_to_install, version = target_bioc_version, ask = FALSE)

    # Verify installation
    for (pkg in packages_to_install) {
      if (!requireNamespace(pkg, quietly = TRUE)) {
        warning("Failed to install '", pkg, "'. Please check the installation logs for errors.")
      } else {
        message("Successfully installed '", pkg, "'.")
      }
    }
  } else {
    message("All required Bioconductor packages are already installed.")
  }

  # Validate the overall Bioconductor environment
  message("\nValidating Bioconductor environment...")
  # This will now validate against the packages installed for target_bioc_version
  valid_check <- BiocManager::valid()
  if (isTRUE(valid_check)) {
    message("Bioconductor environment is valid and all packages are up-to-date for Bioconductor version ", BiocManager::version(), ".")
  } else {
    message("Bioconductor environment has issues. Details from BiocManager::valid():")
    print(valid_check)
    message("Consider running BiocManager::install(version = '", target_bioc_version, "') without arguments to update all out-of-date packages.")
  }
}

# Execute the package installation function
message("Starting Bioconductor package setup...")
install_bioc_packages(required_bioc_packages, enforced_r_major, enforced_r_minor, target_bioc_version)
message("\nBioconductor package setup complete.")


# 3. Plot Top 10 HBV transcripts with coSI
library(trackViewer)
library(Gviz)
library(GenomicFeatures)
options(ucscChromosomeNames=FALSE)

tcons <- txdbmaker::makeTxDbFromGFF("track.gtf")
tcons <- GeneRegionTrack(tcons,transcriptAnnotation="symbol",showId=TRUE,fill="#FF0000")

fig1_track <- file.path(repo_dir, "results", "figures", "fig1", "track.pdf")
pdf(fig1_track, width = 6, height = 1.5)
plotTracks(tcons, chromosome = "HBV")
dev.off()

pos <- read.csv('track.lol.txt', sep='\t')
# Perform min max scaling
scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
pos$cex <- scale_values(pos$coSI)
pos$color <- ifelse(pos$ss == "5ss", "lightcyan", "lavender")
pos$border <- ifelse(pos$ss == "5ss", "#008080", "#4B0082")
gr <- GRanges("HBV genome", IRanges(pos$pos, width=1, names=pos$label, ss=pos$ss, color=pos$color, cex=pos$cex, border=pos$border))
genome <- GRanges("HBV genome", IRanges(c(1, 3246)))

fig1_lol <- file.path(repo_dir, "results", "figures", "fig1", "lol.pdf")
pdf(fig1_lol,width = 7.5, height = 2)
lolliplot(gr, genome)
dev.off()


