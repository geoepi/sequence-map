library(sequencemap)

# Set these paths for your analysis inputs.
fasta_path <- "path/to/alignment.fasta"
metadata_path <- "path/to/metadata.csv"
boundary_path <- "path/to/boundary.gpkg"
output_dir <- "path/to/output-directory"

# Core analysis settings.
sequence_id_col <- "sequence_id"
lon_col <- "longitude"
lat_col <- "latitude"
location_col <- NULL
group_col <- NULL
crs_input <- 4326
crs_projected <- 5070
grid_resolution <- 10000
n_pca <- 20
n_axes_model <- 3
drop_invariant <- TRUE
continue_on_error <- TRUE
overwrite <- FALSE
verbose <- TRUE

# Default workflow: PCA-first, no DAPC.
run_dapc <- FALSE
infer_dapc_groups <- FALSE

# Optional DAPC example:
# run_dapc <- TRUE
# group_col <- "your_group_column"
# infer_dapc_groups <- FALSE

workflow_result <- run_sequence_map_workflow(
  fasta_path = fasta_path,
  metadata_path = metadata_path,
  boundary_path = boundary_path,
  output_dir = output_dir,
  sequence_id_col = sequence_id_col,
  lon_col = lon_col,
  lat_col = lat_col,
  location_col = location_col,
  group_col = group_col,
  crs_input = crs_input,
  crs_projected = crs_projected,
  grid_resolution = grid_resolution,
  n_pca = n_pca,
  n_axes_model = n_axes_model,
  run_dapc = run_dapc,
  infer_dapc_groups = infer_dapc_groups,
  drop_invariant = drop_invariant,
  continue_on_error = continue_on_error,
  overwrite = overwrite,
  verbose = verbose
)

print(workflow_result$output_dir)
