#' Mouse somatosensory cortex seqFISH+ dataset
#'
#' A processed spatial dataset from a single seqFISH+ slide of the mouse somatosensory cortex.
#'
#' @format A list with the following structure:
#' \describe{
#'   \item{samples}{A named list of samples. Each sample contains a data.table with:
#'     \describe{
#'       \item{id}{Cell identifier.}
#'       \item{x}{x-coordinate of the cell.}
#'       \item{y}{y-coordinate of the cell.}
#'       \item{celltype}{Cell type annotation.}
#'       \item{region}{Region annotation.}
#'     }
#'   }
#'   \item{sample_id}{A character vector of sample identifiers.}
#' }
#' @source Eng, C.-H. L. et al. Transcriptome-scale super-resolved imaging in tissues by rna seqfish+. Nature 568, 235–239 (2019).460
"dataset_cortex"

#' Human melanoma imaging mass cytometry dataset
#'
#' A processed spatial dataset from human melanoma imaging mass cytometry (IMC).
#'
#' @format A list with the following structure:
#' \describe{
#'   \item{samples}{A named list of samples. Each sample contains a
#'     \code{data.table} with:
#'     \describe{
#'       \item{id}{Cell identifier.}
#'       \item{x}{x-coordinate of the cell.}
#'       \item{y}{y-coordinate of the cell.}
#'       \item{celltype}{Cell type annotation.}
#'       \item{region}{Region annotation.}
#'     }
#'   }
#'   \item{sample_id}{A character vector of sample identifiers.}
#'   \item{group}{A factor or character vector indicating immune infiltration degree.}
#'   \item{subject_id}{A character vector of subject identifiers.}
#'   \item{des}{A list containing a design matrix specifying variables used in regression analyses.}
#'   \item{coef}{A list containing regression coefficients corresponding to variables tested in the model.}
#' }
#'
#' @source Hoch, T. et al. Multiplexed imaging mass cytometry of the chemokine milieus in melanoma characterizes features of the477
#' response to immunotherapy. Sci. immunology 7, eabk1692 (2022).478
"dataset_melanoma"

#' Human type 1 diabetes imaging mass cytometry dataset
#'
#' A processed spatial dataset from human type 1 diabetes imaging mass cytometry (IMC).
#'
#' @format A list with the following structure:
#' \describe{
#'   \item{samples}{A named list of samples. Each sample contains a
#'     \code{data.table} with:
#'     \describe{
#'       \item{id}{Cell identifier.}
#'       \item{x}{x-coordinate of the cell.}
#'       \item{y}{y-coordinate of the cell.}
#'       \item{celltype}{Cell type annotation.}
#'       \item{region}{Region annotation.}
#'     }
#'   }
#'   \item{sample_id}{A character vector of sample identifiers.}
#'   \item{group}{A character vector indicating disease-related
#'     stratification.}
#'   \item{subject_id}{A character vector of subject identifiers.}
#'   \item{des}{A list containing a design matrix specifying variables used in regression analyses.}
#'   \item{coef}{A list containing regression coefficients corresponding to variables tested in the model.}
#' }
#'
#' @source Damond, N. et al. A map of human type 1 diabetes progression by imaging mass cytometry. Cell metabolism 29, 755–768526
#' (2019).527
"dataset_diabetes"




