#' @title Parallelized Elbow Method for Optimal K in K-means
#'
#' @description Runs k-means clustering in parallel for a range of cluster
#'   numbers (k) and returns the total within-cluster sum of squares (WCSS) for
#'   each k, along with an elbow plot to help determine the optimal number of
#'   clusters. Uses `furrr` for parallelization.
#'
#' @param data A numeric data frame or matrix to cluster.
#' @param k_range Integer vector. The range of cluster numbers to evaluate
#'   (default: 2:15).
#' @param nstart Integer. Number of random sets for k-means (default: 10).
#' @param iter.max Integer. Maximum number of iterations for k-means (default:
#'   100).
#' @param workers Integer. Number of parallel workers to use (default: one less
#'   than the number of available cores).
#'
#' @importFrom furrr future_map
#' @importFrom purrr map_dbl
#' @importFrom ggplot2 ggplot aes geom_line geom_point labs theme_minimal
#' @importFrom future plan multisession sequential
#'
#' @examples
#' \dontrun{
#' result <- find_optimal_k(your_data, k_range = 2:10)
#' print(result$plot)
#' }
#'
#' @return A list with two elements: \code{wcss} (numeric vector of WCSS for
#'   each k) and \code{plot} (a ggplot2 elbow plot).
#'
find_optimal_k <- function(
  data,
  k_range = 2:15,
  nstart = 10,
  iter.max = 100,
  workers = parallel::detectCores() - 1
) {

  # set up parallel processing
  future::plan(
    future::multisession,
    workers = workers
  )

  # define kmeans function for a given k
  kmeans_fun <- function(k) {
    stats::kmeans(
      data,
      centers = k,
      nstart = nstart,
      iter.max = iter.max
    )
  }

  # run kmeans in parallel for each k
  kmeans_results <- furrr::future_map(
    k_range,
    kmeans_fun,
    .progress = TRUE
  )

  # extract total within-cluster sum of squares for each k
  wcss <- purrr::map_dbl(kmeans_results, "tot.withinss")

  # plot the elbow plot
  elbow_plot <- ggplot2::ggplot(
    data.frame(k = k_range, wcss = wcss),
    ggplot2::aes(x = k, y = wcss)
  ) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::labs(
      x = "Number of clusters (k)",
      y = "Total within-cluster sum of squares",
      title = "Elbow Method (Parallel K-means)"
    ) +
    ggplot2::theme_minimal()

  # reset parallel plan
  future::plan(future::sequential)

  return(
    list(
      wcss = wcss,
      plot = elbow_plot
    )
  )

}