% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluation_metrics.R
\name{plotContTable}
\alias{plotContTable}
\title{Plot the confution matrix}
\usage{
plotContTable(
  est_label,
  true_label,
  short.names = NULL,
  y.ord = NULL,
  x.ord = NULL,
  xlab = "Reference",
  ylab = "",
  threshold = NULL
)
}
\arguments{
\item{est_label}{vector of estimated labels}

\item{true_label}{vector of true labels}

\item{short.names}{(optional) rename the reference true labels for plots}

\item{y.ord}{(optional) a vector lenght of the number of rows of table, re-ordering the rows}

\item{x.ord}{(optional) a vector lenght of the number of columns of table, re-ordering the columns}

\item{xlab}{(optional) string, x axis label (default reference)}

\item{ylab}{(optional) string, y axis label (default empty)}

\item{threshold}{(optional) non-negative value, lower threshold, above which the values are shown (default NULL)}
}
\value{
the plot (ggplot object)
}
\description{
Plot the confution matrix
}
