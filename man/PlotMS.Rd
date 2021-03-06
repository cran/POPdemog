% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/popdemog.R
\name{PlotMS}
\alias{PlotMS}
\title{Plot population demographic graph and generate demographic parameters}
\usage{
PlotMS(input.cmd = NULL, input.file = NULL, type, inpos = NULL,
  N4 = 10000, size.scale = "linear", linear.scale = 0.2, log.base = 10,
  time.scale = "4Ne", gen = 25, m.adjust = 0, col.pop = "gray45",
  col.arrow = col.pop, length.arrowtip = 0.15, lwd.arrow = 1,
  angle.arrowtip = 15, pops = NULL, xlab = "Population",
  ylab = paste("Time before present (", time.scale, ")", sep = ""),
  xlim = NULL, ylim = NULL, plot.out = T, demo.out = F, cex.lab = 1,
  cex.axis = 1, axes = T)
}
\arguments{
\item{input.cmd}{An input string containing the simulation program input commands.}

\item{input.file}{A file containing the simulation program input commands (for \code{ms}, \code{msa}, or \code{MaCS}), or parameter files (for \code{cosi}).}

\item{type}{A keyword indicating the type of simulation command: "ms" for \code{ms}, "msa" for \code{msa}, "macs" for \code{MaCS}, "scrm" for \code{scrm}, "cosi" for \code{cosi}, "msprime" for msprime. Please check the \href{https://github.com/YingZhou001/POPdemog/blob/master/doc/POPdemog-tutorial.md}{online tutorial file} to see more supported simulation programs.}

\item{inpos}{Population positions in the plot at time 0.}

\item{N4}{Four times the effective population size. This parameter has the same definition as the 4N0 parameter for the \code{ms} simulation program.}

\item{size.scale}{A keyword to define the scaling of lineage width. "topology" returns only topology structure among simulated populations, ignoring both the population sizes and the length of time between any demographic events. "linear" linearly scales the lineage widths as a function of the population size, with the scale factor defined by the variable \code{linear.scale}; "log" scales the lineage width logarithmically as a function of the population size. The logarithm base is defined by the variable \code{log.base}.}

\item{linear.scale}{Linear scale factor,which will be applied when \code{size.scale} = "linear".}

\item{log.base}{The base of logarithm, which will be applied when \code{size.scale} = "log".}

\item{time.scale}{A keyword to define the time scale used in the plot. It can be "4Ne", "generation", "year", "kyear", and "log10year". When the \code{size.scale} = "topology", this parameter will be ignored.}

\item{gen}{Years per generation. Default value is 25.}

\item{m.adjust}{Migration threshold for plotting migration events. Migration events with strength higher than \code{m.adjust} will be shown. The migration strength is defined as the proportion of the target population being replaced by the source population per generation. Default value is 0. This value should be between 0 and 1.}

\item{col.pop}{Color for each population.}

\item{col.arrow}{Color for each migration arrow.}

\item{length.arrowtip}{Size of arrow tips.}

\item{lwd.arrow}{Controls the width of arrow representing a migration. The arrow width is determined by 0.5+migration strength*\code{lwd.arrow}.}

\item{angle.arrowtip}{Arrow end angle, between 0 and 90.}

\item{pops}{Population name labels. Default as 1:number of populations.}

\item{xlab}{Title for the x-axis.}

\item{ylab}{Title for the y-axis.}

\item{xlim}{Range of x-axis.}

\item{ylim}{Range of y-axis.}

\item{plot.out}{A logical variable that controls the production of the demographic plot. If TRUE, the demographic plot will be produced.}

\item{demo.out}{A logical variable that controls the output of the demographic parameters. If TRUE, all demographic parameters that are used for the graph will be returned.}

\item{cex.lab}{The magnification to be used for x and y labels relative to the current setting of \code{cex}.}

\item{cex.axis}{The magnification to be used for axis annotation relative to the current setting of \code{cex}.}

\item{axes}{A logical value to plot the axes or not.}
}
\value{
if the parameter \code{plot} = F/FALSE, the following three lists will be returned:

\item{demograph_out}{This list contains all demographic details from the input command file:
\itemize{
\item \code{time.series} is a vector of time;
\item \code{Pos} is a numeric matrix of positions for each population at every demographic event;
\item \code{N} is a numeric matrix of population size for each population at every demographic event;
\item \code{m} is a 3-D numeric matrix of migration rates between populations at every demographic event;
\item \code{survive} is a matrix recording the begining and end for each population according to the demographic events;
\item \code{g.rate} is a matrix of exponential growth rates at every demographic event;
\item \code{pop.pos} is a numeric vector of the population positions at time 0;
\item \code{pop.lab} is a vector of population names;
\item \code{mscmd} is the \code{ms} command for the demographic plot. Demographic information from simulation scripts will be turned to ms command format for further extraction;
\item \code{present.pop.num} is the number of populations at present;
\item \code{total.pop.num} is the number of total populations exist in the plot;
\item \code{N4} is 4Ne;
\item \code{gen} is the number of years per generation.
     }
}

\item{evo_par}{This list contains all parameters used to draw the demographic graph, including:
	\code{size.scale},\code{linear.scale}, \code{log.base}, \code{time.scale}, \code{time}, \code{col.pop}, \code{col.arrow}, \code{length.arrowtip}, \code{lwd.arrow}, \code{angle.arrowtip}, \code{lab.pop}, \code{lab.pos}, \code{xlim}, \code{ylim}, \code{xlab}, \code{ylab}, \code{cex.lab}, \code{cex.axis}, \code{axes}.
See more details in the parameter description.
}

\item{mig_par}{This list contains all parameters used to draw the migrations, including:
	\code{size.scale}, \code{linear.scale}, \code{log.base}, \code{time.scale}, \code{time}, \code{lab.pop}, \code{col.pop}, \code{col.arrow}, \code{length.arrowtip}, \code{lwd.arrow}, \code{angle.arrowtip}, \code{xlim},\code{ylim}, \code{events}, \code{cex.lab}
See more details in the parameter description.
}
}
\description{
This is the main function to plot demographic graph for single/multiple populations. The function is named after Hudson's \code{ms} program. It can read the simulation input data used for the \code{ms}, \code{msa}, \code{msHot}, \code{MaCS}, \code{scrm}, and \code{cosi} programs.

The \code{input.file} or \code{imput.cmd} and command \code{type} are required to plot the demographic history. The output graph can be customized.

In the demographic graph, each population has a lineage that stretches back in time. The width of the lineage reflects the population size. Population splits and migrations are represented by arrows.
}
\references{
4Ne: \url{http://home.uchicago.edu/rhudson1/source/mksamples.html}
}
