#' Plot Multiple Migrations
#' @description  This function is used to plot all the migration events  based on the output of \code{PlotMS} with \code{plot.out} = FALSE and \code{demo.out} = TRUE. Plot settings should be customized  in \code{\link{PlotMS}}. Use function \code{\link{PlotMig}} to customize the plot of single migration.
#' @param demograph_out A list of all demographic information. See the return value description of \code{\link{PlotMS}}.
#' @param mig_par A list of all settings for plotting the demographic graph. See the return value description of \code{\link{PlotMS}}.
#' @param m.adjust Migration threshold for plotting migration events. Migration events with strength higher than \code{m.adjust} will be shown. The migration strength is defined as the proportion of the target population being replaced by the source population per generation. Default value is 0.
#' @seealso \code{\link{PlotMS}}
#' @examples
#' mig.cmd <- "./ms 15 100 -t 3.0 -I 6 0 7 0 0 8 0 -m 1 2 2.5 -m 2 1 2.5 -m 2 3 2.5 
#' -m 3 2 2.5 -m 4 5 2.5 -m 5 4 2.5 -m 5 6 2.5 -m 6 5 2.5 -em 2.0 3 4 2.5 
#' -em 2.0 4 3 2.5"
#' out<-PlotMS(input.cmd = mig.cmd, type = "ms", N4 = 10000, 
#' plot.out = FALSE, demo.out = TRUE, col.pop = 1:6, cex.lab = 1.5);
#' PlotMMig(demograph_out = out$demograph_out, mig_par = out$mig_par)
#' @export
PlotMMig <- function(demograph_out, mig_par, m.adjust = 0)
{
	time.scale <- mig_par$time.scale;
	time <- RescaleT(demograph_out, time.scale = mig_par$time.scale);
#count the number of events to plot
	par(mfrow = mig_par$mfrow);
	for(i in mig_par$events){
		x<-PlotMig(event = i, mig_par = mig_par, demograph_out = demograph_out, m.adjust = m.adjust, size.scale = "topology");
		if(i < mig_par$events[length(mig_par$events)]){
			title(paste(time[i]," to ", time[i + 1]," (", time.scale, ")", sep = ""));
		} else {
			title(paste("Time : >= ", time[i], "(", time.scale, ")", sep = ""), cex = mig_par$cex.lab * 1.2)
		}
		legend("topleft", legend = mig_par$lab.pop, col = mig_par$col.pop, pch = 20, bty = "n", cex = mig_par$cex.lab);
	}
}

#' Plot migration event(s) at a particular time
#' @description  This function plots migration events at a particular time point based on the output of \code{PlotMS} with \code{demo.out} = T and \code{plot.out} = F. The \code{time_pt} or \code{event} determines the time points that are plotted. The \code{add} and \code{map.pos} parameters allow the migration graph to be added to other backgrounds such as maps. 
#' @param time_pt A numeric value defining the time point for the migrations which will be plotted. \code{time_pt} should be in the scale defined by \code{time.scale}. For example, \code{time_pt} = 3 when \code{time.scale} = "log10year" corresponds to the migrations 10^3 years ago.
#' @param time.scale A keyword to define the time scale used in the plot. It can be "4Ne", "generation", "year", "kyear", and "log10year". When the \code{size.scale} = "topology", this parameter will be ignored.
#' @param event An index to define at which time to plot migration(s). Every demographic event has an index in the order of time. Demographic changes at the same time are all defined as the same event and share the same index. 
#' @param demograph_out A list containing all demographic information, see \code{\link{PlotMS}}.
#' @param mig_par A list contained all settings for plotting the demographic graph, see \code{\link{PlotMS}}.
#' @param size.scale A keyword to define the scaling of lineage width. "topology" returns only topology structure among simulated populations, ignoring both the population sizes and the length of time between any demographic events. "linear" linearly scales the lineage widths as a function of the population size, with the scale factor defined by the variable \code{linear.scale}; "log" scales the lineage width logarithmically as a function of the population size. The logarithm base is defined by the variable \code{log.base}.
#' @param m.adjust Migration threshold for plotting migration events. Migration events with strength higher than \code{m.adjust} will be shown. The migration strength is defined as the proportion of the target population being replaced by the source population per generation. Default value is 0.
#' @param toposize.scale Control the scaling of the size of circle when the \code{size.scale} = "topology".
#' @param add A logical value allowing one to add the migrations to another background (2-D only). Positions for every population dime should be defined in \code{map.pos}.
#' @param map.pos A matrix with two columns, the ith row is the coordinate for the ith population.
#' @param linear.scale Linear scale magnitude, to be applied when \code{size.scale} = "linear".
#' @param log.base The base of logarithm, to be applied when \code{size.scale} = "log".
#' @param col.pop Population lineage color.
#' @param col.arrow Migration arrow color.
#' @param lwd.arrow Control the line width of arrow representing a migration. The arrow width is defined by 0.5+migration strength*\code{lwd.arrow}.
#' @param length.arrowtip Length of the edges of the arrow tip.
#' @param angle.arrowtip The angle of the arrow tip, between 0 and 90.
#' @param xlim The range of x-axis.
#' @param ylim The range of y-axis.
#' @seealso \code{\link{PlotMS}}
#' @examples 
#' test.mig.cmd <- "./ms 15 100 -t 3.0 -I 6 0 7 0 0 8 0 -m 1 2 2.5 -m 2 1 2.5 
#' -m 2 3 2.5 -m 3 2 2.5 -m 4 5 2.5 -m 5 4 2.5 -m 5 6 2.5 -m 6 5 2.5 
#' -em 2.0 3 4 2.5 -em 2.0 4 3 2.5"
#' out <- PlotMS(input.cmd = test.mig.cmd, type = "ms", 
#' N4 = 10000, plot.out = FALSE, demo.out = TRUE);
#' #check all migration events
#' events <- out$mig_par$events
#' print(events)
#' #check the time for those migration events
#' timeofevents <- out$mig_par$time[events]
#' print(timeofevents)
#' #plot event by event
#' par(mfrow = c(1, 2))
#' PlotMig(event = 1, demograph_out = out$demograph_out, mig_par = out$mig_par)
#' title("Event-1");
#' PlotMig(event = 2, demograph_out = out$demograph_out, mig_par = out$mig_par, 
#' col.pop = 1:6, xlim = c(-5,4))
#' title("Event-2", cex.main = 3);
#' legend("topleft", col = 1:6, pch = 20, bty = "n", cex = 2,
#' legend = c("pop-1", "pop-2", "pop-3", "pop-4", "pop-5", "pop-6"))
#' @export
PlotMig<-function(time_pt = NULL, event = 1, mig_par, demograph_out,
		size.scale = mig_par$size.scale, time.scale = mig_par$time.scale,
		linear.scale = mig_par$linear.scale, log.base = mig_par$log.base,
		col.pop = mig_par$col.pop, col.arrow = mig_par$col.arrow,
		xlim = mig_par$xlim, ylim = mig_par$ylim,
		lwd.arrow = mig_par$lwd.arrow, 
		length.arrowtip = mig_par$length.arrowtip, 
		angle.arrowtip = mig_par$angle.arrowtip,
		toposize.scale = 1, add = FALSE, map.pos = NULL, m.adjust = 0)
{
	if(!is.null(time_pt)){
		if(time.scale!= mig_par$time.scale){
			#print(c(time.scale, mig_par$time.scale))
			time_pt <- TransT(time_pt, scale1 = time.scale, scale2 = mig_par$time.scale, demograph_out);
		}
		for(i in 2:length(mig_par$time)){
#event_pt is used to record the relationship between the real event and the assigned event
			if(time_pt >= mig_par$time[i-1] & time_pt < mig_par$time[i]){
				event_pt = (time_pt - mig_par$time[i-1]) / (mig_par$time[i] - mig_par$time[i - 1]) + i - 1;
				event = i - 1;
				break;
			}
			else {event_pt = i; event = i;}
		}
	} 
	else {
		event_pt = event;
		time_pt = mig_par$time[event];
	}

	m <- demograph_out$m[ , , event];
#adjust the migration matrix
	m[m < m.adjust] = 0;
	total.pop.num = demograph_out$total.pop.num;
	N = demograph_out$N[event, ];

#m is a two-dimension matrix for the migration rates at some time point
	if(add == FALSE){
		plot(NA, NA, xlim = xlim, ylim = ylim, axes = F, xpd = TRUE, xlab = "", ylab = "");
	}
	plot_pos <- c();
	for(i in 1:total.pop.num){
		x <- PlotCircle(ind = i, time_pt = time_pt, event_pt = event_pt, mig_par = mig_par, demograph_out = demograph_out, col = col.pop, size.scale = size.scale, linear.scale = linear.scale, log.base = log.base, toposize.scale = toposize.scale, add = add, map.pos = map.pos, xlim = xlim, ylim = ylim);
		plot_pos <- rbind(plot_pos, x);
	}
	arrow_pos <- c();
	for(ind1 in 1:total.pop.num){
		for(ind2 in 1:total.pop.num){
##from ind2 to ind1
#if(ind1!= ind2&m[ind1, ind2]>0&m[ind1, ind2]<1){
	if(ind1 != ind2 & m[ind1, ind2] > 0 & event_pt <= demograph_out$survive[2, ind1] & event_pt >= demograph_out$survive[1, ind1]){
		x1 = plot_pos[ind2,1];y1 = plot_pos[ind2,2];r1 = plot_pos[ind2,3];
		x2 = plot_pos[ind1,1];y2 = plot_pos[ind1,2];r2 = plot_pos[ind1,3];
		x <- CalArrowPosCircle(x1 = x1, y1 = y1, r1 = r1, x2 = x2, y2 = y2, r2 = r2);
		arrows(x[1], x[2], x[3], x[4], lwd = lwd.arrow, length = length.arrowtip, angle = angle.arrowtip, col = col.arrow[ind2]);
		arrow_pos <- rbind(arrow_pos,x);
	}
}
}
}

#' Output population sizes at a particular time
#' @description Output population sizes at a particular time.
#' @param time_pt A numeric value defining the time point for the migrations which will be plotted. \code{time_pt} should be in the scale defined by \code{time.scale}. For example, \code{time_pt} = 3 when \code{time.scale} = "log10year" corresponds to the migrations 10^3 years ago.
#' @param time.scale A keyword to define the time scale used in the plot. It can be "4Ne", "generation", "year", "kyear", and "log10year". When the \code{size.scale} = "topology", this parameter will be ignored.
#' @param demograph_out A list containing all demographic information, see \code{\link{PlotMS}}.
#' @return A vector of all population sizes for the specified time.
#' @export
NOut<-function(time_pt, time.scale, demograph_out)
{
	times <- demograph_out$time.series;
	if(time.scale!= "4Ne"){
		#print(c(time.scale, mig_par$time.scale))
		time_pt <- TransT(time_pt, scale1 = time.scale, scale2 = "4Ne", demograph_out);
	}
	for(i in 2:length(times)){
		if(time_pt >= times[i-1] & time_pt < times[i]){
			gap <- time_pt - times[i-1];
			event = i - 1;
			break;
		}
		else {gap = 0; event = i;}
	}
	N <- demograph_out$N[event, ]*demograph_out$N4;
	g.rate <- demograph_out$g.rate[event, ];
	N.out <- N * exp(-gap * g.rate);
	N.out[i < demograph_out$survive[1, ] | i >= demograph_out$survive[2, ]] <- NA;
	#print(event)
	floor(N.out)
}

#' Add population size ruler
#' @description Adds a ruler for the population size. This function works like the function \code{legend}, 
#' @param x,y Position of the population size ruler. If y does not have numeric value, x will support keyword input from the list "bottomright", "bottom", "bottomleft","left", "topleft", "top", "topright", "right" and "center".
#' @param Nsize The population sizes of the ticks on the ruler.
#' @param Nlab The labels on the ticks of the ruler. Default labels are the population index.
#' @param N4 Scalar to scale population size.
#' @param size.scale A keyword to define the way population size is scaled and displayed. It supports "log" and "linear".
#' @param linear.scale The scale factor applied to the population size when \code{size.scale} = "linear".
#' @param log.base The base of logarithm used when \code{size.scale} = "log".
#' @param ... Additional arguments can be passed, such as \code{col}, \code{lwd}, \code{lty}.
#' @export

NRuler <- function(x, y = NULL, Nsize, Nlab = Nsize, N4, size.scale, linear.scale = 0.2, log.base = 10, ...){
	if(size.scale == "topology"){
		stop("No size difference among populations in topology plot");
	}
	if(size.scale == "linear"){
		if(linear.scale == 0.2){
			warning("linear.scale is in the default value: 0.2")
		}
		seglen = Nsize / N4 * linear.scale;
	}
	if(size.scale == "log"){
		if(log.base == 10){
			warning("log.base is in the default value: 10")
		}
		seglen = log(Nsize / N4 / 2 + 1, log.base) * 2;
	}
	#print(seglen)
	Nmax <- max(seglen);
	x <- legend(x = x, y = y, legend = c(""), bty = "n")$rect;
	#print(x);
	left <- x$left + x$w/5; right <- left + Nmax; height <- x$top -x$h/2;
	lines(c(left, right), c(height, height), ...);
	lines(c(left,left), c(height, height - x$h / 10), ...);
	text(left, height - x$h / 4, label = "0", ...);
	for(i in 1:length(seglen)){
		lines(c(left,left) + seglen[i], c(height, height + x$h/10), ...);
		text(left + seglen[i], height + x$h/4, label = Nlab[i], ...);
	}
}

#' Plot population demographic graph and generate demographic parameters
#' @description  This is the main function to plot demographic graph for single/multiple populations. The function is named after Hudson's \code{ms} program. It can read the simulation input data used for the \code{ms}, \code{msa}, \code{msHot}, \code{MaCS}, \code{scrm}, and \code{cosi} programs.
#' @description The \code{input.file} or \code{imput.cmd} and command \code{type} are required to plot the demographic history. The output graph can be customized.
#' @description In the demographic graph, each population has a lineage that stretches back in time. The width of the lineage reflects the population size. Population splits and migrations are represented by arrows. 
#' @param input.cmd An input string containing the simulation program input commands.
#' @param input.file A file containing the simulation program input commands (for \code{ms}, \code{msa}, or \code{MaCS}), or parameter files (for \code{cosi}).
#' @param type A keyword indicating the type of simulation command: "ms" for \code{ms}, "msa" for \code{msa}, "macs" for \code{MaCS}, "scrm" for \code{scrm}, "cosi" for \code{cosi}, "msprime" for msprime. Please check the \href{https://github.com/YingZhou001/POPdemog/blob/master/doc/POPdemog-tutorial.md}{online tutorial file} to see more supported simulation programs.
#' @param inpos Population positions in the plot at time 0. 
#' @param N4 Four times the effective population size. This parameter has the same definition as the 4N0 parameter for the \code{ms} simulation program. 
#' @param plot.out A logical variable that controls the production of the demographic plot. If TRUE, the demographic plot will be produced.
#' @param demo.out A logical variable that controls the output of the demographic parameters. If TRUE, all demographic parameters that are used for the graph will be returned.
#' @param size.scale A keyword to define the scaling of lineage width. "topology" returns only topology structure among simulated populations, ignoring both the population sizes and the length of time between any demographic events. "linear" linearly scales the lineage widths as a function of the population size, with the scale factor defined by the variable \code{linear.scale}; "log" scales the lineage width logarithmically as a function of the population size. The logarithm base is defined by the variable \code{log.base}.
#' @param linear.scale Linear scale factor,which will be applied when \code{size.scale} = "linear".
#' @param log.base The base of logarithm, which will be applied when \code{size.scale} = "log".
#' @param time.scale A keyword to define the time scale used in the plot. It can be "4Ne", "generation", "year", "kyear", and "log10year". When the \code{size.scale} = "topology", this parameter will be ignored.
#' @param gen Years per generation. Default value is 25.
#' @param m.adjust Migration threshold for plotting migration events. Migration events with strength higher than \code{m.adjust} will be shown. The migration strength is defined as the proportion of the target population being replaced by the source population per generation. Default value is 0. This value should be between 0 and 1.
#' @param col.pop Color for each population.
#' @param col.arrow Color for each migration arrow.
#' @param length.arrowtip Size of arrow tips.
#' @param lwd.arrow Controls the width of arrow representing a migration. The arrow width is determined by 0.5+migration strength*\code{lwd.arrow}.
#' @param angle.arrowtip Arrow end angle, between 0 and 90.
#' @param pops Population name labels. Default as 1:number of populations.
#' @param xlim Range of x-axis.
#' @param ylim Range of y-axis.
#' @param xlab Title for the x-axis.
#' @param ylab Title for the y-axis.
#' @param cex.lab The magnification to be used for x and y labels relative to the current setting of \code{cex}.
#' @param cex.axis The magnification to be used for axis annotation relative to the current setting of \code{cex}.
#' @param axes A logical value to plot the axes or not.
#' @return if the parameter \code{plot} = F/FALSE, the following three lists will be returned:
#' @return \item{demograph_out}{This list contains all demographic details from the input command file:
#' \itemize{
#'	\item \code{time.series} is a vector of time;
#'	\item \code{Pos} is a numeric matrix of positions for each population at every demographic event;
#'	\item \code{N} is a numeric matrix of population size for each population at every demographic event;
#'	\item \code{m} is a 3-D numeric matrix of migration rates between populations at every demographic event;
#'	\item \code{survive} is a matrix recording the begining and end for each population according to the demographic events;
#'	\item \code{g.rate} is a matrix of exponential growth rates at every demographic event;
#'	\item \code{pop.pos} is a numeric vector of the population positions at time 0;
#'	\item \code{pop.lab} is a vector of population names;
#'	\item \code{mscmd} is the \code{ms} command for the demographic plot. Demographic information from simulation scripts will be turned to ms command format for further extraction;
#'	\item \code{present.pop.num} is the number of populations at present;
#'	\item \code{total.pop.num} is the number of total populations exist in the plot;
#'	\item \code{N4} is 4Ne;
#'	\item \code{gen} is the number of years per generation.
#'      }
#' }
#' @return \item{evo_par}{This list contains all parameters used to draw the demographic graph, including:
#' 	\code{size.scale},\code{linear.scale}, \code{log.base}, \code{time.scale}, \code{time}, \code{col.pop}, \code{col.arrow}, \code{length.arrowtip}, \code{lwd.arrow}, \code{angle.arrowtip}, \code{lab.pop}, \code{lab.pos}, \code{xlim}, \code{ylim}, \code{xlab}, \code{ylab}, \code{cex.lab}, \code{cex.axis}, \code{axes}.
#' See more details in the parameter description.
#' }
#' @return \item{mig_par}{This list contains all parameters used to draw the migrations, including:
#' 	\code{size.scale}, \code{linear.scale}, \code{log.base}, \code{time.scale}, \code{time}, \code{lab.pop}, \code{col.pop}, \code{col.arrow}, \code{length.arrowtip}, \code{lwd.arrow}, \code{angle.arrowtip}, \code{xlim},\code{ylim}, \code{events}, \code{cex.lab}
#' See more details in the parameter description.
#'}
#' @references 4Ne: \url{http://home.uchicago.edu/rhudson1/source/mksamples.html}
#' @examples 
#' #example 1
#' test.1.ms.cmd <- "./ms 44 1 -r 20000 50000000 -t 30000 
#' -I 6 20 20 1 1 1 1 -en 0 1 1 
#' -en 0 2 1 -en 0 3 1e-10 -en 0 4 1e-10 -en 0 5 1e-10 -en 0 6 1e-10 
#' -es 0.0125 2 0.97 -en 0.02500025 7 0.25 -en 0.02500025 2 1 -ej 0.05 4 3 
#' -ej 0.05 6 5 -en 0.05000025 3 0.25 -en 0.05000025 5 0.25 -ej 0.0500025 5 3 
#' -en 0.050005 3 0.25 -ej 0.075 2 1 -en 0.0750025 1 1 -ej 0.1 7 3 
#' -en 0.1000025 3 0.25 -ej 0.3 3 1 -en 0.3000025 1 1"
#' PlotMS(input.cmd = test.1.ms.cmd, type = "ms", N4 = 10000, time.scale = "kyear")
#' #adjust the population position
#' PlotMS(input.cmd = test.1.ms.cmd, type = "ms", N4 = 10000, time.scale = "kyear", 
#' inpos = c(1,2,5,4,6,7,3))
#' #add color for each population
#' PlotMS(input.cmd = test.1.ms.cmd, type = "ms", N4 = 10000, time.scale = "kyear", 
#' inpos = c(1,2,5,4,6,7,3), col.pop = rainbow(10)[3:9])
#' #add population names
#' 
#' #example 2
#' "./ms 15 100 -t 3.0 -I 6 0 7 0 0 8 0 -m 1 2 2.5 -m 2 1 2.5 -m 2 3 2.5 
#' -m 3 2 2.5 -m 4 5 2.5 -m 5 4 2.5 -m 5 6 2.5 -m 6 5 2.5 -em 2.0 3 4 2.5 
#' -em 2.0 4 3 2.5" -> test.2.ms.cmd
#' PlotMS(input.cmd = test.2.ms.cmd, type = "ms", N4 = 10000, col.pop = 3, 
#' col.arrow = 1, lwd.arrow = 1.5)
#' 
#' #example 3
#' "./ms 1 1 -t 1.0 -I 3 10 10 10 -n 1 1.682020 -n 2 3.736830 
#' -n 3 7.292050 -eg 0 2 116.010723 -eg 0 3 160.246047 
#' -ma 0 0.881098 0.561966 0.881098 0 2.797460 0.561966 2.797460 0 
#' -ej 0.028985 3 2 -en 0.028985 2 0.287184 
#' -ema 0.028985 3 0 7.293140 0 7.293140 0 0 0 0 0 -ej 0.197963 2 1 
#' -en 0.303501 1 1" -> Ryan2009.cmd
#' PlotMS(input.cmd = Ryan2009.cmd, type = "ms", N4 = 10000, size.scale = "log", 
#' log.base = 200, pops = c("AFR", "EUR", "ESA"), 
#' col.pop = c("brown", "blue", "yellow"));
#' 
#' #example 4
#' ms 4 1 -t 7156.0000000 -r 2000.0000 10000000 -eN 0 5 
#' -eG 0.000582262 1318.18 -eG 0.00232905 -329.546 -eG 0.00931619 82.3865 
#' -eG 0.0372648 -20.5966 -eG 0.149059 5.14916 -eN 0.596236 0.5 -T" -> zigzag.cmd
#' par(mfrow = c(1,2))
#' PlotMS(input.cmd = zigzag.cmd, type = "ms", N4 = 10000)
#' #change the time scale
#' PlotMS(input.cmd = zigzag.cmd, type = "ms", N4 = 10000, time.scale = "log10year")
#' @export
PlotMS<-function(input.cmd = NULL, input.file = NULL, 
		 type, inpos = NULL, N4 = 10000,
		 size.scale = "linear", linear.scale = 0.2,log.base = 10,
		 time.scale = "4Ne", gen = 25, m.adjust = 0,
		 col.pop = "gray45", col.arrow = col.pop,
		 length.arrowtip = 0.15, lwd.arrow = 1, angle.arrowtip = 15,
		 pops = NULL, xlab = "Population", 
		 ylab = paste("Time before present (", time.scale, ")", sep = ""),
		 xlim = NULL, ylim = NULL, 
		 plot.out = T, demo.out = F,
		 cex.lab = 1, cex.axis = 1, axes = T)
{
	tempfile <- file.path(tempdir(), "PlotMS.cmdtmp");
	if(!is.null(input.cmd)&!is.null(input.file)){stop("Too many inputs!!!")}
	else if(!is.null(input.cmd)){
		cat(input.cmd, file = tempfile);
		input.file = tempfile;
	}

	demograph_out <- ReadDemoGraph(inputfile = input.file, type = type, inpos = inpos, N4 = N4, gen = gen, pops = pops);
	if(!is.null(input.cmd)){unlink(tempfile);}
	message("demographic initiation, done!");
	#print(demograph_out)

	evo_par <- ReadEvoPar(demograph_out = demograph_out,
			      size.scale = size.scale, linear.scale = linear.scale,log.base = log.base,
			      time.scale = time.scale,
			      col.pop = col.pop, col.arrow = col.arrow,
			      length.arrowtip = length.arrowtip, lwd.arrow = lwd.arrow,
			      angle.arrowtip = angle.arrowtip,			
			      pops = demograph_out$pop.lab,
			      xlab = xlab, ylab = ylab,
			      xlim = xlim, ylim = ylim,		
			      cex.lab = cex.lab, cex.axis = cex.axis, axes = axes);

	mig_par <- ReadMigPar(demograph_out = demograph_out,
			      size.scale = size.scale, linear.scale = linear.scale,
			      log.base = log.base,
			      time.scale = time.scale,
			      pops = demograph_out$pop.lab, cex.lab = cex.lab,
			      col.pop = col.pop, col.arrow = col.arrow,
			      length.arrowtip = length.arrowtip, lwd.arrow = lwd.arrow,
			      angle.arrowtip = angle.arrowtip);

	message("plot initiation, done!");
	#########plot region
	#################################################Evolution relationship
	if (plot.out == T) {
		PlotEvo(demograph_out = demograph_out, evo_par = evo_par, m.adjust = m.adjust);
	}
	if (demo.out == T) {
		list(demograph_out = demograph_out, evo_par = evo_par, mig_par = mig_par);
	}
}
