ReadEvoPar<-function(demograph_out,
		size.scale="linear", linear.scale=0.2,log.base=10,
		time.scale="4Ne",
		col.pop="gray45", col.arrow=col.pop,
		length.arrowtip=0.15, lwd.arrow=1,
		angle.arrowtip=15,pops=demograph_out$pop.lab,
		xlim=NULL, ylim=NULL,
		xlab="Population", ylab=paste("Time before present (",time.scale,")", sep=""),
		cex.lab=1, cex.axis=1, axes=T)
{
## three plot modes can be used to adjust the presentation of population size in the plot, the width of each population lineage. size.scale="linear" gives the proportion of linear.scale to original width; size.scale="log" gives the logarithms of base log.base of original width, and size.scale="topology" will treat the width to be zero.
## time.scale gives the unit in the y-axis, has the option of "year", "generation", or with default value as the ms' 4Ne scaled time.

	if(size.scale!="linear"&&size.scale!="log"&&size.scale!="topology")print("Error, you need to redefine the parameter size.scale");
#re-scale time
		time<-RescaleT(demograph_out, time.scale=time.scale);
		if(time.scale!="4Ne"&time.scale!="generation"&time.scale!="year"&time.scale!="kyear"&time.scale!="log10year"){print("Error, you need to redefine the time.scale");}
#color vector
		col.pop.vec<-vector(length=demograph_out$total.pop.num);
		col.arrow.vec<-vector(length=demograph_out$total.pop.num);
		for(i in 1:demograph_out$total.pop.num){
			col.pop.vec[i]=col.pop[(i-1)%%length(col.pop)+1];
			col.arrow.vec[i]=col.arrow[(i-1)%%length(col.arrow)+1];
		}


#plot region

		if(size.scale=="topology"){
#time.series<-time.series+(1:length(time.series));
			time<-1:length(time);
			ylab="";
			plot_N<-0;
		}
		if(is.null(ylim)){
			y_min=min(time);
			y_max=max(time)*1.3;
			ylim=c(y_min-y_max*0.05, y_max);}
		else {ylim=ylim;}
		if(time.scale=="log10year")ylab="Time before present (yr)";
		if(time.scale=="kyear")ylab="Time before present (kyr)";

		if(size.scale=="linear")plot_N<-max(demograph_out$N)*linear.scale;
		if(size.scale=="log")plot_N<-log(max(demograph_out$N)+1, base=log.base);
		if(is.null(xlim)){
			x_max=max(demograph_out$Pos)+max(plot_N, (max(demograph_out$Pos)-min(demograph_out$Pos))*0.05);
			x_min=min(demograph_out$Pos)-max(plot_N, (max(demograph_out$Pos)-min(demograph_out$Pos))*0.05);
			xlim=c(x_min, x_max);}
		else {xlim=xlim;}

#population labels and positions
		if(length(pops)==demograph_out$total.pop.num){
			lab.pop<-pops;
		} else { lab.pop<-1:demograph_out$total.pop.num; }

		lab.pos<-cbind(demograph_out$pop.pos, demograph_out$pop.pos*0-y_max*0.05+y_min);


		list(size.scale=size.scale, linear.scale=linear.scale, log.base=log.base, 
				time.scale=time.scale, time=time,
				col.pop=col.pop.vec,col.arrow=col.arrow.vec,
				length.arrowtip=length.arrowtip, lwd.arrow=lwd.arrow,
				angle.arrowtip=angle.arrowtip,
				lab.pop=lab.pop,lab.pos=lab.pos,
				xlim=xlim, ylim=ylim,
				xlab=xlab,ylab=ylab,
				cex.lab=cex.lab, cex.axis=cex.axis,
				axes=axes)
}


ReadMigPar<-function(demograph_out, 
		size.scale="linear", linear.scale=0.2,log.base=10,
		time.scale="4Ne",
		col.pop="gray45", cex.lab=1.5,
		col.arrow="black",length.arrowtip=0.15, lwd.arrow=1,
		angle.arrowtip=15, pops=demograph_out$pop.lab)
{
## three plot modes can be used to adjust the presentation of population size in the plot, the width of each population lineage. size.scale="linear" gives the proportion of linear.scale to original width; size.scale="log" gives the logarithms of base log.base of original width, and size.scale="topology" will treat the width to be zero.
## time.scale gives the unit in the y-axis, has the option of "year", "generation", or with default value as the ms' 4Ne scaled time.
	if(size.scale!="linear"&&size.scale!="log"&&size.scale!="topology")print("Error, you need to redefine the parameter size.scale");
	if(time.scale!="4Ne"&time.scale!="generation"&time.scale!="year"&time.scale!="kyear"&time.scale!="log10year"){print("Error, you need to redefine the time.scale");}
#re-scale time
	time<-RescaleT(demograph_out, time.scale=time.scale);

#color vector
	col.pop.vec<-vector(length=demograph_out$total.pop.num);
	col.arrow.vec<-vector(length=demograph_out$total.pop.num);
	for(i in 1:demograph_out$total.pop.num){
		col.pop.vec[i]=col.pop[(i-1)%%length(col.pop)+1];
		col.arrow.vec[i]=col.arrow[(i-1)%%length(col.arrow)+1];
	}
#migrations 
	events<-1;
	count<-0;                                
	if(length(time)>=2){
		for(i in 2:length(time)){
			if(abs(sum(demograph_out$m[,,i]-demograph_out$m[,,i-1]))>0)events<-c(events, i);
		}
	}
	count<-length(events);
	if(count>=2){
		n<-floor(sqrt(count));
		if(n*n<count & n*(n+1)>=count){mfrow=c(n, n+1);}
		if(n*(n+1)<count & (n+1)*(n+1)>=count){mfrow=c(n+1, n+1);}
		if(n*n==count){mfrow=c(n, n);}
	}
	else {mfrow=c(1, 1);}
#plot region
	xlim=c(-3.5, 3.5); ylim=c(-3.5, 3.5);

#population labels
	if(length(pops)==demograph_out$total.pop.num){
		lab.pop<-pops;
	} else { lab.pop<-1:demograph_out$total.pop.num; }


	list(size.scale=size.scale, linear.scale=linear.scale, log.base=log.base, 
			time.scale=time.scale, time=time, lab.pop=lab.pop,
			col.pop=col.pop.vec,col.arrow=col.arrow.vec,
			length.arrowtip=length.arrowtip, lwd.arrow=lwd.arrow,
			angle.arrowtip=angle.arrowtip,
			xlim=xlim, ylim=ylim,
			mfrow=mfrow, events=events,
			cex.lab=cex.lab)
}
ReadDemoGraph<-function(inputfile, type, inpos=NULL, N4, gen, pops)
{
## inputfile is the language file describing the demograph history
## type is the type of language in file, currently support ms, msa, MaCs, and cosi.
## inpos is a vector positions for each population on the plot at time 0.
## N4 is the effective population size of 4Ne, used to rescale the time and geneflow strength


	demograph_out<-list(time.series=c(), Pos=c(), N=c(), m=c(), survive=c(), g.rate=c(), Pos_map=c(), pop.pos=inpos, pop.lab=pops, mscmd=c(), present.pop.num=0, total.pop.num=0, N4=N4, gen=gen);

#translate the coomand to ms languiage
	x<-readcmd(cmdfile=inputfile, type=type, N4=N4);
	demograph_out$mscmd<-x$mscmd;
	if(is.null(pops)&type=="cosi")demograph_out$pop.lab=x$pop.lab;
#	L=length(mscmd);

#the number of populations for simulation
	x<-readpopnum(demograph_out);
	demograph_out$present.pop.num=x$present.pop.num;
	demograph_out$total.pop.num=x$total.pop.num;


#read the time series
	demograph_out$time.series<-readtime(demograph_out);
	time_num=length(demograph_out$time.series);

	message(paste("There are ",time_num,"time events for ",demograph_out$total.pop.num," populations"));

	x<-readNg(demograph_out);
	demograph_out$N<-x$N;
	demograph_out$g.rate<-x$g.rate;
	message("read N and g, done!");
	demograph_out$m<-readm(demograph_out);
	message("read m, done!");

#population positions
	if(is.null(inpos))demograph_out$pop.pos<-(1:demograph_out$total.pop.num) * max(demograph_out$N, na.rm=T);

#population positions and all
	x<-readpos(demograph_out);
	message("read pos and update, done!");
	demograph_out$time.series<-x$time.series;
	demograph_out$survive<-x$survive;
	demograph_out$Pos<-x$Pos;
	demograph_out$g.rate<-x$g.rate;
	demograph_out$m<-x$m;
	demograph_out$N<-x$N;
	demograph_out$Pos_map<-x$Pos_map;

#update all parameters
	x<-update(demograph_out);
	demograph_out$g.rate<-x$g.rate;
	demograph_out$m<-x$m;
	demograph_out$N<-x$N;
	if(type=="scrm"){demograph_out$m<-updatescrm(demograph_out);}##for the command of -eps


	demograph_out
}
