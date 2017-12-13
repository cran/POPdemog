iscmd<-function(str)
{
##to determine whether str is a ms command
	cmdlist<-c("-I", "-m", "-ma", "-eM", "-ema", "-em", "-G", "-g", "-n", "-eG", "-eg", "-eN", "-en", "-ej", "-es", "-eA", "-f", "-seeds", "-t", "-s", "-T", "-L", "-p", "-r", "-c");
	for(cmd in cmdlist){
		if(str==cmd)return(str==cmd);
	}
	return(FALSE);
}

fmatch<-function(x, y)
{
#print(c(x,y));
	if(abs(x-y)<0.000000000001){
		return(TRUE);
	} else {
		return(FALSE);
	}
}

TransT<-function(time, scale1, scale2, demograph_out)
{
	#print(paste("translate time from", scale1,time,"to", scale2))
##scale1 to years
	if(scale1=="4Ne"){
		x<-time*demograph_out$N4*demograph_out$gen;
	}
	else if(scale1=="generation"){
		x<-time*demograph_out$gen;
	}
	else if(scale1=="year"){
		x<-time;
	}
	else if(scale1=="kyear"){
		x<-time*1000;
	}
	else if(scale1=="log10year"){
		x<-10^time;
	}
	else {
		stop(paste("Error, you need to redefine the time.scale:", scale1));
	}
	## years to scale2
	if(scale2=="4Ne"){
                time<-x/demograph_out$N4/demograph_out$gen;
        }
        else if(scale2=="generation"){
                time<-x/demograph_out$gen;
        }
        else if(scale2=="year"){
                time<-x;
        }
        else if(scale2=="kyear"){
                time<-x/1000;
        }
        else if(scale2=="log10year"){
                time<-log(x, base=10);
        }
        else {
                stop(paste("Error, you need to redefine the time.scale:", scale2));
        }
	time
}

RescaleT<-function(demograph_out, time.scale, time=demograph_out$time.series)
{
	if(time.scale=="4Ne"){
		time<-time;
	}
	else if(time.scale=="generation"){
		time<-time*demograph_out$N4;
	}
	else if(time.scale=="year"){
		time<-time*demograph_out$N4*demograph_out$gen;
	}
	else if(time.scale=="kyear"){
		time<-time*demograph_out$N4*demograph_out$gen/1000;
	}
	else if(time.scale=="log10year"){
		time<-log10(time*demograph_out$N4*demograph_out$gen+1);
	}
	else {print("Error, you need to redefine the time.scale");}
	time
}

#time, Pos, Pos_map, N, survive, g.rate, ind, col, smooth=TRUE, plotmode=plotmode)
plotpath<-function(demograph_out, evo_par, ind, 
		time.scale=evo_par$time.scale, size.scale=evo_par$size.scale,
		col=evo_par$col.pop, smooth=TRUE)
{

	time<-demograph_out$time;
	g.rate<-demograph_out$g.rate;
	if(size.scale=="topology")time<-1:length(time);

	Pos=demograph_out$Pos;
	Pos_map=demograph_out$Pos_map;
	N=demograph_out$N;
	survive=demograph_out$survive;



	scale<-0.85;
	dif1<-max(time)*0.005;
	L_time<-length(time);
	life<-survive[1, ind]:survive[2, ind];
	L_life<-length(life);	
	width<-N[life[1], ind]/2;
	pos<-Pos[life[1], ind];
	lne<-c(time[life[1]], pos, width, -1, 1);
	if(L_life==1){
		lne<-c(time[life[1]], pos, 0, -1, 1);
		return(rbind(lne, lne));
	}
	for(i in life[2:L_life])
	{
		if(!fmatch(Pos[i, ind], Pos[i-1, ind])|!fmatch(N[i, ind], N[i-1, ind])){
			if(!smooth | (smooth & fmatch(N[i, ind], N[i-1, ind]))){
				pos<-Pos[i-1, ind];
				dif2<-(time[i]-time[i-1])*(1-scale);
				if(dif1<dif2){time_dif<-dif1;}
				else {time_dif<-dif2;}
				newtime<-time[i]-time_dif;

				width<-N[i-1, ind]*exp(-(newtime-time[i-1])*g.rate[i-1,ind])/2;
				newlne<-c(newtime, pos, width, -1, 1);
				lne<-rbind(lne, newlne);

				pos<-Pos[i, ind];
				width<-N[i, ind]/2;
				newlne<-c(time[i], pos, width, -1, 1);
				lne<-rbind(lne, newlne);
			}
			else {
				pos<-Pos[i-1, ind];
				dif1<-time[i]-time[i-1];
				time_steps<-seq(time[i-1], time[i]-0.01*dif1, 0.01*dif1);
				for(newtime in time_steps){
					width<-N[i-1, ind]*exp(-(newtime-time[i-1])*g.rate[i-1,ind])/2;
					newlne<-c(newtime, pos, width, -1, 1);
					lne<-rbind(lne, newlne);
				}
				pos<-Pos[i, ind];
				width<-N[i, ind]/2;
				newlne<-c(time[i], pos, width, -1, 1);
				lne<-rbind(lne, newlne);
			}

		} else {
			pos<-Pos[i, ind];
			width<-N[i, ind]/2;
			newlne<-c(time[i], pos, width, -1, 1);
			lne<-rbind(lne, newlne);
		}

	}
	if(survive[2, ind]==L_time){
		i<-L_time;
		newtime<-time[i]*1.01;
#pos<-Pos[i+1, ind];
		pos<-Pos[i+1, Pos_map[i+1, ind]];
#width<-N[i, ind]*exp(-time[i]*0.3*g.rate[i,ind])/2;
		width<-N[i, Pos_map[i+1, ind]]*exp(-time[i]*0.01*g.rate[i,Pos_map[i+1, ind]])/2;
		old_pos<-Pos[i,ind];
#if(pos>old_pos)newlne<-c(newtime, pos, width, -1, -1);
#if(pos<old_pos)newlne<-c(newtime, pos, width, 1, 1);
		if(fmatch(pos, old_pos)){
			newlne<-c(newtime, pos, width, -1, 1);
			lne<-rbind(lne, newlne);
			newlne<-c(newtime*1.3, pos, width, -1, 1);

		}

#newlne<-c(newtime, pos-width, pos+width);
		lne<-rbind(lne, newlne);
	}

	if(life[1]>1){
		i<-life[1];

		dif2<-(time[i]-time[i-1])*(1-scale);
		if(dif1<dif2){time_dif<-dif1;}
		else {time_dif<-dif2;}
		newtime<-time[i]-time_dif;


#pos<-Pos[i-1, ind];
		pos<-Pos[i-1, Pos_map[i-1, ind]];
		width<-N[i, Pos_map[i-1, ind]]/2;
		old_pos<-Pos[i,ind];
		if(pos>old_pos)newlne<-c(newtime, pos, width, -1, -1);
		if(pos<old_pos)newlne<-c(newtime, pos, width, 1, 1);
#		lne<-rbind(newlne, lne);
	}
#print(lne[1,])
	if(size.scale!="topology")lne[,1]<-RescaleT(demograph_out, time.scale=time.scale, time=lne[,1]);

	if(life[L_life]<L_time){
		i<-life[L_life];

		dif2<-(time[i+1]-time[i])*(1-scale);
		if(dif1<dif2){time_dif<-dif1;}
		else {time_dif<-dif2;}
		newtime<-time[i]+time_dif;

		pos<-Pos[i+1, Pos_map[i+1, ind]];
		width<-N[i, Pos_map[i+1, ind]]/2;
		old_pos<-Pos[i,ind];		
		if(pos>old_pos)newlne<-c(newtime, pos, width, -1, -1);
		if(pos<old_pos)newlne<-c(newtime, pos, width, 1, 1);
#		lne<-rbind(lne, newlne);
	}
	if(size.scale=="topology"){
		l1<-lne[,c(2,1)];
		l2<-lne[,c(2,1)];
		lines(l1, col=col[ind], lwd=1.5);
		lines(l2, col=col[ind], lwd=1.5);
		lne[,3]<-0;
	}
	if(size.scale=="linear"){
		l1<-cbind(lne[,2]+lne[,3]*lne[,4]*evo_par$linear.scale,lne[,1]);
		l2<-cbind(lne[,2]+lne[,3]*lne[,5]*evo_par$linear.scale,lne[,1]);
		lne[,3]<-lne[,3]*evo_par$linear.scale;
		lines(l1, col=col[ind]);
		lines(l2, col=col[ind]);

	}
	if(size.scale=="log"){
		l1<-cbind(lne[,2]+log(lne[,3]+1, base=evo_par$log.base)*lne[,4],lne[,1]);
		l2<-cbind(lne[,2]+log(lne[,3]+1, base=evo_par$log.base)*lne[,5],lne[,1]);
		lne[,3]<-log(lne[,3]+1, base=evo_par$log.base);
		lines(l1, col=col[ind]);
		lines(l2, col=col[ind]);

	}
	ll1<-c(l1[,1], rev(l2[,1]),l1[1,1]);
	ll2<- c(l1[,2], rev(l2[,2]),l1[1,2]);
	polygon(ll1, ll2, col=col[ind], border=col[ind]);
	return(lne);
}

ScalePlotRegion<-function()
{
	pin<-par("pin");#plot region ratio
        usr<-par("usr");#numeric ratio
        ra1<-pin[1]/pin[2];
        ra2<-(usr[4]-usr[3])/(usr[2]-usr[1]);
	rx<-1;
	ry<-ra1*ra2;
	if(ry>rx){
		rx <- 1/ry;
		ry <- 1;
	}
	c(rx, ry)
}

addcircle<-function(x, y, col, radius){
	r<-ScalePlotRegion();
	rx<-r[1];ry<-r[2];
	theta<-seq(0, 360, 1);
	l1<-radius*cos(theta)*rx+x;
	l2<-radius*sin(theta)*ry+y;
	polygon(l1, l2, col=col, border=col);
}

PlotCircle<-function(ind, time_pt, event_pt, mig_par, demograph_out,
		col=mig_par$col.pop, size.scale=mig_par$size.scale, 
		linear.scale=mig_par$linear.scale, log.base=mig_par$log.base,
		toposize.scale=1, xlim = xlim, ylim = ylim,
		add=FALSE, map.pos=NULL)
{
	r<-ScalePlotRegion();
        rx<-r[1];ry<-r[2];
	#center of the plot
	lx<-abs(xlim[2]-xlim[1]);
	ly<-abs(ylim[2]-ylim[1]);
	if(lx>=ly){
		cx<-xlim[2]-ly/2;
		cy<-mean(ylim);
	}
	else {
		cx<-mean(xlim);
		cy<-ylim[1]+lx/2;
	}
	#circle radius to put population deme
	Rad<-min(lx, ly)*3/10/max(rx,ry);
	#print(c(Rad, xlim, ylim, cx, cy))
	#N=demograph_out$N[event, ];
	N<-NOut(time_pt, demograph_out, time.scale=mig_par$time.scale);
	total.pop.num=demograph_out$total.pop.num;
#lab=mig_par$lab.pop[ind];	
	if(add==FALSE){
		x=Rad*cos((ind-1)/total.pop.num*2*pi)*rx+cx;
		y=Rad*sin((ind-1)/total.pop.num*2*pi)*ry+cy;
	}
	else {
		x=map.pos[ind, 1];y=map.pos[ind, 2];
	}	
	if(size.scale=="topology")radius=sqrt(2/total.pop.num*toposize.scale);
	if(size.scale=="linear")radius=sqrt(N[ind]*linear.scale);
	if(size.scale=="log")radius=sqrt(log(1+N[ind]/2,base=log.base));
	if(event_pt<=demograph_out$survive[2, ind]&event_pt>=demograph_out$survive[1, ind])addcircle(x=x, y=y, col=col[ind], radius=radius);
	return(c(x, y, radius));
}

CalArrowPosCircle<-function(x1,y1,r1, x2,y2,r2){
#from pos 1 to pos 2
	r<-ScalePlotRegion();
        rx<-r[1];ry<-r[2];
	dis<-sqrt((x1-x2)^2/(rx^2)+(y1-y2)^2/(ry^2));
	#points(c(x1,x2), c(y1, y2), pch=19)
	ratio_x<-(x2-x1)/dis;
	ratio_y<-(y2-y1)/dis;
	x_f<-x1+ratio_x*r1;
	y_f<-y1+ratio_y*r1;
	x_t<-x2-ratio_x*r2;
	y_t<-y2-ratio_y*r2;
	return(c(x_f, y_f, x_t, y_t));
}

cal_arrow_pos<-function(ind2, ind1, time, list_pos, m)
{
#from ind2 to ind1
	time_num<-length(time);
#step<-max(time[min(survive[2,ind1], survive[2,ind2])]-time[max(survive[1,ind1], survive[1,ind2])])/2;
	data1<-list_pos[[ind1]];
	data2<-list_pos[[ind2]];
	time_ps<-c();
	f_strength<-c();

	if(FALSE)
	{
		for(i in 2:time_num){
			if(m[ind1, ind2, i-1]>0){
				if(m[ind1, ind2, i-1]<1){step=(time[i]-time[i-1])/2.1;}
				else {step=time[i]-time[i-1]+1;}
				time_ps<-c(time_ps, seq(time[i-1], time[i], step));
				f_strength<-c(f_strength, seq(time[i-1], time[i], step)*0+m[ind1, ind2, i-1]);
			}
		}

	}
##generate the time points of migration rate changes
	time_sparse<-c(1);
	i=2;
	while(i<=time_num){
		L<-time_sparse[length(time_sparse)];
		if(!fmatch(m[ind1, ind2, i], m[ind1, ind2, L]))time_sparse<-c(time_sparse, i);
		i=i+1;
	}
	if(time_sparse[length(time_sparse)]!=time_num)time_sparse<-c(time_sparse, time_num);

#print(time_sparse);
#print(m[ind1, ind2, time_sparse]);

	for(i in 2:length(time_sparse)){
		a<-time_sparse[i-1];
		b<-time_sparse[i];
		if(m[ind1, ind2, a]>0){
			if(m[ind1, ind2, a]<1){
step=max((time[b]-time[a])/2, 1.001*time[time_num]/time_num);
}
			else {step=time[b]-time[a]+1;}

			time_ps<-c(time_ps, seq(time[a], time[b], step));
			f_strength<-c(f_strength, seq(time[a], time[b], step)*0+m[ind1, ind2, a]);
		}
	}
#print(time_ps);
	if(m[ind1, ind2, time_num]>0){
#print(c(time[time_num], min(max(data1[,1]), max(data2[,1])), step))
		step=(min(max(data1[,1]), max(data2[,1]))-time[time_num])/3;
		if(time[time_num]<=min(max(data1[,1]), max(data2[,1]))){
			time_ps<-c(time_ps,seq(time[time_num], min(max(data1[,1]), max(data2[,1])), step));
			f_strength<-c(f_strength, seq(time[time_num], min(max(data1[,1]), max(data2[,1])), step)*0+m[ind1, ind2, time_num]);
		}
	}
#print(time_ps);
	time_ps<-time_ps[time_ps<=max(data1[,1])&time_ps<max(data2[,1])];
	f_strength<-f_strength[time_ps<=max(data1[,1])&time_ps<max(data2[,1])];
	L<-length(time_ps);
	if(L>0){
		out<-matrix(NA, nrow=L, ncol=4);
		for(i in 1:L){
			out[i,1]<-time_ps[i];
			out[i,4]<-f_strength[i];
			x<-findpos(p=time_ps[i], data=data1);
			pos1<-x[1];wid1<-x[2];
			x<-findpos(p=time_ps[i], data=data2);
			pos2<-x[1];wid2<-x[2];
#print(data1);
#print(data2);
#print(c(pos1, wid1, pos2, wid2))
			if(pos1>pos2){out[i,2]=pos2+wid2;out[i,3]=pos1-wid1;}
			if(pos1<pos2){out[i,2]=pos2-wid2;out[i,3]=pos1+wid1;}
			if(pos1==pos2){out[i,2]=NA;out[i,3]=NA;}
		}
	} else {
		out<-matrix(NA, nrow=1, ncol=4);
	}
	return(out);
}

findpos<-function(p, data){
#print(p)
#	print(data)
	L<-dim(data)[1];
	for(i in 2:L){
		if(p>=data[i-1,1]&p<data[i,1]){
			ratio<-(p-data[i-1,1])/(data[i,1]-data[i-1,1]);
			pos<-ratio*(data[i,2]-data[i-1,2])+data[i-1,2];
			wid<-ratio*(data[i,3]-data[i-1,3])+data[i-1,3];
			return(c(pos, wid));
		}
	}
	if(p>=data[L,1]){pos<-data[L,2];wid<-data[L,3];}
	if(p<data[1,1]){pos<-data[1,2];wid<-data[1,3];}
	return(c(pos, wid));
}

PlotEvo<-function(demograph_out, evo_par,
		size.scale=evo_par$size.scale, time.scale=evo_par$time.scale,
		col.pop=evo_par$col.pop,col.arrow=evo_par$col.arrow,
		xlim=evo_par$xlim, ylim=evo_par$ylim,
		xlab=evo_par$xlab,ylab=evo_par$ylab,
		cex.lab=evo_par$cex.lab, cex.axis=evo_par$cex.axis,
		length.arrowtip=evo_par$length.arrowtip, angle.arrowtip = evo_par$angle.arrowtip,
		lwd.arrow=evo_par$lwd.arrow,
		lab.pos=evo_par$lab.pos, lab.pop=evo_par$lab.pop, axes=evo_par$axes,
		m.adjust=0
		)
{
	m<-demograph_out$m;
	m[m<m.adjust]=0;

	time_num=length(evo_par$time);
#print(col.arrow)
#lwd.m<-vector(length=total.pop.num)*0+lwd.arrow;
	if(time_num>1){

		plot(NA,NA, xlim=xlim, ylim=ylim, axes=F, ylab=ylab, xlab="", xpd=TRUE, cex.lab=cex.lab);
		title(xlab=xlab, line=1, cex.lab=cex.lab);
		list_pos<-vector("list", demograph_out$total.pop.num);
		for(i in 1:demograph_out$total.pop.num){
#x<-plotpath(time=time.series, Pos=Pos, Pos_map=Pos_map, N=N, g.rate=g.rate, survive=survive, ind=i, col=col.tree, plotmode=plotmode);
			x<-plotpath(demograph_out=demograph_out, evo_par=evo_par, ind=i, col=col.pop, size.scale=size.scale);
			list_pos[[i]]<-x;
		}
		for(ind1 in 1:demograph_out$total.pop.num){
			for(ind2 in 1:demograph_out$total.pop.num){
##from ind2 to ind1
				out<-cal_arrow_pos(ind2=ind2, ind1=ind1, time=evo_par$time, list_pos=list_pos, m=m);
				out<-na.omit(out);
				if(dim(out)[1]>0){arrows(out[,2], out[,1], out[,3], out[,1], lwd=0.5+out[,4]*lwd.arrow, length=length.arrowtip, angle = angle.arrowtip, col=col.arrow[ind2], code=2);
#print(out);
#print(length.arrowtip)
				}
#print(dim(out))
			}
		}
#lab_pos<-cbind(inpos, inpos*0);
		if(axes==TRUE){
			if(size.scale!="topology"){
				if(time.scale!="log10year")axis(2, cex.axis=cex.axis);
				if(time.scale=="log10year"){
					time_log_L<-floor(max(ylim));
					y.log.ticks<-0:time_log_L;
					y.log.label<-parse(text=paste("10^",y.log.ticks, sep=""));
					y.log.label[1]="0";
					axis(side=2, at=y.log.ticks, labels=y.log.label, cex.axis=cex.axis);
				}
			}
			text(lab.pos, labels=lab.pop, cex=cex.lab);
		}
	}
	if(time_num==1){
		print("You do not need to plot a evolutionary graph, please try function 'plotmig'.");
	}
}


