findpopind<-function(pop, pop_index)
{
	for(i in 1:length(pop_index)){
		if(pop==pop_index[i])return(i);
	}
	return(NA); 
}

readcosi<-function(cosifile, N4)
{
	tag_n<-0;
	npop<-0;
	pop_index<-c();
	pop.lab<-c();
	con<-file(cosifile, "r");
	while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
# do stuff
		input_vec <- (strsplit(oneLine, " "));
		if(length(input_vec[[1]])>1){
			if(input_vec[[1]][1]=="pop_define"){
				pop_index<-c(pop_index, input_vec[[1]][2]);
				pop.lab<-c(pop.lab, input_vec[[1]][3]);
				npop<-npop+1;
			}
			if(input_vec[[1]][1]=="sample_size"){
				popid<-findpopind(pop=input_vec[[1]][2], pop_index=pop_index);
			}
		}
	}
	close(con);


	con<-file(cosifile, "r");
	sample_size<-vector(length=npop);
	while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
# do stuff
		input_vec <- (strsplit(oneLine, " "));
		if(length(input_vec[[1]])>1){
			if(input_vec[[1]][1]=="sample_size"){
				popid<-findpopind(pop=input_vec[[1]][2], pop_index=pop_index);
				sample_size[popid]<-input_vec[[1]][3];
			}
		}
	}
	close(con);



	#mscmd<-c("-I", npop, rep(10, npop));
	mscmd<-c("-I", npop, sample_size);

	con<-file(cosifile, "r");
	while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
# do stuff
		input_vec <- (strsplit(oneLine, " "));
		if(length(input_vec[[1]])>1){
			if(input_vec[[1]][1]=="pop_size"){
				popid<-findpopind(pop=input_vec[[1]][2], pop_index=pop_index);
				size<-as.numeric(input_vec[[1]][3])/N4;
				mscmd<-c(mscmd, "-n", popid, size);
			}
			if(input_vec[[1]][1]=="pop_event"){
				if(input_vec[[1]][2]=="change_size"){
					L<-length(input_vec[[1]]);
					t<-as.numeric(input_vec[[1]][L-1])/N4;
					size<-as.numeric(input_vec[[1]][L])/N4;
#popid<-input_vec[[1]][L-2];
					popid<-findpopind(pop=input_vec[[1]][L-2], pop_index=pop_index);
					mscmd<-c(mscmd, "-en", t, popid, size);}
				if(input_vec[[1]][2]=="exp_change_size"){
					L<-length(input_vec[[1]]);
#popid<-input_vec[[1]][L-4];
					popid<-findpopind(pop=input_vec[[1]][L-4], pop_index=pop_index);
					t_start<-as.numeric(input_vec[[1]][L-2])/N4;
					t_end<-as.numeric(input_vec[[1]][L-3])/N4;
					size_start<-as.numeric(input_vec[[1]][L])/N4;
					size_end<-as.numeric(input_vec[[1]][L-1])/N4;
					g<- log(size_end/size_start)/(t_start-t_end);
					mscmd<-c(mscmd, "-en", t_end, popid, size_end);
					mscmd<-c(mscmd, "-eg", t_end, popid, g);
					mscmd<-c(mscmd, "-eg", t_start, popid, 0);
				}
				if(input_vec[[1]][2]=="bottleneck"){
					L<-length(input_vec[[1]]);
#popid<-input_vec[[1]][L-2];
					popid<-findpopind(pop=input_vec[[1]][L-2], pop_index=pop_index);
					size<-1/as.numeric(input_vec[[1]][L])/N4;
					t<-as.numeric(input_vec[[1]][L-1])/N4;
					mscmd<-c(mscmd, "-en", t, popid, size);
				}
				if(input_vec[[1]][2]=="migration_rate"){
					L<-length(input_vec[[1]]);
#popid1<-input_vec[[1]][L-3];#1->2
					popid1<-findpopind(pop=input_vec[[1]][L-3], pop_index=pop_index);
#popid2<-input_vec[[1]][L-2];
					popid2<-findpopind(pop=input_vec[[1]][L-2], pop_index=pop_index);
					t<-as.numeric(input_vec[[1]][L-1])/N4;
					m<-as.numeric(input_vec[[1]][L])*N4;
					mscmd<-c(mscmd, "-em", t, popid2, popid1, m);
				}
				if(input_vec[[1]][2]=="split"){
					L<-length(input_vec[[1]]);
#popid1<-input_vec[[1]][L-2];#2->1;
					popid1<-findpopind(pop=input_vec[[1]][L-2], pop_index=pop_index);
#popid2<-input_vec[[1]][L-1];
					popid2<-findpopind(pop=input_vec[[1]][L-1], pop_index=pop_index);
					t<-as.numeric(input_vec[[1]][L])/N4;
					mscmd<-c(mscmd, "-ej", t, popid2, popid1);
				}
				if(input_vec[[1]][2]=="admix"){
					L<-length(input_vec[[1]]);
#popid1<-input_vec[[1]][L-3];#2->1;
					popid1<-findpopind(pop=input_vec[[1]][L-3], pop_index=pop_index);
#popid2<-input_vec[[1]][L-2];
					popid2<-findpopind(pop=input_vec[[1]][L-2], pop_index=pop_index);
					t<-as.numeric(input_vec[[1]][L-1])/N4;
					m<-as.numeric(input_vec[[1]][L])*N4;
					mscmd<-c(mscmd, "-em", t, popid1, popid2, m);
					mscmd<-c(mscmd, "-em", t+1/N4, popid1, popid2, 0);
				}
			}
		}
	} 
	close(con);
	list(mscmd=mscmd, pop.lab=pop.lab);
}

cleanmsHot<-function(imput.cmd)
{
out<-c()
i=1
while(i<=length(imput.cmd)){
if(imput.cmd[i]!="-v"&&imput.cmd[i]!="-V"){out<-c(out, imput.cmd[i]);i=i+1;}
else {
n<-as.numeric(imput.cmd[i+1]);
i=i+n*3+2;
}
}
out
}

#cmd<-readcosi("/projects/browning/brwnlab/joe/projects/plotms/plotms.master/examples/cosi2/software/cosi-2.0/examples/bestfit/params");
readcmd<-function(cmdfile, type, N4)
{
	types="ms, msa, mshot, msprime, macs, scrm, cosi"
	if(type=="msprime"){
		input.cmd<-scan(cmdfile, what=character(),  comment.char = "\\");
		input.cmd[input.cmd=="--structure"]<-"-I";
		input.cmd[input.cmd=="--migration-matrix-entry"]<-"-m";
		input.cmd[input.cmd=="--migration-matrix"]<-"-ma";
		input.cmd[input.cmd=="--migration-rate-change"]<-"-eM";
		input.cmd[input.cmd=="--migration-matrix-change"]<-"-ema";
		input.cmd[input.cmd=="--migration-matrix-entry-change"]<-"-em";
		input.cmd[input.cmd=="--growth-rate"]<-"-G";
		input.cmd[input.cmd=="--population-growth-rate"]<-"-g";
		input.cmd[input.cmd=="--population-size"]<-"-n";
		input.cmd[input.cmd=="--growth-rate-change"]<-"-eG";
		input.cmd[input.cmd=="--population-growth-rate-change"]<-"-eg";
		input.cmd[input.cmd=="--size-change"]<-"-eN";
		input.cmd[input.cmd=="--population-size-change"]<-"-en";
		input.cmd[input.cmd=="--population-split"]<-"-ej";
		input.cmd[input.cmd=="--admixture"]<-"-es";
		list(mscmd=input.cmd);
	}
	else if(type=="macs"|type=="msa"){
		input.cmd<-scan(cmdfile, what=character(),  comment.char = "\\");
		list(mscmd=input.cmd)
	}
	else if(type=="ms"|type=="scrm")
	{
		input.cmd<-scan(cmdfile, what=character(),  comment.char = "\\");
		input.cmd[input.cmd=="-eA"]<-"-XX";
		list(mscmd=input.cmd)
	}
	else if(type=="cosi")
	{
		x<-readcosi(cmdfile, N4);
		list(mscmd=x$mscmd, pop.lab=x$pop.lab)
	}
	else if(type =="mshot")
	{
		input.cmd<-scan(cmdfile, what=character(),  comment.char = "\\");
		input.cmd<-cleanmsHot(input.cmd);
		list(mscmd=input.cmd)
	}
	else
	{
		print("Error: unsupported command type");
		print(paste("currently support:", types));
	}

}


readpopnum<-function(demograph_out)
{
	ms.cmd=demograph_out$mscmd;
	present.pop.num=1;
	total.pop.num=0;;
	L=length(ms.cmd);
	for(i in 1:L){
		if(ms.cmd[i]=="-I"){
			present.pop.num=as.numeric(ms.cmd[i+1]);
			total.pop.num=total.pop.num+present.pop.num;
		} else if(ms.cmd[i]=="-es"){
			total.pop.num=total.pop.num+1;
		}
	}
	if(total.pop.num==0)total.pop.num=total.pop.num+present.pop.num;
#print(c(present.pop.num, total.pop.num));
	return(list(present.pop.num=present.pop.num, total.pop.num=total.pop.num));
}


readtime<-function(demograph_out)
{
##read the time points in the ms.cmd
	ms.cmd=demograph_out$mscmd;
	time.series<-c(0);
	L=length(ms.cmd);
	for(i in 1:L){
		if(ms.cmd[i]=="-eG"){time.series=c(time.series, ms.cmd[i+1]);}
		else if(ms.cmd[i]=="-eg"){time.series=c(time.series, ms.cmd[i+1]);}
		else if(ms.cmd[i]=="-eN"){time.series=c(time.series, ms.cmd[i+1]);}
		else if(ms.cmd[i]=="-en"){time.series=c(time.series, ms.cmd[i+1]);}
		else if(ms.cmd[i]=="-eM"){time.series=c(time.series, ms.cmd[i+1]);}
		else if(ms.cmd[i]=="-em"){time.series=c(time.series, ms.cmd[i+1]);}
		else if(ms.cmd[i]=="-ema"){time.series=c(time.series, ms.cmd[i+1]);}
		else if(ms.cmd[i]=="-es"){time.series=c(time.series, ms.cmd[i+1]);}
		else if(ms.cmd[i]=="-eps"){time.series=c(time.series, ms.cmd[i+1]);}
		else if(ms.cmd[i]=="-ej"){time.series=c(time.series, ms.cmd[i+1]);}
		else if(ms.cmd[i]=="-eA"){time.series=c(time.series, ms.cmd[i+1]);}
	}
	time.series<-unique(as.numeric(time.series));
	time.series<-sort(time.series);
	return(time.series);
}

readNg<-function(demograph_out)
{
	N0=1;ind=0;
	mscmd=demograph_out$mscmd;
	time=demograph_out$time.series;
	total.pop.num=demograph_out$total.pop.num;
	present.pop.num=demograph_out$present.pop.num;

	L<-length(mscmd);
	time_num<-length(time);
	N<-matrix(-9, nrow=time_num, ncol=total.pop.num);
	for(i in 1:present.pop.num){N[1, i]=N0;}	
	g.rate<-matrix(-9, nrow=time_num, ncol=total.pop.num);
	for(i in 1:present.pop.num){g.rate[1, i]=0;}	
	for(iter in 1:time_num){
		for(i in 1:L){
			if(mscmd[i]=="-G"){
				tmpd=as.numeric(mscmd[i+1]);#i=i+1;
				for(k in 1:total.pop.num){
					g.rate[1,k]=tmpd;
				}
			}
			else if(mscmd[i]=="-n"){
				ind=as.integer(mscmd[i+1]);i=i+1;
				tmpd=as.numeric(mscmd[i+1]);#i=i+1;
				N[1,ind]=tmpd*N0;
			}
			else if(mscmd[i]=="-g"){
				ind=as.integer(mscmd[i+1]);i=i+1;
				tmpd=as.numeric(mscmd[i+1]);#i=i+1;
				g.rate[1,ind]=tmpd;
			}
			else if(mscmd[i]=="-eG"){
				tmpd=as.numeric(mscmd[i+1]);i=i+1;
				if(fmatch(x=tmpd,y=time[iter])){
#print(mscmd[i-1]);
					tmpd=as.numeric(mscmd[i+1]);#i=i+1;
					for(k in 1:total.pop.num){
						g.rate[iter,k]=tmpd;
					}
				}
			}
			else if(mscmd[i]=="-eg"){
				tmpd=as.numeric(mscmd[i+1]);i=i+1;
				if(fmatch(x=tmpd,y=time[iter])){
#print(mscmd[i-1]);
					ind=as.integer(mscmd[i+1]);i=i+1;
					tmpd=as.numeric(mscmd[i+1]);#i=i+1;
					g.rate[iter, ind]=tmpd;
				}
			}
			else if(mscmd[i]=="-eN"){
				tmpd=as.numeric(mscmd[i+1]);i=i+1;
				if(fmatch(x=tmpd,y=time[iter])){
#print(mscmd[i-1]);
					tmpd=as.numeric(mscmd[i+1]);#i=i+1;
					for(k in 1:total.pop.num){
						N[iter,k]=tmpd*N0;
						g.rate[iter,k]=0;
					}
				}
			}
			else if(mscmd[i]=="-en"){
#print(mscmd[i:(i+3)])
				tmpd=as.numeric(mscmd[i+1]);i=i+1;
				if(fmatch(x=tmpd,y=time[iter])){
#print(mscmd[i-1]);
					ind=as.integer(mscmd[i+1]);i=i+1;
					tmpd=as.numeric(mscmd[i+1]);#i=i+1;
					N[iter, ind]=tmpd*N0;
					g.rate[iter, ind]=0;
				}
			}
		}

	}
	out=list(N=N, g.rate=g.rate);
	#print(N)
	return(out);
}

readm<-function(demograph_out)
{
	mscmd=demograph_out$mscmd;
	time=demograph_out$time.series;
	total.pop.num=demograph_out$total.pop.num;
	present.pop.num=demograph_out$present.pop.num;
	N4=demograph_out$N4;

	L<-length(mscmd);
	time_num<-length(time);
	init.M=-9;
	for(i in 1:L){
		if(mscmd[i]=="-I"){
			if(i+2+present.pop.num<= L){
				if(!iscmd(mscmd[i+2+present.pop.num]))init.M=as.numeric(mscmd[i+2+present.pop.num]);
			}
		}
	}
	m<-array(-9, dim=c(total.pop.num,total.pop.num,time_num));
	for(k in 1:present.pop.num){
		for(l in 1:present.pop.num){
			if(init.M>0){
				m[k,l,1]=init.M;
			}			else {
				m[k,l,1]=0;
			}
		}
	}

	for(iter in 1:time_num){
		for(i in 1:L){
			if(mscmd[i]=="-m"){
				ind1=as.integer(mscmd[i+1]);i=i+1;
				ind2=as.integer(mscmd[i+1]);i=i+1;
				tmpd=as.numeric(mscmd[i+1]);i=i+1;
				m[ind1,ind2,1]=tmpd/N4;
			}
			else if(mscmd[i]=="-ma"){
				for(k in 1:present.pop.num){
					for(l in 1:present.pop.num){
						if(mscmd[i+1]=="x"){tmpd=0;}
						else {tmpd=as.numeric(mscmd[i+1]);}
						i=i+1;
						m[k,l,1]=tmpd/(present.pop.num-1)/N4;
					}
				}
			}

			else if(mscmd[i]=="-eM"){
				tmpd=as.numeric(mscmd[i+1]);i=i+1;
				if(fmatch(x=tmpd,y=time[iter])){
					tmpd=as.numeric(mscmd[i+1]);i=i+1;
					for(k in 1:total.pop.num){
						for(l in 1:total.pop.num){
							m[k,l,iter]=tmpd/(present.pop.num-1)/N4;
						}
						m[k,k,1]=0;
					}
				}
			}
			else if(mscmd[i]=="-em"){
				tmpd=as.numeric(mscmd[i+1]);i=i+1;
				if(fmatch(x=tmpd,y=time[iter])){
					ind1=as.integer(mscmd[i+1]);i=i+1;
					ind2=as.integer(mscmd[i+1]);i=i+1;
					m[ind1, ind2, iter]=as.numeric(mscmd[i+1])/N4;i=i+1;
				}
			}
			else if(mscmd[i]=="-ema"){
				tmpd=as.numeric(mscmd[i+1]);i=i+1;
				if(fmatch(x=tmpd,y=time[iter])){
					tmp_num=as.integer(mscmd[i+1]);i=i+1;
					for(k in 1:total.pop.num){
						for(l in 1:total.pop.num){
							if(mscmd[i+1]=="x"){tmpd=0;}
							else {tmpd=as.numeric(mscmd[i+1]);}
							m[k,l,iter]=tmpd/N4;i=i+1;
						}
						m[k,k,1]=0;
					}
				}
			}


		}
	}
	return(m);
}

updatescrm<-function(demograph_out)
{
	m=demograph_out$m;
	mscmd=demograph_out$mscmd;
	time=demograph_out$time.series;
	L<-length(mscmd);
	time_num<-length(time);
	for(iter in 1:time_num){
		for(i in 1:L){
			if(mscmd[i]=="-eps"){
				tmpd=as.numeric(mscmd[i+1]);i=i+1;
				if(fmatch(x=tmpd,y=time[iter])){
					ind1=as.integer(mscmd[i+1]);i=i+1;
					ind2=as.integer(mscmd[i+1]);i=i+1;
					m[ind1, ind2, iter]=1-as.numeric(mscmd[i+1]);i=i+1;
				}
			}
		}
	}

	return(m);
}

readpos<-function(demograph_out)
{
	N0=1;
	time=demograph_out$time.series;
	total.pop.num=demograph_out$total.pop.num;
	present.pop.num=demograph_out$present.pop.num;
	mscmd=demograph_out$mscmd;
	inpos=demograph_out$pop.pos;
	g.rate=demograph_out$g.rate;
	m=demograph_out$m;
	N=demograph_out$N;


	L<-length(mscmd);
	time_num<-length(time);
	inPos<-vector(length=total.pop.num);
	if(length(inpos)<1){
		print("Use the default positions");
		for(i in 1:total.pop.num){inPos[i]=i;}
	} else {
		for(i in 1:total.pop.num){inPos[i]=inpos[i];}
	}
	Pos<-matrix(-9, nrow=time_num+1, ncol=total.pop.num);
	Pos_map<-matrix(-9, nrow=time_num+1, ncol=total.pop.num);
	survive<-matrix(-9, nrow=2, ncol=total.pop.num);
	npop=present.pop.num;
	for(i in 1:present.pop.num){
		survive[1,i]=1;survive[2,i]=time_num;
		Pos[1,i]=inPos[i];Pos_map[1,i]=i;
	}
	for(iter in 1:time_num){
		if(iter>1){
			for(i in 1:total.pop.num){
				if(Pos[iter, i]<0)Pos[iter, i]=Pos[iter-1, i];
				if(Pos_map[iter, i]<0)Pos_map[iter, i]=Pos_map[iter-1, i];
			}
		}
		for(i in 1:L){
			if(mscmd[i]=="-es"){
				tmpd=as.numeric(mscmd[i+1]);i=i+1;
				if(fmatch(x=tmpd,y=time[iter])){
					ind=as.integer(mscmd[i+1]);i=i+1;
					tmpd=as.numeric(mscmd[i+1]);i=i+1;
					npop=npop+1;
					for(k in 1:iter-1){
						Pos[k, npop]=Pos[k, ind];
						Pos_map[k, npop]=Pos_map[k, ind];
					}
					Pos[iter, npop]=inPos[npop];
					Pos_map[iter, npop]=npop;
					survive[1,npop]=iter;
					g.rate[iter, npop]=0;
					for(k in 1:npop){
						m[k, npop, iter]=0;m[npop,k, iter]=0;
					}
					#print("there1")
					m[ind, npop, iter]=1;
					if(iter<time_num){if(m[ind, npop, iter+1]<0)m[ind, npop, iter+1]=0;}
					N[iter, npop]=N0;
				}
			}
			else if(mscmd[i]=="-eA"){
				if(FALSE){
					tmpd=as.numeric(mscmd[i+1]);i=i+1;
					if(fmatch(x=tmpd,y=time[iter])){
						ind=as.integer(mscmd[i+1]);i=i+1;
						tmpd=as.numeric(mscmd[i+1]);i=i+1;
						npop=npop+1;
						for(k in 1:iter-1){
							Pos[k, npop]=Pos[k, ind];
							Pos_map[k, npop]=Pos_map[k, ind];
						}
						Pos[iter, npop]=inPos[npop];
						Pos_map[iter, npop]=npop;
						survive[1,npop]=iter;
						survive[2,npop]=iter;
						g.rate[iter, npop]=0;
						for(k in 1:npop){
							m[k, npop, iter]=0;m[npop,k, iter]=0;
						}
						#print("there1")
						m[npop, ind, iter]=1;
						if(iter<time_num){if(m[npop, ind, iter+1]<0)m[npop, ind, iter+1]=0;}
						N[iter, npop]=N0;
					}
				}
			}
			else if(mscmd[i]=="-ej"){
				tmpd=as.numeric(mscmd[i+1]);i=i+1;
				if(fmatch(x=tmpd,y=time[iter])){
					ind1=as.integer(mscmd[i+1]);i=i+1;
					ind2=as.integer(mscmd[i+1]);i=i+1;
					survive[2,ind1]=iter;
					for(k in 1:npop){
						m[ind1, k, iter]=0;
					}
					m[ind1, ind2, iter]=1;
					if(iter<time_num){if(m[ind1, ind2, iter+1]<0)m[ind1, ind2, iter+1]=0;}
					#g.rate[iter+1, ind1]=0;
					Pos[iter+1, ind1]=Pos[iter, ind2];
					Pos_map[iter+1, ind1]=Pos_map[iter, ind2];
				}
			}
		}
	}
	for(k in 1:total.pop.num){
		if(Pos[time_num+1, k]<0)Pos[time_num+1, k]=Pos[time_num, k];
		if(Pos_map[time_num+1, k]<0)Pos_map[time_num+1, k]=Pos_map[time_num, k];
	}

	out<-list(time.series=time, survive=survive, Pos=Pos, g.rate=g.rate, m=m, N=N, Pos_map=Pos_map);
	return(out);
}

update<-function(demograph_out)
{
	total.pop.num=demograph_out$total.pop.num;
	time.series=demograph_out$time.series;
	g.rate=demograph_out$g.rate;
	m=demograph_out$m;
	N=demograph_out$N;
	L=length(time.series);
	survive=demograph_out$survive;	
	if(L==1){
		for(i in 1:total.pop.num){
			if(N[1,i]<0)N[1,i]=1;
			if(g.rate[1,i]<0)g.rate[1,i]=0;
		}
	}
	if(L>1){
		for(i in 2:L){
			for(k in 1:total.pop.num){
				#update population growth rate and Ne
				if(g.rate[i, k]<0)g.rate[i, k]=g.rate[i-1, k];
				if(N[i,k]<0&N[i-1,k]>0)N[i,k]=N[i-1, k]*exp(-g.rate[i-1, k]*(time.series[i]-time.series[i-1]));
				#update the m
				for(l in 1:total.pop.num){
					if(m[k,l,i]<0&m[k,l,i-1]>-1)m[k,l,i]=m[k,l,i-1];
					if(i>survive[2,k]|i>survive[2,l])m[k,l,i]=0;
				}
			}
		}
	}


	out<-list(g.rate=g.rate, m=m, N=N);
	return(out);
}
