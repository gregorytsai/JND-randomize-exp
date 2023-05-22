##Random Experiment process
#Copyright: GregTsai 2010 
#gregorytsai@gmail.com
#���w�D���]�l�ƥءA�ƦC����]�l�A�H���ƦC
#�٨S�gblock.level=0����

rm(list=ls())

##########################################################
#################### ----��ƿ�J��---- ####################
#input(a.level)(b.level)(c.level)(exp.rpt)(random.method)(a)(b) 
#(block.level)

id="b98207078"

##��J�]�l��(�D�]�l,���]�l�έn���Ū��ƥ�)
a.level=c(3,9)
b.level=c(2,2)
c.level=c(0,0)

##��J�]�l�W��
a.name=c("2cm","8cm","32cm")
b.name=c("����","����")
c.name=c(NA,NA)
asub=c("����1","����2","����3","����4","����5","����6","����7","����8","����9")
bsub=c("�зǨ�E��m","�A�˨�E��m")
csub=c(NA,NA)

##��J���筫�Ʀ���
exp.rpt=4

##���w�H���覡(0,1,2):0.�L�H�� 1.�����H���k 2.�����H���k
random.method=2 

##���w�����H�����϶���� 0:6
block.level=4


#################### ----��J�ϵ���---- ####################
##########################################################
##########################################################
#################### ----�{���϶}�l---- ####################


####�ƦC�Ҧ�����զX####
#used variables:(a.level)(b.level)(c.level)(exp.rpt)
#new variables:(col)(trial.num)(odr.all)(odr.level)(odr.matrix)
#(odr.list)(col.num)(times)(trial.matrix)(onetrial.matrix)

col=6 #��K�p���`�@���禸��(trial.num�����ۭ�)�A�ó]�w�������
if (c.level[2]==0 || c.level[2]==1) {c.level[2]=1; col=col-1} 
if (c.level[1]==0) {c.level[1]=1; col=col-1} 
if (b.level[2]==0 || b.level[2]==1) {b.level[2]=1; col=col-1}
if (b.level[1]==0) {b.level[1]=1; col=col-1}
if (a.level[2]==0 || a.level[2]==1) {a.level[2]=1; col=col-1}

#�X�֩Ҧ�LEVEL�����@�x�}
level.matrix=matrix(NA, nrow=3, ncol=2) 
level.matrix[1,]=a.level
level.matrix[2,]=b.level
level.matrix[3,]=c.level

trial.num=exp.rpt
for (i in 1:6) (trial.num=trial.num*level.matrix[i])
trial.matrix=matrix(data=NA, nrow=trial.num, ncol=col)
onetrial.matrix=matrix(NA, nrow=(trial.num/exp.rpt), ncol=col)


#�����H���A�̫�@��A�ƦC�H���]�l�϶���K�\Ū
if (random.method==0) (block.level=6)
if (block.level==0) (block.level=6)

{
	odr.all=c(1:6)             #�Ψӭp��ƦC����
	for (i in block.level:6)   #��block.level���̫�
	{
		if (i==6) (odr.all[i]=block.level)
		else (odr.all[i]=odr.all[i+1])
	}

	odr.list=rep(NA,col)       #��z���禸�ǹ����M��
	col.num=1
	for (i in 1:6)
	{
		if (level.matrix[odr.all[i]]!=1) {odr.list[col.num]=odr.all[i]; col.num=col.num+1}
	}


	odr.matrix=matrix(NA, nrow=3, ncol=2 )
	for (i in 1:6) (odr.matrix[i]=level.matrix[odr.all[i]])


	odr.level=rep(NA,col)      #��z���禸��level�ƥ�
	col.num=1
	for (i in 1:6)
	{
		if (odr.matrix[i]!=1) {odr.level[col.num]=odr.matrix[i]; col.num=col.num+1}
	}

	
	for (i in 1:col)
	{
		times.post=1
		if (i!=col) {for (t in (i+1):col) (times.post=times.post*odr.level[t])} 
		times.pre=1
		if (i!=1) {for (t in 1:(i-1)) (times.pre=times.pre*odr.level[t])}
		for (j in 1:times.pre)
		{
			for (k in 1:odr.level[i])
			{
				for (l in 1:times.post)
				{
					onetrial.matrix[((j-1)*times.post*odr.level[i]+(k-1)*times.post+l),i]=k
				}
			}
		}
	}
}


#���ƹ��禸��


for (i in 1:exp.rpt)
{
	for (j in 1:(trial.num/exp.rpt))
	{
		trial.matrix[((i-1)*(trial.num/exp.rpt)+j),]=onetrial.matrix[j,]
	}
}





####�B���H���x�}####
#used variables:(odr.level)(exp.rpt)(trial.num)(col)
#new variables:(random.matrix)(rnd.times)(rnd.other)

rnd.times=exp.rpt
for (i in 1:(col-1)) (rnd.times=rnd.times*odr.level[i])
random.other=sample(1:rnd.times, size=rnd.times, replace=F)


if (block.level!=0)
{
	random.matrix=matrix(1:trial.num, ncol=1)
	for (i in 1:rnd.times)
	{
		random.block=sample(1:(odr.level[col]), size=(odr.level[col]), replace=F)
		for (j in 1:(odr.level[col]))
		{
			random.matrix[((i-1)*odr.level[col]+j),1]=((random.other[i]-1)*odr.level[col]+random.block[j])  
		}
	}
}

if (random.method==1)
{
	random.matrix=matrix(sample(1:trial.num,size=trial.num),ncol=1)
}


save.image(file="temp.Rdata")



####�}�l����{�ǡA�������####
#used variables:(trial.num)(random.matrix)(trial.matrix)(odr.list)
#new variables:(id)(result.matrix)(name)

#�������G�x�}�A�������H���ƦC

name=list(a=a.name,b=b.name,c=c.name,as=asub,bs=bsub,cs=csub)

start=function(random.matrix,trial.matrix,odr.list,name)
{
	save.image(file="temp.Rdata")

	result.matrix=matrix(NA, nrow=trial.num, ncol=1)
	
	attach(name)
	x=random.matrix[1]
	id=""
	for (i in 1:trial.num)
	{
		x.last=x           #�p�G���~��K�ק�e�@��
		x=random.matrix[i] #�̷��H�����i�����
		if ((i-1)%%odr.level[col]==0) (cat("[New Group]"))
		cat("����", i,"�A")
		for (j in 1:col)   #��ܹ���]�l���e
		{
			switch(odr.list[j],"1"=cat(a[trial.matrix[x,j]],"�A"),
			                   "2"=cat(b[trial.matrix[x,j]],"�A"), 
			                   "3"=cat(c[trial.matrix[x,j]],"�A"),
			                   "4"=cat(as[trial.matrix[x,j]],"�A"),
			                   "5"=cat(bs[trial.matrix[x,j]],"�A"),
			                   "6"=cat(cs[trial.matrix[x,j]],"�A"), )                   
		}
		ans=readline("(a���,r�ק�e�@��):")  #��������
		if (ans=="a") (break)    #quitter!
		else if (ans=="r")       #�ק�e�@�Ӹ��
		{
	
			cat("[�ץ�]����", i-1,"�A")
			for (j in 1:col)
			{
				switch(odr.list[j],"1"=cat(a[trial.matrix[x.last,j]],"�A"),
			                  	 "2"=cat(b[trial.matrix[x.last,j]],"�A"), 
				                   "3"=cat(c[trial.matrix[x.last,j]],"�A"),
				                   "4"=cat(as[trial.matrix[x.last,j]],"�A"),
				                   "5"=cat(bs[trial.matrix[x.last,j]],"�A"),
				                   "6"=cat(cs[trial.matrix[x.last,j]],"�A"), )                   
			}
			ans=readline("�Э��s��J:")
			result.matrix[x.last]=ans
			cat("����", i,"�A")
			for (j in 1:col)
			{
				switch(odr.list[j],"1"=cat(a[trial.matrix[x,j]],"�A"),
				                   "2"=cat(b[trial.matrix[x,j]],"�A"), 
				                   "3"=cat(c[trial.matrix[x,j]],"�A"),
				                   "4"=cat(as[trial.matrix[x,j]],"�A"),
				                   "5"=cat(bs[trial.matrix[x,j]],"�A"),
				                   "6"=cat(cs[trial.matrix[x,j]],"�A"), )                   
			}
			ans=readline(":")
		}
		result.matrix[x]=ans
		#�w�������C���x�s���絲�G�Ȧs��
		write.table(result.matrix, file="exptemp.txt", col.name=F, row.name=F, quote=F)
	}
	detach(name)
	cat("\n")
	cat("---���絲��---","\n")
	return(result.matrix)	
}




result=start(random.matrix,trial.matrix,odr.list,name)

result.matrix=matrix(as.integer(result), ncol=1)



#####�N���G�x�s��TXT�ɮ�####
#used variables:(result.matrix)(random.matrix)...
#new variables:(exp2.txt)(nrow)(exp.matrix)...

#�@�����榡
exp.matrix=cbind(random.matrix,trial.matrix,result.matrix)

#���ন�߹����G�榡
nrow=(trial.num/odr.level[col])*2
exp2.txt=matrix(NA,ncol=(col+odr.level[col]),nrow=nrow)  

#�аO���禸��
for (i in 1:(nrow/2)) (exp2.txt[(2*i-1)]=i) 

#�аO��l����
for (i in 1:trial.num) 
{
	#�p��s�榡��row
	exp2.row=((random.matrix[i]-1)%/%odr.level[col])+1
	k=2*exp2.row-1
	#�p��s�榡��column
	exp2.col=((random.matrix[i]-1)%%odr.level[col])+1
	#�аO����e�T�]�l
	for (j in 1:(ncol(trial.matrix)-1))
	{	
		exp2.txt[k,1+j]=trial.matrix[i,j]
	}
	#�аO��V���̫�B�z���]�l
	exp2.txt[k,(ncol(trial.matrix)+exp2.col)]=trial.matrix[i,col]
	#�аO���G
	exp2.txt[k+1,(ncol(trial.matrix)+exp2.col)]=result.matrix[i]
}


#�����z��X�ɮ�
exp2.fmt=matrix(NA,ncol=(col+odr.level[col]),nrow=nrow)
exp2.fmt=exp2.txt
exp2.fmt[,1]=format(exp2.txt[,1], digits=2)
for (i in 1:nrow(exp2.fmt)) (if (exp2.fmt[i,1]=="NA") (exp2.fmt[i,1]="  "))

#�x�s�߹�榡���絲�G�A�w�]��b�ڪ����
fmt=".txt"
filename=paste(id,fmt, sep="")
write.table(exp2.fmt,file=filename,na=" ", col.name=F, row.name=F, quote=F)

#�x�s���`�榡���絲�G�A�w�]��b�ڪ����
fmt=".txt"
middle="_nornal"
filename=paste(id,middle,fmt, sep="")
exp.matrix=exp.matrix[unlist(sapply(sort(unique(exp.matrix[,1])), function(x) which(exp.matrix[,1]==x)[order(exp.matrix[(exp.matrix[,1]==x),2], decreasing=TRUE)])),]

exp.fmt=matrix(NA,ncol=(col+odr.level[col]),nrow=nrow)
exp.fmt=exp.matrix
exp.fmt[,1]=format(exp.matrix[,1], digits=2)
for (i in 1:nrow(exp.fmt)) (if (exp.fmt[i,1]=="NA") (exp.fmt[i,1]="  "))

write.table(exp.matrix,file=filename,na=" ", col.name=F, row.name=F, quote=F)





