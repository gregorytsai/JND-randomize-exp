##Random Experiment process
#Copyright: GregTsai 2010 
#gregorytsai@gmail.com
#指定主次因子數目，排列實驗因子，隨機排列
#還沒寫block.level=0怎麼排

rm(list=ls())

##########################################################
#################### ----資料輸入區---- ####################
#input(a.level)(b.level)(c.level)(exp.rpt)(random.method)(a)(b) 
#(block.level)

id="b98207078"

##輸入因子數(主因子,次因子或要平衡的數目)
a.level=c(3,9)
b.level=c(2,2)
c.level=c(0,0)

##輸入因子名稱
a.name=c("2cm","8cm","32cm")
b.name=c("水平","垂直")
c.name=c(NA,NA)
asub=c("長度1","長度2","長度3","長度4","長度5","長度6","長度7","長度8","長度9")
bsub=c("標準刺激位置","顛倒刺激位置")
csub=c(NA,NA)

##輸入實驗重複次數
exp.rpt=4

##指定隨機方式(0,1,2):0.無隨機 1.完全隨機法 2.分區隨機法
random.method=2 

##指定分區隨機的區塊單位 0:6
block.level=4


#################### ----輸入區結束---- ####################
##########################################################
##########################################################
#################### ----程式區開始---- ####################


####排列所有實驗組合####
#used variables:(a.level)(b.level)(c.level)(exp.rpt)
#new variables:(col)(trial.num)(odr.all)(odr.level)(odr.matrix)
#(odr.list)(col.num)(times)(trial.matrix)(onetrial.matrix)

col=6 #方便計算總共實驗次數(trial.num直接相乘)，並設定資料欄位數
if (c.level[2]==0 || c.level[2]==1) {c.level[2]=1; col=col-1} 
if (c.level[1]==0) {c.level[1]=1; col=col-1} 
if (b.level[2]==0 || b.level[2]==1) {b.level[2]=1; col=col-1}
if (b.level[1]==0) {b.level[1]=1; col=col-1}
if (a.level[2]==0 || a.level[2]==1) {a.level[2]=1; col=col-1}

#合併所有LEVEL成為一矩陣
level.matrix=matrix(NA, nrow=3, ncol=2) 
level.matrix[1,]=a.level
level.matrix[2,]=b.level
level.matrix[3,]=c.level

trial.num=exp.rpt
for (i in 1:6) (trial.num=trial.num*level.matrix[i])
trial.matrix=matrix(data=NA, nrow=trial.num, ncol=col)
onetrial.matrix=matrix(NA, nrow=(trial.num/exp.rpt), ncol=col)


#分區隨機，最後一欄再排列隨機因子區塊方便閱讀
if (random.method==0) (block.level=6)
if (block.level==0) (block.level=6)

{
	odr.all=c(1:6)             #用來計算排列順序
	for (i in block.level:6)   #把block.level丟到最後
	{
		if (i==6) (odr.all[i]=block.level)
		else (odr.all[i]=odr.all[i+1])
	}

	odr.list=rep(NA,col)       #整理實驗次序對應清單
	col.num=1
	for (i in 1:6)
	{
		if (level.matrix[odr.all[i]]!=1) {odr.list[col.num]=odr.all[i]; col.num=col.num+1}
	}


	odr.matrix=matrix(NA, nrow=3, ncol=2 )
	for (i in 1:6) (odr.matrix[i]=level.matrix[odr.all[i]])


	odr.level=rep(NA,col)      #整理實驗次序level數目
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


#重複實驗次數


for (i in 1:exp.rpt)
{
	for (j in 1:(trial.num/exp.rpt))
	{
		trial.matrix[((i-1)*(trial.num/exp.rpt)+j),]=onetrial.matrix[j,]
	}
}





####運算隨機矩陣####
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



####開始實驗程序，紀錄資料####
#used variables:(trial.num)(random.matrix)(trial.matrix)(odr.list)
#new variables:(id)(result.matrix)(name)

#紀錄結果矩陣，不按照隨機排列

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
		x.last=x           #如果有誤方便修改前一個
		x=random.matrix[i] #依照隨機表進行實驗
		if ((i-1)%%odr.level[col]==0) (cat("[New Group]"))
		cat("實驗", i,"，")
		for (j in 1:col)   #顯示實驗因子內容
		{
			switch(odr.list[j],"1"=cat(a[trial.matrix[x,j]],"，"),
			                   "2"=cat(b[trial.matrix[x,j]],"，"), 
			                   "3"=cat(c[trial.matrix[x,j]],"，"),
			                   "4"=cat(as[trial.matrix[x,j]],"，"),
			                   "5"=cat(bs[trial.matrix[x,j]],"，"),
			                   "6"=cat(cs[trial.matrix[x,j]],"，"), )                   
		}
		ans=readline("(a放棄,r修改前一個):")  #紀錄答案
		if (ans=="a") (break)    #quitter!
		else if (ans=="r")       #修改前一個資料
		{
	
			cat("[修正]實驗", i-1,"，")
			for (j in 1:col)
			{
				switch(odr.list[j],"1"=cat(a[trial.matrix[x.last,j]],"，"),
			                  	 "2"=cat(b[trial.matrix[x.last,j]],"，"), 
				                   "3"=cat(c[trial.matrix[x.last,j]],"，"),
				                   "4"=cat(as[trial.matrix[x.last,j]],"，"),
				                   "5"=cat(bs[trial.matrix[x.last,j]],"，"),
				                   "6"=cat(cs[trial.matrix[x.last,j]],"，"), )                   
			}
			ans=readline("請重新輸入:")
			result.matrix[x.last]=ans
			cat("實驗", i,"，")
			for (j in 1:col)
			{
				switch(odr.list[j],"1"=cat(a[trial.matrix[x,j]],"，"),
				                   "2"=cat(b[trial.matrix[x,j]],"，"), 
				                   "3"=cat(c[trial.matrix[x,j]],"，"),
				                   "4"=cat(as[trial.matrix[x,j]],"，"),
				                   "5"=cat(bs[trial.matrix[x,j]],"，"),
				                   "6"=cat(cs[trial.matrix[x,j]],"，"), )                   
			}
			ans=readline(":")
		}
		result.matrix[x]=ans
		#預防當機每次儲存實驗結果暫存檔
		write.table(result.matrix, file="exptemp.txt", col.name=F, row.name=F, quote=F)
	}
	detach(name)
	cat("\n")
	cat("---實驗結束---","\n")
	return(result.matrix)	
}




result=start(random.matrix,trial.matrix,odr.list,name)

result.matrix=matrix(as.integer(result), ncol=1)



#####將結果儲存成TXT檔案####
#used variables:(result.matrix)(random.matrix)...
#new variables:(exp2.txt)(nrow)(exp.matrix)...

#一般實驗格式
exp.matrix=cbind(random.matrix,trial.matrix,result.matrix)

#先轉成心實實驗二格式
nrow=(trial.num/odr.level[col])*2
exp2.txt=matrix(NA,ncol=(col+odr.level[col]),nrow=nrow)  

#標記實驗次序
for (i in 1:(nrow/2)) (exp2.txt[(2*i-1)]=i) 

#標記其餘項目
for (i in 1:trial.num) 
{
	#計算新格式的row
	exp2.row=((random.matrix[i]-1)%/%odr.level[col])+1
	k=2*exp2.row-1
	#計算新格式的column
	exp2.col=((random.matrix[i]-1)%%odr.level[col])+1
	#標記實驗前三因子
	for (j in 1:(ncol(trial.matrix)-1))
	{	
		exp2.txt[k,1+j]=trial.matrix[i,j]
	}
	#標記橫向的最後處理次因子
	exp2.txt[k,(ncol(trial.matrix)+exp2.col)]=trial.matrix[i,col]
	#標記結果
	exp2.txt[k+1,(ncol(trial.matrix)+exp2.col)]=result.matrix[i]
}


#對齊整理輸出檔案
exp2.fmt=matrix(NA,ncol=(col+odr.level[col]),nrow=nrow)
exp2.fmt=exp2.txt
exp2.fmt[,1]=format(exp2.txt[,1], digits=2)
for (i in 1:nrow(exp2.fmt)) (if (exp2.fmt[i,1]=="NA") (exp2.fmt[i,1]="  "))

#儲存心實格式實驗結果，預設放在我的文件
fmt=".txt"
filename=paste(id,fmt, sep="")
write.table(exp2.fmt,file=filename,na=" ", col.name=F, row.name=F, quote=F)

#儲存正常格式實驗結果，預設放在我的文件
fmt=".txt"
middle="_nornal"
filename=paste(id,middle,fmt, sep="")
exp.matrix=exp.matrix[unlist(sapply(sort(unique(exp.matrix[,1])), function(x) which(exp.matrix[,1]==x)[order(exp.matrix[(exp.matrix[,1]==x),2], decreasing=TRUE)])),]

exp.fmt=matrix(NA,ncol=(col+odr.level[col]),nrow=nrow)
exp.fmt=exp.matrix
exp.fmt[,1]=format(exp.matrix[,1], digits=2)
for (i in 1:nrow(exp.fmt)) (if (exp.fmt[i,1]=="NA") (exp.fmt[i,1]="  "))

write.table(exp.matrix,file=filename,na=" ", col.name=F, row.name=F, quote=F)






