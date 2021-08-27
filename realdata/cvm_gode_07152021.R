library(np)
library(nlsrk)
library(splines2)
library(grpreg)
library(mvtnorm)
library(igraph)
library(ggplot2)
setwd("/Users/chixiangchen/Desktop/phd_in_psu/research/Prof._Wu/Network/tissue/blood_vessel/")




#t=seq(0,20,length=n)
#t=round(t,6)
set.seed(4351)
t<-exp_index
data_observe<-xsis
alpha<-1
#data_observe<-cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,(all_nodes[,12:pn]+0.01*t))
data_observe<-data_observe[c(1,diff(t,lag=1))>=0.0001,]
dim(data_observe)
x_cov<-as.numeric(index_data$smoke)
x_cov<-x_cov[c(1,diff(t,lag=1))>=0.0001]
t<-t[c(1,diff(t,lag=1))>=0.0001]
#n<-nrow(data_observe)
#plot(data_observe[,9])

#get fitted values of observations using varying coefficient kernel regression
data_fitted<-rep()
data_fitted_cov<-rep()
#bws<-0.01
for(i in 1:ncol(data_observe))
{
    bw <- npscoefbw(formula=data_observe[,i]~x_cov|t,betas=T,regtype="ll")
    bw <- npscoef(bw,betas=T,exdat=data.frame(x_cov=x_cov),ezdat=data.frame(t=t),regtype="ll")
    fitted<-coef(bw)[,1]
    data_fitted<-cbind(data_fitted,fitted)
    fitted_cov<-coef(bw)[,2]
    data_fitted_cov<-cbind(data_fitted_cov,fitted_cov)
}
id_sub<-8
par(mfrow = c(1, 2))
plot(x=t[x_cov==1],y=data_observe[x_cov==1,id_sub],col="red")
lines(data_fitted[,id_sub]+data_fitted_cov[,id_sub],x=t,col="red")
plot(x=t[x_cov==0],y=data_observe[x_cov==0,id_sub],col="blue")
lines(data_fitted[,id_sub],x=t,col="blue")

#plot(y=mean_node[,7],x=t)
#lines(data_fitted[,id_sub]+x_cov*data_fitted_cov[,id_sub],x=t,col="red")

#spline basis and its integration for baseline
library(splines2)
library(Matrix)
library(grpreg)
#some useful fct
step_fct<-function(x,step)
{
    t<-length(x)
    x_new<-x[1]
    for (i in 1:(t-1))
    {
        x_int<-seq(x[i],x[i+1],by=step)[-1]
        x_new<-c(x_new,x_int)
    }
    x_new
}

#generate B-spline and its integration for baseline
#generate B-spline and its integration for cov
X_big<-rep()
X_big_cov<-rep()
X_big_int<-rep()
X_big_int_cov<-rep()
degree<-3
round_num<-attr(regexpr("(?<=\\.)0+", format(min(diff(t,lag=1)),scientific = FALSE), perl = TRUE), "match.length")+1
t<-round(t,round_num)
if(min(diff(t,lag=1))>=2*10^(-round_num))
{
    step<-10^(-round_num)
}else
{
    step<-5*10^(-round_num-1)
    t<-round(t,round_num+1)
}
cluster<-ncol(data_observe)
lambda<-seq()
for (i in 1:cluster)
{
    knots <- sort(data_fitted[,i])[c(round(nrow(data_fitted)/4),round(nrow(data_fitted)*2/4),round(nrow(data_fitted)*3/4))]
    knots_cov <-sort(data_fitted_cov[,i])[c(round(nrow(data_fitted)/4),round(nrow(data_fitted)*2/4),round(nrow(data_fitted)*3/4))]
    bsMat <- bSpline(data_fitted[,i], knots = knots, degree = degree, intercept=F)
    bsMat_cov <- bSpline(data_fitted_cov[,i], knots = knots_cov, degree = degree, intercept=F)
    x_int<-step_fct(as.vector(t),step)
    if(min(diff(t,lag=1))<2*10^(-round_num))
    {x_int<-round(x_int,round_num+1)}else
    {x_int<-round(x_int,round_num)}
    x_cov_int<-rep(0,length(x_int))
    bw <- npscoefbw(formula=data_observe[,i]~x_cov|t,betas=T,regtype="ll")
    bw <- npscoef(bw,betas=T,exdat=data.frame(x_cov=c(x_cov_int)),ezdat=data.frame(t=c(x_int)),regtype="ll")
    x_hat<-coef(bw)[,1]
    x_hat_cov<-coef(bw)[,2]
    
    #x_hat<-exp(cbind(1,log(x_int)) %*% beta[,i])
    basis_int<-bSpline(x_hat, knots = knots, degree = degree, intercept=F)
    basis_int_cov<-bSpline(x_hat_cov, knots = knots_cov, degree = degree, intercept=F)
    base_int<-0
    base_int_cov<-0
    for (j in 1:(length(t)-1))
    {
        int_row<-apply(basis_int[x_int>=t[j] & x_int<=t[j+1],][-1,],2,function (x) sum(x)*step)
        base_int<-rbind(base_int,int_row)
        
        int_row_cov<-apply(basis_int_cov[x_int>=t[j] & x_int<=t[j+1],][-1,],2,function (x) sum(x)*step)
        base_int_cov<-rbind(base_int_cov,int_row_cov)
    }
    base_int<-apply(base_int,2,cumsum)
    base_int_cov<-apply(base_int_cov,2,cumsum)
    #base_int<-ibs(data_fitted[,i], knots = knots, degree = degree, intercept=F)*(max(t)-min(t)) #not correct
    X_big <- cbind(X_big,bsMat)
    X_big_cov <- cbind(X_big_cov,bsMat_cov)
    
    X_big_int<-cbind(X_big_int,base_int)
    X_big_int_cov<-cbind(X_big_int_cov,base_int_cov)
}
dim(X_big)
dim(X_big_int)

dim(X_big_cov)
dim(X_big_int_cov)

#X_big<-apply(X_big,2,function(x) x/sum(x))
#X_big_int<-apply(X_big_int,2,function(x) x/sum(x))
num_cov<-degree+length(knots)
#X_big_int_exp_intcep<-cbind(1,t,X_big_int,x_cov)

#some setup for final analysis
X_big_int_exp<-cbind(t,X_big_int,x_cov*X_big_int_cov,t*x_cov,x_cov)
X_big_int_exp_intcep<-cbind(1,t,X_big_int,x_cov*X_big_int_cov,t*x_cov,x_cov)
X_big_int_exp_intcep_0<-cbind(1,t,X_big_int,0*X_big_int_cov,t*0,0)
X_big_int_exp_intcep_1<-cbind(1,t,X_big_int,1*X_big_int_cov,t*1,1)

x_cov<-matrix(x_cov,ncol=1)
group<-c(0,rep(1:(2*ncol(data_observe)),each=num_cov),rep(0,ncol(x_cov)),rep(0,ncol(x_cov)))
self_size<-rep()
self_size_cov<-rep()
self_size_all<-rep()
gene_whole<-rep()
gene_whole_cov<-rep()
gene_whole_all<-rep()
fitted_1_all<-rep()
fitted_0_all<-rep()
smoking_effect_all<-rep()
trend_base_all<-rep()
trend_smoking_effect_all<-rep()

for (j in 1:ncol(data_observe))  #
{
    #cvfit <- cv.grpreg(X_big_int_exp, data_observe[,j], group=group, penalty="grLasso",nfolds=5,alpha=alpha1)
    #beta_ini<-grpreg(X_big_int_exp, data_observe[,j], group=group, penalty="grLasso",lambda = cvfit$lambda.min,alpha=alpha1)$beta
    #group_id<-c(1,2,rep(3:(2+ncol(data_observe)),each=num_cov))
    #beta_each_g<-unlist(tapply(beta_ini,group_id,function(x) x))
    #weight_ini<-rep()
    #for(i in 1:ncol(data_observe))
    #{
    #  m<-sum((beta_each_g[(3+(i-1)*num_cov):(3+(i-1)*num_cov+num_cov-1)])^2)
    #  weight_ini<-c(weight_ini,m)
    #}
    #weight<-(1/(weight_ini+0.0000001))^0.3
    cvfit <- cv.grpreg(X_big_int_exp, data_observe[,j], group=group, penalty="grLasso",nfolds=20,alpha=alpha,max.iter=100000,seed=23556)
    fit<-grpreg(X_big_int_exp, data_observe[,j], group=group, penalty="grLasso",lambda =cvfit$lambda.min,alpha=alpha,max.iter=100000)
    #beta_group<-fit$beta
    #fit<-grpreg(X_big_int_exp, data_observe[,j], group=group, penalty="grLasso",alpha=alpha,lambda.min = lambda.min)
    #fit<-grpreg(X_big_int_exp, data_observe[,j], group=group, penalty="grLasso",alpha=alpha,lambda = cvfit$lambda.min)
    #beta_group<-select(fit,"EBIC")$beta
    #cvfit <- cv.grpreg(X_big_int_exp, data_observe[,4], group=group, penalty="grLasso",nfolds=50)
    #fit<-grpreg(X_big_int_exp, data_observe[,4], group=group, penalty="grLasso",lambda = cvfit$lambda.min)
    #fit<-grpreg(X_big_int_exp, data_observe[,2], group=group, penalty="grLasso")
    #beta_group<-select(fit,"EBIC")$beta
    #weight
    
    beta_select<-fit$beta
    index_nonzero<-which(beta_select!=0)
    beta_select[index_nonzero]
    
    find_index<-index_nonzero[-c(1:2,length(index_nonzero)-1,length(index_nonzero))]-2
    find_index_comb<-unique(ifelse(find_index>num_cov*cluster,find_index-num_cov*cluster,find_index))
    #influ<-find_index[ which(find_index%%num_cov==0)]/num_cov
    #m<-which(influ==j)
    #if (length(m)==0)
    #{
    fitted<-X_big_int_exp_intcep[,c(index_nonzero)]%*%beta_select[index_nonzero]
    fitted_1<-X_big_int_exp_intcep_1[,c(index_nonzero)]%*%beta_select[index_nonzero]
    fitted_0<-X_big_int_exp_intcep_0[,c(index_nonzero)]%*%beta_select[index_nonzero]
    fitted_1_all<-cbind(fitted_1_all,fitted_1)
    fitted_0_all<-cbind(fitted_0_all,fitted_0)
    fitted_self<-X_big_int_exp_intcep[,1:2]%*%beta_select[1:2]
    fitted_self<-t(fitted_self)
    self_size<-rbind(self_size,fitted_self)
    smoking_effect<-(X_big_int_exp_intcep_1[,ncol(X_big_int_exp_intcep_1)]*beta_select[ncol(X_big_int_exp_intcep_1)])[1]
    smoking_effect_all<-c(smoking_effect_all,smoking_effect)
    trend_base<-X_big_int_exp_intcep_0[,2]*beta_select[2]
    trend_base_all<-cbind(trend_base_all,trend_base)
    trend_smoking_effect<-X_big_int_exp_intcep_1[,ncol(X_big_int_exp_intcep_1)-1]*beta_select[ncol(X_big_int_exp_intcep_1)-1]
    trend_smoking_effect_all<-cbind(trend_smoking_effect_all,trend_smoking_effect)
    
    fitted_self_cov<-X_big_int_exp_intcep_1[,c(length(beta_select)-1,length(beta_select))]%*%beta_select[c(length(beta_select)-1,length(beta_select))]
    fitted_self_cov<-t(fitted_self_cov)
    self_size_cov<-rbind(self_size_cov,fitted_self_cov)
    
    fitted_self_all<-X_big_int_exp_intcep_1[,c(1,2,length(beta_select)-1,length(beta_select))]%*%beta_select[c(1,2,length(beta_select)-1,length(beta_select))]
    fitted_self_all<-t(fitted_self_all)
    self_size_all<-rbind(self_size_all,fitted_self_all)
    
    
    #}
    #else
    #{
    #  fitted<-X_big_int_exp_intcep[,c(index_nonzero)]%*%beta_select[index_nonzero]
    #  fitted_self<-X_big_int_exp_intcep[,c(1:2,(2+(j-1)*num_cov+1):(2+(j-1)*num_cov+num_cov))]%*%beta_select[c(1:2,(2+(j-1)*num_cov+1):(2+(j-1)*num_cov+num_cov))]
    #  fitted_self<-t(fitted_self)
    #  self_size<-rbind(self_size,fitted_self)
    #  find_self<-which(influ==j)
    #  find_index<-find_index[-c((num_cov*(find_self-1)+1):(num_cov*(find_self-1)+num_cov))]
    #}
    par(mfrow = c(1, 1))
    gene_index<-j
    plot(x=t,y=data_observe[,gene_index],ylim=c(-1,max(data_observe[,gene_index])+0.2))
    lines(x=t,y=fitted_self,col="black")
    lines(x=t,y=fitted,col="red")
    
    num_gene<-length(find_index)/num_cov
    num_gene_comb<-length(find_index_comb)/num_cov
    gene_cluster<-rep()
    gene_cluster_cov<-rep()
    gene_cluster_all<-rep()
    fitted_gene_matrix<-rep()
    fitted_gene_matrix_cov<-rep()
    fitted_gene_matrix_all<-rep()
    if (num_gene>0)
    {
        for (i in 1:num_gene)
        {
            gene_relate<-find_index[(i-1)*num_cov+1]%/%num_cov+1
            if (gene_relate<=ncol(data_observe))
            {
                gene_cluster<-c(gene_cluster,gene_relate)
                fitted_gene<-X_big_int_exp_intcep[,index_nonzero[(3+(i-1)*num_cov):(3+i*num_cov-1)]]%*%beta_select[index_nonzero[(3+(i-1)*num_cov):(3+i*num_cov-1)]]
                fitted_gene_matrix<-cbind(fitted_gene_matrix,fitted_gene)
                lines(x=t,y=fitted_gene,ylab="y",col=i+2)
            }else
            {
                gene_cluster_cov<-c(gene_cluster_cov,gene_relate)
                fitted_gene_cov<-X_big_int_exp_intcep_1[,index_nonzero[(3+(i-1)*num_cov):(3+i*num_cov-1)]]%*%beta_select[index_nonzero[(3+(i-1)*num_cov):(3+i*num_cov-1)]]
                fitted_gene_matrix_cov<-cbind(fitted_gene_matrix_cov,fitted_gene_cov)
                lines(x=t,y=fitted_gene_cov,ylab="y",col=i+2)
            }
        }
        if(!is.null(fitted_gene_matrix))
        {
            gene_cluster<-cbind(gene_cluster,j,t(fitted_gene_matrix))
            gene_whole<-rbind(gene_whole,gene_cluster)
        }
        if(!is.null(fitted_gene_matrix_cov))
        {
            gene_cluster_cov<-cbind(gene_cluster_cov,j,t(fitted_gene_matrix_cov))
            gene_whole_cov<-rbind(gene_whole_cov,gene_cluster_cov)
        }
        for (i in 1:num_gene_comb)
        {
            gene_relate_comb<-find_index_comb[(i-1)*num_cov+1]%/%num_cov+1
            gene_cluster_all<-c(gene_cluster_all,gene_relate_comb)
            fitted_gene_all<-X_big_int_exp_intcep_1[,c(find_index_comb[(1+(i-1)*num_cov):(1+i*num_cov-1)]+2,find_index_comb[(1+(i-1)*num_cov):(1+i*num_cov-1)]+2+cluster*num_cov)]%*%beta_select[c(find_index_comb[(1+(i-1)*num_cov):(1+i*num_cov-1)]+2,find_index_comb[(1+(i-1)*num_cov):(1+i*num_cov-1)]+2+cluster*num_cov)]
            fitted_gene_matrix_all<-cbind(fitted_gene_matrix_all,fitted_gene_all)
        }
        gene_cluster_all<-cbind(gene_cluster_all,j,t(fitted_gene_matrix_all))
        gene_whole_all<-rbind(gene_whole_all,gene_cluster_all)
    }
    print(j)
    
    ###capture the beta_cov and overall 
    
}
self_size<-cbind(rep(1:cluster),self_size)
dim(self_size)
self_size_cov<-cbind(rep(1:cluster),self_size_cov)
dim(self_size_cov)
self_size_all<-cbind(rep(1:cluster),self_size_all)
dim(self_size_all)
gene_whole_cov[,1]<-gene_whole_cov[,1]-ncol(data_observe)


##fitting plot smoke or not together in one figure
library(scales)
par(mfrow = c(2, 2))
par(cex=1)
par(mar = c(3.5, 3.5, 1, 1)) #bottom, left, top, right
#index_fitting<-c(8,16,50,51,61,68)
index_fitting<-c(13,50,28,62)
geneselect60<-colnames(data_observe)
geneselect60[index_fitting]
geneselect60_sub<-c("TOX2","FMNL3","GLRB","CBLL1")
zz<-0
for (i in index_fitting)
{
    zz<-zz+1
    
    plot(y=data_observe[x_cov==0,i],x=t[x_cov==0],ylab="",
         xlab="",ylim = c(min(data_observe[,i]),max(data_observe[,i])),col = alpha("red",0.6))
    lines(x=t,y=data_fitted[,i],col="red",lwd=2)
    
    points(y=data_observe[x_cov==1,i],x=t[x_cov==1],col=alpha("blue",0.6))
    lines(x=t,y=data_fitted[,i]+data_fitted_cov[,i],col="blue",lwd=2)
    text(x=0.6, y=min(data_observe[,i]), pos=4, labels=geneselect60_sub[zz],cex=1.3)
    #text(x=0.6, y=min(data_observe[,i])+0.3, pos=4)
    
    if (zz==1)
    {
        title(ylab="Gene Expression (log scale)", line=2.2, cex.lab=1.4)
        #title(xlab="Hypertension Probability", line=2.2, cex.lab=1.2)
    }else if(zz==3){
        title(ylab="Gene Expression (log scale)", line=2.2, cex.lab=1.4)
        title(xlab="Hypertension Probability", line=2.2, cex.lab=1.4)
    }else if(zz==4){
        title(xlab="Hypertension Probability", line=2.2, cex.lab=1.4)
    }
    
}

##network for x=0 at i=250 after some tuning
par(mfrow = c(1, 1))
par(cex=1)
par(mar = c(0, 0, 0, 0), oma = c(6, 6, 3, 3))
tunned_value<-0.03
for (i in c(250))
{
    tunning_index<-apply(gene_whole_all,1,function (x) ifelse((abs(x[i+2]))>=tunned_value,1,0))
    gene_whole_all_tuned<-gene_whole_all[which(tunning_index==1),]
    tunning_index<-apply(gene_whole,1,function (x) ifelse((abs(x[i+2]))>=tunned_value,1,0))
    gene_whole_tuned<-gene_whole[which(tunning_index==1),]
    dim(gene_whole_tuned)
    unique(gene_whole_tuned[,2])
    table(gene_whole_tuned[,2])
    #geneate link data with 25th tissue
    links_gene<-gene_whole_tuned[,c(1,2,i+2)]
    colnames(links_gene)<-c("from","to","tissue_name")
    links_gene<-as.data.frame(links_gene)
    links_gene$edge_type<-ifelse(gene_whole_tuned[,i+2]>0,1,2)
    rownames(links_gene)<-rep(1:nrow(links_gene))
    
    self_index_all<-unique(c(gene_whole_all_tuned[,1],gene_whole_all_tuned[,2],gene_whole_tuned[,1],gene_whole_tuned[,2]))
    self_index<-unique(c(gene_whole_tuned[,1],gene_whole_tuned[,2]))
    #nodes_gene<-data.frame(cbind(self_index,self_size[self_index,(i+1)]))
    nodes_gene<-data.frame(self_size[self_size[,1] %in% self_index,c(1,(i+1))])
    colnames(nodes_gene)<-c("id","self_weight")
    #nodes_gene<-data.frame(self_size[,c((i+1))])
    #colnames(nodes_gene)<-c("self_weight")
    nodes_gene$node_type<-ifelse(nodes_gene$id %in% c(39,40,48,64,76),1, ifelse(nodes_gene$id %in% (self_index),2,3))
    
    #nodes_gene$node_type<-ifelse(nodes_gene$id %in% self_index,3,4)
    net <- graph.data.frame(links_gene, nodes_gene, directed=T) 
    class(net)
    net 
    
    # It's easy to access nodes, edges, and their attributes:
    E(net)
    V(net)$self_weight
    E(net)$tissue_name
    V(net) #total 47 clusters included in the network
    
    # Now we should be able to do this:
    #plot(net, edge.arrow.size=.2,vertex.label.family="Arial Black" )
    
    #adjust size of nodes and width of edges
    colrs_nodes <- c("red","orange","gray")
    V(net)$color <- colrs_nodes[V(net)$node_type]
    colrs_edge<-c("tomato","gray40")
    E(net)$color <- colrs_edge[E(net)$edge_type]
    V(net)$size <- (V(net)$self_weight)^0.7*4
    E(net)$width <- (abs(E(net)$tissue_name))^0.4*7
    
    #change coordinate to separate the nodes
    
    #set.seed(41414343)
    net.bg<-net
    set.seed(4151320)
    l<-layout_nicely(net.bg)
    l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
    l[20,]<-c(-0.9,0.1)
    l[16,]<-c(0.7,-0.5)
    l[46,]<-c(-0.7,-0.7)
    l[8,]<-c(0.7,-0.7)
    #l<-layout_in_circle(net.bg)
    #l<-layout_as_star(net.bg,center=V(net)[1:3])
    #l<-layout_as_bipartite(net.bg)
    #l<-layout_on_grid(net.bg)
    #l<-layout_on_sphere(net.bg)
    #l<-layout_with_sugiyama(net.bg)
    #l<-layout_with_dh(net.bg)
    # Normalize them so that they are in the -1, 1 interval:
    #l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
    gene_names_used<-as.character(gene_exact_names$gene.name[nodes_gene$id])
    #gene_names_used[which(gene_names_used=="CTSD\xca")]<-"CTSD"
    plot.igraph(net.bg, rescale=F, layout=l,edge.arrow.size=.1, edge.curved=.2,vertex.label=gene_names_used,
                vertex.label.color="blue",vertex.label.cex=1,mark.border=NA)
}

#generate a list of genes in the network
genelist_network<-gene_exact_names$gene.name[self_index_all]

##network for x=0 at i=10 after some tuning
par(mfrow = c(1, 1))
par(cex=1)
par(mar = c(0, 0, 0, 0), oma = c(6, 6, 3, 3))
tunned_value<-0.03
for (i in c(10))
{
    tunning_index<-apply(gene_whole_all,1,function (x) ifelse((abs(x[i+2]))>=tunned_value,1,0))
    gene_whole_all_tuned<-gene_whole_all[which(tunning_index==1),]
    tunning_index<-apply(gene_whole,1,function (x) ifelse((abs(x[i+2]))>=tunned_value,1,0))
    gene_whole_tuned<-gene_whole[which(tunning_index==1),]
    dim(gene_whole_tuned)
    unique(gene_whole_tuned[,2])
    table(gene_whole_tuned[,2])
    #geneate link data with 25th tissue
    links_gene<-gene_whole_tuned[,c(1,2,i+2)]
    colnames(links_gene)<-c("from","to","tissue_name")
    links_gene<-as.data.frame(links_gene)
    links_gene$edge_type<-ifelse(gene_whole_tuned[,i+2]>0,1,2)
    rownames(links_gene)<-rep(1:nrow(links_gene))
    
    self_index_all<-unique(c(gene_whole_all_tuned[,1],gene_whole_all_tuned[,2],gene_whole_tuned[,1],gene_whole_tuned[,2]))
    self_index<-unique(c(gene_whole_tuned[,1],gene_whole_tuned[,2]))
    #nodes_gene<-data.frame(cbind(self_index,self_size[self_index,(i+1)]))
    nodes_gene<-data.frame(self_size[self_size[,1] %in% self_index,c(1,(i+1))])
    colnames(nodes_gene)<-c("id","self_weight")
    #nodes_gene<-data.frame(self_size[,c((i+1))])
    #colnames(nodes_gene)<-c("self_weight")
    nodes_gene$node_type<-ifelse(nodes_gene$id %in% c(39,40,76),1, ifelse(nodes_gene$id %in% (self_index),2,3))
    
    #nodes_gene$node_type<-ifelse(nodes_gene$id %in% self_index,3,4)
    net <- graph.data.frame(links_gene, nodes_gene, directed=T) 
    class(net)
    net 
    
    # It's easy to access nodes, edges, and their attributes:
    E(net)
    V(net)$self_weight
    E(net)$tissue_name
    V(net) #total 47 clusters included in the network
    
    # Now we should be able to do this:
    #plot(net, edge.arrow.size=.2,vertex.label.family="Arial Black" )
    
    #adjust size of nodes and width of edges
    colrs_nodes <- c("red","orange","gray")
    V(net)$color <- colrs_nodes[V(net)$node_type]
    colrs_edge<-c("tomato","gray40")
    E(net)$color <- colrs_edge[E(net)$edge_type]
    V(net)$size <- (V(net)$self_weight)^0.7*4
    E(net)$width <- (abs(E(net)$tissue_name))^0.4*7
    
    #change coordinate to separate the nodes
    
    #set.seed(41414343)
    net.bg<-net
    set.seed(4151320)
    l<-layout_nicely(net.bg)
    l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
    #l[20,]<-c(-0.9,0.1)
    #l[16,]<-c(0.7,-0.5)
    #l[46,]<-c(-0.7,-0.7)
    #l[8,]<-c(0.7,-0.7)
    #l<-layout_in_circle(net.bg)
    #l<-layout_as_star(net.bg,center=V(net)[1:3])
    #l<-layout_as_bipartite(net.bg)
    #l<-layout_on_grid(net.bg)
    #l<-layout_on_sphere(net.bg)
    #l<-layout_with_sugiyama(net.bg)
    #l<-layout_with_dh(net.bg)
    # Normalize them so that they are in the -1, 1 interval:
    #l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
    gene_names_used<-as.character(gene_exact_names$gene.name[nodes_gene$id])
    #gene_names_used[which(gene_names_used=="CTSD\xca")]<-"CTSD"
    plot.igraph(net.bg, rescale=F, layout=l*0.99,edge.arrow.size=.15, edge.curved=.2,vertex.label=gene_names_used,
                vertex.label.color="blue",vertex.label.cex=1,mark.border=NA)
}

##network for x=0 at i=100 after some tuning
par(mfrow = c(1, 1))
par(cex=1)
par(mar = c(0, 0, 0, 0), oma = c(6, 6, 3, 3))
tunned_value<-0.03
for (i in c(100))
{
    tunning_index<-apply(gene_whole_all,1,function (x) ifelse((abs(x[i+2]))>=tunned_value,1,0))
    gene_whole_all_tuned<-gene_whole_all[which(tunning_index==1),]
    tunning_index<-apply(gene_whole,1,function (x) ifelse((abs(x[i+2]))>=tunned_value,1,0))
    gene_whole_tuned<-gene_whole[which(tunning_index==1),]
    dim(gene_whole_tuned)
    unique(gene_whole_tuned[,2])
    table(gene_whole_tuned[,2])
    #geneate link data with 25th tissue
    links_gene<-gene_whole_tuned[,c(1,2,i+2)]
    colnames(links_gene)<-c("from","to","tissue_name")
    links_gene<-as.data.frame(links_gene)
    links_gene$edge_type<-ifelse(gene_whole_tuned[,i+2]>0,1,2)
    rownames(links_gene)<-rep(1:nrow(links_gene))
    
    self_index_all<-unique(c(gene_whole_all_tuned[,1],gene_whole_all_tuned[,2],gene_whole_tuned[,1],gene_whole_tuned[,2]))
    self_index<-unique(c(gene_whole_tuned[,1],gene_whole_tuned[,2]))
    #nodes_gene<-data.frame(cbind(self_index,self_size[self_index,(i+1)]))
    nodes_gene<-data.frame(self_size[self_size[,1] %in% self_index,c(1,(i+1))])
    colnames(nodes_gene)<-c("id","self_weight")
    #nodes_gene<-data.frame(self_size[,c((i+1))])
    #colnames(nodes_gene)<-c("self_weight")
    nodes_gene$node_type<-ifelse(nodes_gene$id %in% c(39,40,48,76),1, ifelse(nodes_gene$id %in% (self_index),2,3))
    
    #nodes_gene$node_type<-ifelse(nodes_gene$id %in% self_index,3,4)
    net <- graph.data.frame(links_gene, nodes_gene, directed=T) 
    class(net)
    net 
    
    # It's easy to access nodes, edges, and their attributes:
    E(net)
    V(net)$self_weight
    E(net)$tissue_name
    V(net) #total 47 clusters included in the network
    
    # Now we should be able to do this:
    #plot(net, edge.arrow.size=.2,vertex.label.family="Arial Black" )
    
    #adjust size of nodes and width of edges
    colrs_nodes <- c("red","orange","gray")
    V(net)$color <- colrs_nodes[V(net)$node_type]
    colrs_edge<-c("tomato","gray40")
    E(net)$color <- colrs_edge[E(net)$edge_type]
    V(net)$size <- (V(net)$self_weight)^0.7*4
    E(net)$width <- (abs(E(net)$tissue_name))^0.4*7
    
    #change coordinate to separate the nodes
    
    #set.seed(41414343)
    net.bg<-net
    set.seed(4151320)
    l<-layout_nicely(net.bg)
    l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
    l[9,]<-c(0,0.1)
    l[8,]<-c(-0.9,-0.3)
    #l[46,]<-c(-0.7,-0.7)
    #l[8,]<-c(0.7,-0.7)
    #l<-layout_in_circle(net.bg)
    #l<-layout_as_star(net.bg,center=V(net)[1:3])
    #l<-layout_as_bipartite(net.bg)
    #l<-layout_on_grid(net.bg)
    #l<-layout_on_sphere(net.bg)
    #l<-layout_with_sugiyama(net.bg)
    #l<-layout_with_dh(net.bg)
    # Normalize them so that they are in the -1, 1 interval:
    #l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
    gene_names_used<-as.character(gene_exact_names$gene.name[nodes_gene$id])
    #gene_names_used[which(gene_names_used=="CTSD\xca")]<-"CTSD"
    plot.igraph(net.bg, rescale=F, layout=l*0.99,edge.arrow.size=.12, edge.curved=.2,vertex.label=gene_names_used,
                vertex.label.color="blue",vertex.label.cex=1,mark.border=NA)
}

par(mfrow = c(1, 3))
tunned_value<-0.01
for (i in c(10,100,250))
{
    tunning_index<-apply(gene_whole_cov,1,function (x) ifelse((abs(x[i+2]))>=tunned_value,1,0))
    gene_whole_cov_tuned<-gene_whole_cov[which(tunning_index==1),]
    dim(gene_whole_cov_tuned)
    unique(gene_whole_cov_tuned[,2])
    table(gene_whole_cov_tuned[,2])
    #geneate link data with 25th tissue
    links_gene<-gene_whole_cov_tuned[,c(1,2,i+2)]
    colnames(links_gene)<-c("from","to","tissue_name")
    links_gene<-as.data.frame(links_gene)
    links_gene$edge_type<-ifelse(gene_whole_cov_tuned[,i+2]>0,1,2)
    rownames(links_gene)<-rep(1:nrow(links_gene))
    
    #self_index_all<-unique(c(gene_whole_tuned[,1],gene_whole_tuned[,2],gene_whole_cov_tuned[,1],gene_whole_cov_tuned[,2]))
    self_index_cov<-unique(c(gene_whole_cov_tuned[,1],gene_whole_cov_tuned[,2]))
    self_index_all<-self_index_cov
    #nodes_gene<-data.frame(cbind(self_index,self_size[self_index,(i+1)]))
    nodes_gene<-abs(data.frame(self_size_cov[self_size_cov[,1] %in% self_index_all,c(1,(i+1))]))
    colnames(nodes_gene)<-c("id","self_weight")
    #nodes_gene<-data.frame(self_size[,c((i+1))])
    #colnames(nodes_gene)<-c("self_weight")
    nodes_gene$node_type<-ifelse(nodes_gene$id %in% (self_index_cov),3,4)
    net <- graph.data.frame(links_gene, nodes_gene, directed=T) 
    class(net)
    net 
    
    # It's easy to access nodes, edges, and their attributes:
    E(net)
    V(net)$self_weight
    E(net)$tissue_name
    V(net) #total 34 clusters included in the network
    
    # Now we should be able to do this:
    #plot(net, edge.arrow.size=.2,vertex.label.family="Arial Black" )
    
    #adjust size of nodes and width of edges
    colrs_nodes <- c("red", "green","orange","gray")
    V(net)$color <- colrs_nodes[V(net)$node_type]
    colrs_edge<-c("tomato","gray40")
    E(net)$color <- colrs_edge[E(net)$edge_type]
    V(net)$size <- (V(net)$self_weight)^0.7*4
    E(net)$width <- (abs(E(net)$tissue_name))^0.4*7
    net.bg<-net
    set.seed(4151320)
    l<-layout_nicely(net.bg)
    l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
    #l[1,2]<-0.3
    #l[3,2]<-0.7
    #l[25,2]<--0.5
    #l[23,2]<-0
    #l[14,2]<--0.6
    #l[8,]<-c(-0.9,-0.35)
    #l[16,1]<-0.4
    #l[15,1]<--0.2
    #l<-layout_in_circle(net.bg)
    #l<-layout_as_star(net.bg,center=V(net)[1:3])
    #l<-layout_as_bipartite(net.bg)
    #l<-layout_on_grid(net.bg)
    #l<-layout_on_sphere(net.bg)
    #l<-layout_with_sugiyama(net.bg)
    #l<-layout_with_dh(net.bg)
    # Normalize them so that they are in the -1, 1 interval:
    l <- layout.norm(l, ymin=-1, ymax=1, xmin=-1, xmax=1)
    gene_names_used<-as.character(gene_exact_names$gene.name[nodes_gene$id])
    gene_names_used[which(gene_names_used=="CTSD\xca")]<-"CTSD"
    plot.igraph(net.bg, rescale=F, layout=l*0.9,edge.arrow.size=.12, edge.curved=.2,vertex.label=gene_names_used,
                vertex.label.color="blue",vertex.label.cex=1.2,mark.border=NA)
}
