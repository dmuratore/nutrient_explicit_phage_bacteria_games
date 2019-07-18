## The purpose of this script is to generate a figure showing
## the potential paths of heteroclinic cycles given
## certain parameter regimes in the environmental feedback
## bimatrix game model of host siderophore/viral Fe tail ecological
## dynamics

## Loading required libraries
## If you need to install igraph first
## install.packages('igraph')
library(igraph)


## Initializing payoff matrices
## These qs respond to the qs in get_pars.m, but in terms of original
## dynamical model as described in main text section 2.2:
## q1: a-c
## q2: b-d
## q3: alpha-beta
## q4: gamma-delta
## q1p: a'-c'
## q2p: b'-d'
## q3p: alpha'-beta'
## q4p: gamma'-delta'
## thetax,thetay: environmental restoration parameters
q1<-1
q2<-2
q3<--1
q4<--3
q1p<--1
q2p<--1.1
q3p<-1
q4p<-2
thetax<-4
thetay<-1


## Make characteristic matrix (method from Hofbauer 1994)
c_mat<-matrix(data=c(q2,q2p,q1,q1p,0,0,0,0,
                     0,0,0,0,-q2,-q2p,-q1,-q1p,
                     q4,q4p,0,0,q3,q3p,0,0,
                     0,0,-q4,-q4p,0,0,-q3,-q3p,
                     -1,0,-(1+thetay),0,thetax-1,0,thetax-(1+thetay),0,
                     0,1,0,1+thetay,0,-(thetax-1),0,-(thetax-(1+thetay))),
              ncol=6)

## Determine which elements are positive 
## (indicate direction of travel for dynamics)
c_plus<-apply(c_mat,1,function(x) which(x>0))

## Now we will analyze this matrix as a directed graph to determine 
## heteroclinic cycles, any cycles in the graph are heteroclinic
## cycles in the network

## Create the first part of the edgelist for the network
first_edges<-c()
for(i in 1:length(c_plus)){
  first_edges<-c(first_edges,rep(i,length(c_plus[[i]])))
}

## Now we have to figure out which nodes each node identified in the
## previous command connects to

## We will do this with a complicated if loop sorry
links<-do.call('c',c_plus) #Find all outward directed edges
second_edges<-c()
for(i in 1:length(links)){
## The edge assignment is based off of how the characteristic matrix is
## organized, please see appendix B.0.4 for details
  if(first_edges[i]==1){
    if(links[i]==1){
      second_edges[i]<-5
    }
    if(links[i]==3){
      second_edges[i]<-3
    }
    if(links[i]==5){
      second_edges[i]<-2
    }
  }
  if(first_edges[i]==2){
    if(links[i]==1){
      second_edges[i]<-6
    }
    if(links[i]==3){
      second_edges[i]<-4
    }
    if(links[i]==6){
      second_edges[i]<-1
    }
  }
  if(first_edges[i]==3){
    if(links[i]==1){
      second_edges[i]<-7
    }
    if(links[i]==4){
      second_edges[i]<-1
    }
    if(links[i]==5){
      second_edges[i]<-4
    }
  }
  if(first_edges[i]==4){
    if(links[i]==1){
      second_edges[i]<-8
    }
    if(links[i]==4){
      second_edges[i]<-2
    }
    if(links[i]==6){
      second_edges[i]<-3
    }
  }
  if(first_edges[i]==5){
    if(links[i]==2){
      second_edges[i]<-1
    }
    if(links[i]==3){
      second_edges[i]<-7
    }
    if(links[i]==5){
      second_edges[i]<-6
    }
  }
  if(first_edges[i]==6){
    if(links[i]==2){
      second_edges[i]<-2
    }
    if(links[i]==3){
      second_edges[i]<-8
    }
    if(links[i]==6){
      second_edges[i]<-5
    }
  }
  if(first_edges[i]==7){
    if(links[i]==2){
      second_edges[i]<-3
    }
    if(links[i]==4){
      second_edges[i]<-5
    }
    if(links[i]==5){
      second_edges[i]<-8
    }
  }
  if(first_edges[i]==8){
    if(links[i]==2){
      second_edges[i]<-4
    }
    if(links[i]==4){
      second_edges[i]<-6
    }
    if(links[i]==6){
      second_edges[i]<-7
    }
  }
}

## Combining into edgelist for graph visualization
edgelist<-c(rbind(first_edges,second_edges))

## Make directed graph
hc_graph<-make_directed_graph(edgelist)
## Naming nodes
vertex.attributes(hc_graph)<-list(name=c('x=0,y=0,n=0',
                                  'x=0,y=0,n=1',
                                  'x=0,y=1,n=0',
                                  'x=0,y=1,n=1',
                                  'x=1,y=0,n=0',
                                  'x=1,y=0,n=1',
                                  'x=1,y=1,n=0',
                                  'x=1,y=1,n=1'))

## Identifying cycles
cycle_function<-function(g,vert){
  paths<-all_simple_paths(g,from=vert)
  final_nodes<-sapply(paths,function(x) tail(x,1))
  cycle_ids<-sapply(final_nodes,function(x) are_adjacent(g,x,vert))
  output_paths<-paths[which(cycle_ids==TRUE)]
  output_paths<-lapply(output_paths,function(x) c(x,x[1]))
  return(output_paths)
}
all_cycles<-do.call('c',sapply(1:vcount(hc_graph),cycle_function,g=hc_graph))
numeric_paths<-lapply(all_cycles,function(x) sort(unique(as.numeric(x))))
unique_paths<-all_cycles[which(duplicated(numeric_paths)==FALSE)]


## Changing cycle aesthetic attributes for plotting
emph_cycle<-function(g,p){
  path_graph<-g
  edges_to_color<-c(head(p,-1),tail(p,-1))
  color_ids<-get.edge.ids(path_graph,edges_to_color)
  E(path_graph)$color<-NA
  E(path_graph)$color[color_ids]<-'black'
  return(path_graph)
}
## Implementing aesthetic changes
graphs_for_plotting<-lapply(unique_paths,emph_cycle,g=hc_graph)
## Initializing figure
jpeg('../figures/figure_4.jpg',units='in',width=8.5,height=11,res=600)
par(mfrow=c(3,2))
## Setting vertex colors
vertex_fills<-list(c(0,0,0,1,1,1),
                   c(0,0,1,1,1,0),
                   c(0,1,0,1,0,1),
                   c(0,1,1,1,0,0),
                   c(1,0,0,0,1,1),
                   c(1,0,1,0,1,0),
                   c(1,1,0,0,0,1),
                   c(1,1,1,0,0,0))
vertex_colset<-list(c('green','magenta','blue','white','white','white'))
## Making subplots for each heteroclinic cycle in the network
set.seed(888)
coords<-layout_with_fr(hc_graph)
for(i in 1:length(graphs_for_plotting)){
plot(graphs_for_plotting[[i]],
     layout=coords,
     vertex.size=35,
     vertex.shape='pie',
     vertex.pie=vertex_fills,
     vertex.pie.color=vertex_colset,
     vertex.label.font=2,
     vertex.label.dist=0,
     edge.arrow.size=0.75,
     margin=rep(0,4),
     vertex.label=NA)
text(-0.8,-0.8,'n',cex=2.2,srt=-40)
text(0.6,-1.2,'y',cex=2.2,srt=15)
text(-1.2,0.2,'x',cex=2.2,srt=90)
#legend(x=-1.55,y=1.75,fill=c('green','magenta','blue'),legend=c('All Hosts Produce Siderophores (x=1)',
#                                                                'All Viruses Have Iron Tails (y=1)',
#                                                                'Iron Replete (n=1)'),
#       bty='n',border='black',x.intersp=0.2,cex=1.4)
}
dev.off()

## Preparing a figure legend to use the colors in bottom legend for F4. 
jpeg('../figures/figure_4_legend.jpg',units='in',width=8.5,height=4,res=600)
plot(hc_graph,
     layout=coords,
     vertex.size=0,
     vertex.shape='pie',
     vertex.pie=vertex_fills,
     vertex.pie.color=vertex_colset,
     vertex.label.font=2,
     vertex.label.dist=0,
     edge.arrow.size=0,
     margin=rep(0,4),
     vertex.label=NA,
    edge.color='white')
#text(-0.8,-0.8,'n',cex=2.2,srt=-40)
#text(0.6,-1.2,'y',cex=2.2,srt=15)
#text(-1.2,0.2,'x',cex=2.2,srt=90)
legend(x=-2,y=1,fill=c('green','magenta','blue'),legend=c('All Hosts Produce Siderophores (x=1)',
                                                                'All Viruses Have Iron Tails (y=1)',
                                                                'Iron Replete (n=1)'),
       bty='n',border='black',x.intersp=0.2,cex=1.4)
dev.off()
