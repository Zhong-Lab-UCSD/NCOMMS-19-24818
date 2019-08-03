### Script with the function used to plot super enhancer network

mes_markers = c("SMAD3","VIM","CTGF","SERPINE1","LINC00607","FN1","VWF")

# super enhancers for mesenchymal markers
SE_mes_markers = c()
for (i in mes_markers){
  SE_mes_markers=rbind(SE_mes_markers,super_enhancers_sort_annotated_new[which(super_enhancers_sort_annotated_new[,"SE_genes"] %like% paste0("%",i,"%")),])
}
SE_mes_markers_index = as.numeric(rownames(SE_mes_markers))

##############################################

plot_network_NEW <- function(indexes, output_path, vertex_size=1, v_label_cex=1, vertex_color = "#0000ff",
                             vertex_to_increase=list(), 
                             increase_v_size=list(), 
                             increase_v_color = c(),
                             labels_to_increase=c(), 
                             increase_label_size=1.5,
                             edge_width = 0.5, arrow_size = 0.01, arrow_width = 0.01,
                             edges_to_increase=list(), # matrix of two columns where each row represent two vertexes connected by the edge that should be increased
                             increase_e_size=c(), increase_arrow_size=c(), increase_arrow_width=c(), increase_e_color=c(),
                             mes_markers_names = FALSE,
                             my_layout = layout.sphere,
                             plot_chr_label = FALSE,
                             add_single_nodes = c()) 

{
  
  SE_mes_markers = c()
  SE_mes_markers_gene = c()
  for (i in mes_markers){
    SE_mes_markers=rbind(SE_mes_markers,super_enhancers_sort_annotated_new[which(super_enhancers_sort_annotated_new[,"SE_genes"] %like% paste0("%",i,"%")),])
    SE_mes_markers_gene = c(SE_mes_markers_gene,rep(i,nrow(super_enhancers_sort_annotated_new[which(super_enhancers_sort_annotated_new[,"SE_genes"] %like% paste0("%",i,"%")),])))
  }
  SE_mes_markers_index = as.numeric(rownames(SE_mes_markers))
  SE_mes_markers_index_full = cbind(SE_mes_markers_index,SE_mes_markers_gene,as.character(SE_mes_markers$labels_new))
  
  SE_network_info = matrix(seq(1:nrow(super_enhancers_sort_annotated_new)),nrow(super_enhancers_sort_annotated_new),1)
  colnames(SE_network_info) = "SE_index"
  SE_network_info = data.frame(SE_network_info)
  SE_network_info$type <- ifelse(SE_network_info$SE_index %in% SE_mes_markers[,"SE_index_new"],'MES marker',"Non marker")
  SE_network_info$SE_index = super_enhancers_sort_annotated_new$labels_new
  SE_network_info$SE_seq = seq(1:nrow(super_enhancers_sort_annotated_new))
  SE_network_info$SE_color = vertex_color
  SE_network_info$label_dist = 0
  SE_network_info$label_degree = -pi/4
  # Adding this column to plot chromosome label on intrachromosomal matrix (plot_chr_label=TRUE)
  temp = gsub("-.*","",SE_network_info$SE_index)
  temp = gsub("C","chr",temp)
  SE_network_info$SE_chr = temp
  
  if (mes_markers_names == TRUE){
    SE_network_info$SE_index = mapvalues(SE_network_info$SE_index,SE_mes_markers_index_full[,3],SE_mes_markers_index_full[,2],warn_missing = FALSE) # replace mes markers with gene name
  }
  
  indexes_temp_row = as.character(super_enhancers_sort_annotated_new[indexes[,"row"],"labels_new"])
  indexes_temp_col = as.character(super_enhancers_sort_annotated_new[indexes[,"col"],"labels_new"])
  
  indexes_temp_label = cbind(indexes_temp_row,indexes_temp_col)
  d = data.frame(indexes_temp_label)
  d$e_width = edge_width
  d$arr_size = arrow_size
  d$arr_width = arrow_width
  d$e_color = "darkgrey"
  
  if (length(edges_to_increase) > 0){
    for (j in 1:length(edges_to_increase)){
      edges_to_increase_j = edges_to_increase[[j]]
      edges_to_increase_temp = cbind(as.character(super_enhancers_sort_annotated_new[edges_to_increase_j[,1],"labels_new"]),as.character(super_enhancers_sort_annotated_new[edges_to_increase_j[,2],"labels_new"]))
      for (i in 1:nrow(edges_to_increase_temp)){
        d[which(d[,"indexes_temp_row"]==edges_to_increase_temp[i,1] & d[,"indexes_temp_col"]==edges_to_increase_temp[i,2]),"e_width"] = increase_e_size[j]
        d[which(d[,"indexes_temp_row"]==edges_to_increase_temp[i,2] & d[,"indexes_temp_col"]==edges_to_increase_temp[i,1]),"e_width"] = increase_e_size[j]
        
        d[which(d[,"indexes_temp_row"]==edges_to_increase_temp[i,1] & d[,"indexes_temp_col"]==edges_to_increase_temp[i,2]),"e_color"] = increase_e_color[j]
        d[which(d[,"indexes_temp_row"]==edges_to_increase_temp[i,2] & d[,"indexes_temp_col"]==edges_to_increase_temp[i,1]),"e_color"] = increase_e_color[j]
        
        if (length(increase_arrow_size) > 0){
          d[which(d[,"indexes_temp_row"]==edges_to_increase_temp[i,1] & d[,"indexes_temp_col"]==edges_to_increase_temp[i,2]),"arr_size"] = increase_arrow_size[j]
          d[which(d[,"indexes_temp_row"]==edges_to_increase_temp[i,2] & d[,"indexes_temp_col"]==edges_to_increase_temp[i,1]),"arr_size"] = increase_arrow_size[j]
        }
        
        if (length(increase_arrow_width) > 0){
          d[which(d[,"indexes_temp_row"]==edges_to_increase_temp[i,1] & d[,"indexes_temp_col"]==edges_to_increase_temp[i,2]),"arr_width"] = increase_arrow_width[j]
          d[which(d[,"indexes_temp_row"]==edges_to_increase_temp[i,2] & d[,"indexes_temp_col"]==edges_to_increase_temp[i,1]),"arr_width"] = increase_arrow_width[j]
        }
      }
    }
  }
  
  
  if (mes_markers_names == TRUE){
    d$indexes_temp_row = mapvalues(d$indexes_temp_row,SE_mes_markers_index_full[,3],SE_mes_markers_index_full[,2],warn_missing = FALSE)
    d$indexes_temp_col = mapvalues(d$indexes_temp_col,SE_mes_markers_index_full[,3],SE_mes_markers_index_full[,2],warn_missing = FALSE)
  }
  
  graph_object = graph_from_data_frame(d, directed = TRUE, vertices = NULL)
  
  if (length(add_single_nodes)>0){
    single_nodes_label = as.character(super_enhancers_sort_annotated_new[add_single_nodes,"labels_new"])
    graph_object = add_vertices(graph_object, length(add_single_nodes), name = single_nodes_label)
  } 
  
  if (length(vertex_to_increase)>0){
    SE_network_info$v_size = vertex_size
    # iterate on every group of vertexes
    for (i in 1:length(vertex_to_increase)){
      vertex_to_increase_label = as.character(super_enhancers_sort_annotated_new[vertex_to_increase[[i]],"labels_new"])
      incr_v_size = increase_v_size[[i]]
      
      if (length(incr_v_size)==2){
        hubs = as.numeric(degree(graph_object)[vertex_to_increase_label])
        increase_v_size_norm = ((hubs - min(hubs))/(max(hubs)-min(hubs)))*(incr_v_size[2]-incr_v_size[1]) + incr_v_size[1]
        SE_network_info[which(SE_network_info[,"SE_index"] %in% vertex_to_increase_label),"v_size"] = increase_v_size_norm
        ### Set label distance parameter depending on the dimension of the vertex
        if (incr_v_size[1] < 4){
          temp = SE_network_info[which(SE_network_info[,"v_size"] >= incr_v_size[1] & SE_network_info[,"v_size"] < 4),]
          # Normalize label distance for vertexes between 1 and 4
          min_label_dist = 0.2
          max_label_dist = 0.6
          temp_label_dist = ((temp$v_size - 1)/(4-1))*(max_label_dist-min_label_dist) + min_label_dist
          SE_network_info[rownames(temp),"label_dist"] = temp_label_dist
        }
      }
      else {
        SE_network_info[which(SE_network_info[,"SE_index"] %in% vertex_to_increase_label),"v_size"] = incr_v_size
        ### Set label distance parameter depending on the dimension of the vertex
        if (incr_v_size < 4){
          temp = SE_network_info[which(SE_network_info[,"v_size"] == incr_v_size),]
          # Normalize label distance for vertexes between 1 and 4
          min_label_dist = 0.5
          max_label_dist = 2
          temp_label_dist = ((temp$v_size - 1)/(4-1))*(max_label_dist-min_label_dist) + min_label_dist
          SE_network_info[rownames(temp),"label_dist"] = temp_label_dist
        }
      }
      
      # Change color of hubs
      if (length(increase_v_color) > 0){
        SE_network_info[which(SE_network_info[,"SE_index"] %in% vertex_to_increase_label),"SE_color"] = increase_v_color[i]
      }
      
      V(graph_object)$size=as.numeric(SE_network_info$v_size[match(V(graph_object)$name,SE_network_info$SE_index)])
      V(graph_object)$color=as.character(SE_network_info$SE_color[match(V(graph_object)$name,SE_network_info$SE_index)])
      V(graph_object)$label.dist=as.numeric(SE_network_info$label_dist[match(V(graph_object)$name,SE_network_info$SE_index)])
    }
    

  }
  else {
    V(graph_object)$size = vertex_size
    V(graph_object)$color = vertex_color
  }
  
  
  if (length(labels_to_increase)>0){
    SE_network_info$label_size = 1
    SE_network_info[which(SE_network_info[,"SE_index"] %in% labels_to_increase),"label_size"] = increase_label_size
    V(graph_object)$label.cex=as.numeric(SE_network_info$label_size[match(V(graph_object)$name,SE_network_info$SE_index)])
  }
  else {
    V(graph_object)$label.cex = v_label_cex
  }
  
  if (length(vertex_to_increase)>0){
    total_vertex_to_increase = c()
    for (i in 1:length(vertex_to_increase)){
      total_vertex_to_increase = c(total_vertex_to_increase,vertex_to_increase[[i]])
    }
  }
  
  V(graph_object)$label.degree=as.numeric(SE_network_info$label_degree[match(V(graph_object)$name,SE_network_info$SE_index)])
  
  # To have SE index as label
  V(graph_object)$name = mapvalues(V(graph_object)$name,SE_network_info[,1],SE_network_info[,3],warn_missing = FALSE)
  
  if (plot_chr_label == TRUE){
    total_vertex_to_increase = hg38_chromosomes
    # Select one node per chromosome that will carry the chromosome label: we select the one with the highest DON
    for (i in hg38_chromosomes){
      se_ind = SE_network_info[which(SE_network_info$SE_chr == i),3]
      temp_ind = V(graph_object)$name[which(V(graph_object)$name %in% se_ind)]
      if (length(temp_ind) > 0){
        DONs = degree(graph_object, v = temp_ind)
        x = names(DONs[which(DONs==max(DONs))])[1] # index that will be replaced with i ([1] in case there are multiple nodes with the same DON)
        V(graph_object)$name[which(V(graph_object)$name == x)] = i
      }
    }
  }
  png(output_path,width = 1500, height = 1500, units = "px")
  if (plot_chr_label == TRUE){
    p<-plot(graph_object, vertex.label = ifelse(V(graph_object)$name %in% total_vertex_to_increase, V(graph_object)$name, NA), layout=my_layout,
            vertex.label.color = "black", vertex.label.cex = 5, vertex.label.family = "Helvetica",
            edge.arrow.size = E(graph_object)$arr_size, edge.arrow.width = E(graph_object)$arr_width, edge.width = E(graph_object)$e_width, edge.color = E(graph_object)$e_color)
  } else {
    p<-plot(graph_object, vertex.label = ifelse(V(graph_object)$name %in% total_vertex_to_increase, V(graph_object)$name, NA), layout=my_layout,
            vertex.label.color = "black", vertex.label.family = "Helvetica", #vertex.label.degree = pi/4,
            edge.arrow.size = E(graph_object)$arr_size, edge.arrow.width = E(graph_object)$arr_width, edge.width = E(graph_object)$e_width, edge.color = E(graph_object)$e_color)
  } 
  if (color_markers == TRUE){
    legend("bottomleft", legend=c("MES marker", "Non marker") , col = c("red","blue") , bty = "n", pch=20 , pt.cex = 2, cex = 1.2, text.col=c("red","blue") , horiz = FALSE, inset = c(0.01, 0.01))
  }
  print(p)
  dev.off()
}
