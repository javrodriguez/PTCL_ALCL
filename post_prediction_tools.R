# FUNCTIONS
point_score = function(locus,radius=500000,M, pseudocount=0.1){
  l_edge = max((locus - radius),1)
  r_edge = min((locus + radius),nrow(M))
  l_mask = M[,l_edge:locus]
  r_mask = M[,locus:r_edge]
  center_mask = M[l_edge:locus,locus:r_edge]
  score = (max(mean(l_mask,na.rm=T),mean(r_mask,na.rm=T)) +  pseudocount) / (mean(center_mask) +  pseudocount)
  return(score)
}

region_score = function(M,resolution=8192,radius=500000 ,pseudocount=0.1){
  M = as.matrix(M)
  pixel_radius = ceiling(radius/resolution)
  scores = unlist(lapply(1:nrow(M),point_score,radius=pixel_radius,M=M,pseudocount=pseudocount))
  return(scores)
}

save_scores = function(prediction_file,inpdir,window=2097152,resolution=8192,radius=500000){
  # Data reading
  coord = gsub(x=prediction_file,pattern = ".npy",replacement = "")
  chr_start=unlist(strsplit(coord,"_"))
  outfile= paste0("contact_maps_",coord,".pdf")
  m_pred <- as.data.frame(np$load(file.path(inpdir,prediction_file)))
  rscl = abs(min(m_pred))
  m_pred = m_pred+rscl
  
  # HiC values
  preds = m_pred[upper.tri(m_pred, diag = F)]
  
  # Insulation values
  prediction_scores = region_score(M=m_pred,resolution=resolution,radius=radius,pseudocount=pseudocount)
  
  # Bin coordinates
  n_bins=nrow(m_pred)
  resolution=window/(n_bins)
  chr=chr_start[1]
  start=as.numeric(chr_start[2])
  end=start+window-1
  start_seq = seq(start,end,resolution)
  end_seq = seq(start_seq[2]-1,end,resolution)
  
  return(list(chr=rep(chr,n_bins),start=start_seq,end=end_seq,
              pred_scores=prediction_scores,window=rep(coord,n_bins)))
}

make_scores_bed = function(results,inpdir,outname="prediction_scores.csv"){
  chr=vector(mode = "character")
  start=vector(mode = "numeric")
  end=vector(mode = "numeric")
  scores=vector(mode = "numeric")
  window=vector(mode = "character")

  for(i in 1:length(results)){
    xi=results[[i]]
    for(j in 1:length(xi)){
      xij=unlist(xi[j])
      if(j==1){ chr = c(chr,xij) }
      if(j==2){ start = c(start,xij) }
      if(j==3){ end = c(end,xij) }
      if(j==4){ scores = c(scores,xij) }
      if(j==5){ window = c(window,xij) }
    }
  }
  scores_bed = data.frame(chr=chr,start=start,end=end,score=scores,window=window,stringsAsFactors = F)
  write.csv(scores_bed,file.path(inpdir,outname),row.names = F)
}

compute_scores_parallel = function(sample_name,save_target=F,n_cores=6){
  print(sample_name)
  inpdir = file.path(wd,sample_name,"prediction/npy/")
  matrix_files = list.files(inpdir,pattern = ".npy")
  target_files = matrix_files[grep("_target.npy",matrix_files,invert = F)]
  prediction_files = matrix_files[grep("_target.npy",matrix_files,invert = T)]
  if(length(prediction_files) != n_windows){ print(paste0("WARNING: ",sample_name)," does not have ",n_windows," matrices.") }

  print("Prediction scores...")
  results = mclapply(prediction_files,save_scores,inpdir=inpdir,mc.cores = n_cores)
  make_scores_bed(results,inpdir=inpdir)

  if(save_target){
      print("Target scores...")
      results = mclapply(target_files,save_scores,inpdir=inpdir,mc.cores = n_cores)
      make_scores_bed(results,inpdir=inpdir,outname="target_scores.csv")
  }
}

save_hic = function(matrix_file,inpdir,window=2097152,resolution=8192,radius=500000){
  # Data reading
  coord = gsub(x=matrix_file,pattern = "_target.npy",replacement = "")
  coord = gsub(x=coord,pattern = ".npy",replacement = "")
  chr_start=unlist(strsplit(coord,"_"))
  outfile= paste0("contact_maps_",coord,".pdf")
  m <- as.data.frame(np$load(file.path(inpdir,matrix_file)))
  rscl = abs(min(m))
  m = m+rscl
  
  # Bin coordinates
  n_bins=nrow(m)
  resolution=window/(n_bins)
  chr=chr_start[1]
  start=as.numeric(chr_start[2])
  end=start+window-1
  start_seq = seq(start,end,resolution)

  # Create matrix with coordinate index
  m_idx=as.data.frame(matrix(nrow=n_bins,ncol=n_bins))
  for (i in 1:n_bins){
      for (j in 1:n_bins){
        m_idx[i,j]=paste0(start_seq[i],":",start_seq[j])
      }
    }
  # Get HiC values and coords
  hic = m[upper.tri(m, diag = F)]
  coords = m_idx[upper.tri(m_idx, diag = F)]
  return(list(pred=hic,coord=coords))
}

make_hic_bed = function(results,inpdir,chr,outname="prediction_hic"){
  hic=vector(mode = "numeric")
  coords=vector(mode = "character")

  for(i in 1:length(results)){
    xi=results[[i]]
    for(j in 1:length(xi)){
      xij=unlist(xi[j])
      if(j==1){ hic = c(hic,xij) }
      if(j==2){ coords = c(coords,xij) }
    }
  }
  hic_bed = data.frame(coords=coords,hic=hic,stringsAsFactors = F)
  write.csv(hic_bed,file.path(inpdir,paste0(outname,"_",chr,".csv")),row.names = F)
}

extract_hic_parallel = function(sample_name,save_target=F,n_cores=6){
  print(sample_name)
  inpdir = file.path(wd,sample_name,"prediction/npy/")
  matrix_files = list.files(inpdir,pattern = ".npy")
  target_files = matrix_files[grep("_target.npy",matrix_files,invert = F)]
  prediction_files = matrix_files[grep("_target.npy",matrix_files,invert = T)]
  if(length(prediction_files) != n_windows){ print(paste0("WARNING: ",sample_name," does not have ",n_windows," matrices.")) }
  
  chromosomes = unique(gsub(pattern = "_.*.npy",replacement = "",x = prediction_files))
  for (chr in chromosomes){
    print(chr)
    prediction_files_chr = prediction_files[grep(paste0("^",chr,"_"),prediction_files,fixed = F)]
    target_files_chr = target_files[grep(paste0("^",chr,"_"),target_files,fixed = F)]
    
    print("Prediction HiC...")
    results = mclapply(prediction_files_chr,save_hic,inpdir=inpdir,window=window,mc.cores = n_cores)
    make_hic_bed(results,inpdir=inpdir,chr=chr)
    
    if(save_target){
      print("Target HiC...")
      results = mclapply(target_files_chr,save_hic,inpdir=inpdir,window=window,mc.cores = n_cores)
      make_hic_bed(results,inpdir=inpdir,chr=chr,outname="target_hic")
    }
  }
}

main = function(n_cores=6,extract_hic=T,compute_scores=T,save_target=F){
  if(extract_hic) { 
    print("Extracting HiC data...")
    lapply(sample_names,extract_hic_parallel,save_target=save_target,n_cores=n_cores) }
  if(compute_scores) { 
    print("Computing insulation scores...")
    lapply(sample_names,compute_scores_parallel,save_target=save_target,n_cores=n_cores) }
}

### RUN #########################
library(reticulate)
library(parallel)
np <- import("numpy")

wd=getwd()
window=2097152
n_windows=1359		# genome_size(chr1:chr22 + chrX) / window
resolution=8192		# 2097152 / 256
radius=500000		# insulation param
pseudocount=0.1
n_cores=8
extract_hic=T
compute_scores=T
save_target=T

#sample_names = list.dirs(wd,full.names = F,recursive = F)
sample_names = c("PANC1","GM12878","HCT116","HEPG2","IMR90","K562","MCF7")
main(n_cores=n_cores,extract_hic=extract_hic,compute_scores=compute_scores,save_target=save_target)
