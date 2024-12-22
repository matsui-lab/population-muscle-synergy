require(biwavelet)
require(parallel)
require(doParallel)
require(pracma)
require(psych)
require(GPArotation)
require(pcaMethods)
require(grDevices)
require(RColorBrewer)

data_check <- function(file_path){
  
  cat("Starting data validity check\n\n")
  
  dims <- matrix(NA,ncol=2,nrow=length(file_path))
  colnames(dims) <- c("nrow","ncol")
  rownames(dims) <- basename(file_path)
  
  data_type <- data_view <- vector("list",length(file_path))
  names(data_type) <- names(data_view) <- basename(file_path)
  
  for(i in seq_along(file_path)){
    dat <- fread(file_path[i],header=FALSE,data.table = FALSE)
    dims[i,] <- dim(dat)
    data_type[[i]] <- apply(dat,2,class)
    data_view[[i]] <- head(dat)
    
    cat("Extract information from: ",basename(file_path)[i],"\n")
  }
  cat("\n\n")
  
  unique_ncol <- unique(dims[,2])
  cat("Data sizes of columns are being checked ... ")
  if(length(unique_ncol)!=1){
    print(dims[,2,drop=FALSE])
    stop("Different size, please check and fix the data")
  }else{
    cat("O.K.\n")
  }
  cat("\n\n")
  
  
  class_check <- sapply(data_type,function(x)all(x=="numeric"))
  cat("Data type of each columns are being checked ...")
  if(any(!class_check)){
    print(do.call(rbind,data_type))
    stop(names(which(!class_check)),"include invalid type. Please confirm whether all data is numeric.")
  }else{
    cat("O.K. \n")
  }
  cat("\n\n")
  
  cat("Checking data completeness...\n\n")
  
  completeness_check <- matrix(NA,nrow = length(file_path),ncol= 3)
  colnames(completeness_check) <- c("na","inf","nan")
  rownames(completeness_check) <- basename(file_path)
  for(i in seq_along(file_path)){
    dat <- fread(file_path[i],header=FALSE,data.table = FALSE)
    
    nas <- apply(dat,2,is.na)
    infs <- apply(dat,2,is.infinite)
    nans <- apply(dat,2,is.nan)
    
    is_nas <- apply(nas,2,any)
    is_infs <- apply(infs,2,any)
    is_nans <- apply(nans,2,any)
    
    is_incomplete <- rbind(na = is_nas, inf = is_infs, nan = is_nans)
    is_incomplete_column <- apply(is_incomplete,2,any)
    completeness_check[i,] <- any(is_incomplete_column)
    
    cat("Checking ",basename(file_path)[i],"\n")
  }
  
  if(!all(completeness_check)){
    cat("\nO.K.\n")
  }else{
    warning("Some data set include incomplete data")
    print(completeness_check)
  }
  
  cat("\n\n")
  
  # unique_nrow <- unique(dims[,1])
  # cat("Data sizes of rows are being checked ... ")
  # if(length(unique_ncol)!=1){
  #   print(dims[,1,drop = FALSE])
  #   isTrim <- readline("Trimming into same row size? (yes or no): ")
  #   if(isTrim=="yes"){
  #     lower_nrow <- 1
  #     upper_nrow <- c(1,min(dims[,1]))
  #   }
  # }
  TRUE
  
}


filt <- function(obj,thres,method="btw",is_env = TRUE){
  for(i in seq_along(obj)){
    obj_i <- obj[[i]]
    
    freq <- obj_i$frequency
    
    if(method=="bandpass"){
      
      signal_filt <- apply(obj_i$signal,2,function(x)btw(x,freq,order = 4,bandpass = thres))
      
    }else if(method=="lowpass"){
      
      signal_filt <- apply(obj_i$signal,2,function(x)lowpass(x,freq,order = 2,thres = thres))
      
    }
    
    if(is_env){
      
      signal_filt <- apply(signal_filt,2,function(x)envelope(x,freq = freq,thres = thres))
    }
    
    
    obj_i$signal <- signal_filt
    obj[[i]] <- obj_i
  }
  return(obj)
}

detect_event <- function(obj = NULL,pad = 100,time = NULL,event = NULL,interactive = FALSE){
  
  if(is.null(obj) & (is.null(time) & is.null(event))){
    stop("Specifying any of arguments; 'obj' or ('time' / 'event')")  
  }
  
  #  output <- vector("list",length(obj))
  
  for(j in seq_along(obj)){
    
    obj_j <- obj[[j]]
    time <- obj_j$time
    event <- obj_j$event
    
    cat("Detecting gait cycle using event annotation table ... \n")
    mat_j <- matrix(0,nrow=nrow(event),ncol=ncol(event))
    
    for(i in 1:ncol(event)){
      
      x <- event[,i]
      xts <- cbind(time$time_stamp,x)
      
      #xxwt <- biwavelet::wt(xts,dt = 1/obj_j$frequency,dj = 5)
      xxwt <- biwavelet::wt(xts,dj = 5)
      cat("Select appropriate frequency (period) ... ")
      view_idx <- 1:50000
      
      filt <- t(xxwt$wave)
      filt_y <- as.matrix(Re(filt))
      filt_y[filt_y < 0] <- 0
      
      
      # if(interactive){
      #   matplot(filt_y[view_idx,],type="l",ylab="filtered wave",xlab="time",main=colnames(event)[i])
      #   lines(x[view_idx],col="#f18d08",lty=1, lwd = 3)
      #   legend("topright",
      #          legend = c(paste("period",1:ncol(filt),":",round(xxwt$period,2)),"origin"),
      #          col = c(1:ncol(filt),"#f18d08"),
      #          lty = c(1:ncol(filt),1),
      #          lwd = c(rep(1,ncol(filt)),3))
      #   
      #   cidx <- readline("select a period: ")
      #   
      #   cidx <- as.numeric(cidx)
      #   
      #   cat("Peak finding along with band filtered wave ... ")
      #   fp <- findpeaks(filt_y[,cidx])
      #   peak_idx <- fp[,2]
      #   cat("done \n")
      #   
      #   cat("Check the peak position ... \n")
      #   view_peak_idx <- peak_idx[peak_idx <= max(view_idx)]
      #   points(view_peak_idx,filt_y[view_peak_idx,cidx],pch=16,col="red")
      #   points(view_peak_idx,x[view_peak_idx],pch=16,col="blue")
      #   
      # }else{
      #   cidx <- 3
      #   cat("Peak finding along with band filtered wave ... ")
      #   fp <- findpeaks(filt_y[,cidx])
      #   peak_idx <- fp[,2]
      #   cat("done \n")
      # }
      # mat_j[fp[,2],i] <- 1
      
      
      if(interactive){
        matplot(filt_y[view_idx,],type="l",ylab="filtered wave",xlab="time",main=colnames(event)[i])
        lines(x[view_idx],col="#f18d08",lty=1, lwd = 3)
        legend("topright",
               legend = c(paste("period",1:ncol(filt),":",round(xxwt$period,2)),"origin"),
               col = c(1:ncol(filt),"#f18d08"),
               lty = c(1:ncol(filt),1),
               lwd = c(rep(1,ncol(filt)),3))
        
        cidx <- readline("select a period: ")
        
        cidx <- as.numeric(cidx)
        
        zeropoint <- filt_y[,cidx] == 0
        nextpoint <- which(zeropoint) + 1
        if(max(nextpoint) > nrow(filt_y)){
          zeropoint <- zeropoint[-length(zeropoint)]
          nextpoint <- nextpoint[-length(nextpoint)]
        }
        
        strikepoint <- nextpoint[filt_y[nextpoint,cidx] > 0] - 1 - pad
        strikepoint <- pmax(0,strikepoint)
        
        cat("done \n")
        
        cat("Check the peak position ... \n")
        
        view_strike_idx <- strikepoint[strikepoint <= max(view_idx)]
        points(view_strike_idx,filt_y[view_strike_idx,cidx],pch=16,col="red")
        points(view_strike_idx,x[view_strike_idx],pch=16,col="blue")
        
      }else{
        
        cidx <- 3
        cat("Peak finding along with band filtered wave ... ")
        #fp <- findpeaks(filt_y[,cidx])
        #peak_idx <- fp[,2]
        zeropoint <- filt_y[,3] == 0
        nextpoint <- which(zeropoint) + 1
        if(max(nextpoint) > nrow(filt_y)){
          zeropoint <- zeropoint[-length(zeropoint)]
          nextpoint <- nextpoint[-length(nextpoint)]
        }
        
        strikepoint <- nextpoint[filt_y[nextpoint,3] > 0] - 1 - pad
        strikepoint <- pmax(0,strikepoint)
        
        cat("done \n")
      }
      
      mat_j[strikepoint,i] <- 1
      cat("Finish event",i,"\n\n")
    }
    
    colnames(mat_j) <- paste(".event",colnames(event),sep="_")
    mat_j <- cbind(time,event,mat_j)
    obj[[j]]$event <- mat_j
    
    cat("Finish detecting gait cycles of","file",j,"!\n")
    
  }
  obj
}

time_normalize = function(signal,time,event = NULL,n){
  
  xx <- signal
  tt <- time
  
  
  approx_out <- lapply(1:ncol(xx),
                       function(ii)approx(tt,
                                          xx[,ii],
                                          method = "linear",
                                          n = n))
  
  
  ttout <- sapply(approx_out,function(x)x$x)[,1]
  xxout <- sapply(approx_out,function(x)x$y)
  
  rownames(xxout) <- 1:nrow(xxout)
  colnames(xxout) <- colnames(xx)
  
  if(!is.null(event)){
    zz <- event[,!(colnames(event) %in% c("time_stamp","time_step"))]
    
    approx_z <- lapply(1:ncol(zz),
                       function(ii)approx(tt,
                                          zz[,ii],
                                          method = "linear",
                                          n = n))
    zzout <- sapply(approx_z,function(x)x$y)
    colnames(zzout) <- colnames(zz)
    
    new_time <- data.frame(time_stamp = seq(0,100,length.out = length(ttout)),time_step = seq_along(ttout),time_stamp_raw = ttout)
    zzout <- cbind(new_time,zzout)
    
    output <- list(signal_approx = xxout, time_approx = new_time,event_approx = zzout)
    
  }else{
    new_time <- data.frame(time_stamp = seq(0,100,length.out = length(ttout)),time_step = seq_along(ttout),time_stamp_raw = ttout)
    output <- list(signal_approx = xxout, time_approx = new_time)
  }
  
  return(output)
}

# create_gait_epoch <- function(signal,freqency = 1000,event_anno,number_of_cycle = 5){
# 
#   cat("Segmentation of signals into gait cycle ... \n")
#   cat("Check number of event ... ")
#   
#   is_event <- event_anno[,grepl(".event",colnames(event_anno))]
#   n_event0 <- apply(is_event,2,sum)
#   
#   if(n_event0[1] != n_event0[2]){
#     warning("Number of event is different. Automatically select smaller one.\n")
#     n_event <- min(n_event0)
#     cat("Trimming event annotation adjusting to smaller number of event ... ")
#     trim_idx <- apply(is_event,2,function(x)cumsum(x) <= n_event)
#     is_event <- is_event[trim_idx,]
#     event_anno <- event_anno[trim_idx,]
#     cat("done\n")
#   }else{
#     cat("Number of event is",n_event0,"\n")
#   }
#   
#   event_min_indice <- apply(is_event,2,function(xx)min(which(xx==1)))
#   event_min_idx <- which.min(event_min_indice)
#   
#   cat("Specified number of gait cycle per epoch is: ",number_of_cycle,"\n")
#   
#   epoch_idx <- which(is_event[,event_min_idx] == 1)
#   if(number_of_cycle > 1){
#     sti <- seq(1,length(epoch_idx),by = number_of_cycle - 1)
#     edi <- seq(number_of_cycle,length(epoch_idx),by = number_of_cycle - 1)
#   }else{
#     sti <- seq(1,length(epoch_idx),by = number_of_cycle)
#     edi <- seq(1 + number_of_cycle,length(epoch_idx),by = number_of_cycle)
#   }
#   if(length(sti)!=length(edi)){sti <- sti[-length(sti)]}
#   st <- epoch_idx[sti]
#   ed <- epoch_idx[edi]
#   
#   n_epoch <- length(st)
#   
#   cat("Number of epoch is",n_epoch,"\n")
#   
#   cat("Segmenting ... ")
#   
#   epoch_li <- vector("list",length(st))
#   for(i in seq_along(epoch_li)){
#     epoch_i <- signal[st[i]:ed[i],]
#     time_i <- time[st[i]:ed[i],]
#     event_i <- event[st[i]:ed[i],]
#     
#     time_i$time_step <- 1:nrow(time_i)
#     time_i$time_stamp <- time_i$time_stamp - time_i$time_stamp[1] + 1/freqency
#     time_i$time_stamp = seq(1 / frequency, nrow(time_i) / frequency,by = 1 / frequency)
#     
#     epoch_li[[i]]$signal <- epoch_i
#     epoch_li[[i]]$time <- time_i
#     epoch_li[[i]]$event <- event_i
#   }
#   cat("done\n")
#   epoch_li
# }


create_gait_epoch <- function(obj,freq = 2000,number_of_cycle = 1,epoch_size = 1000){
  
  cat("Segmentation of signals into gait cycle ... \n")
  cat("Check number of event ... ")
  
  for(j in seq_along(obj)){
    
    obj_j <- obj[[j]]
    signal <- obj_j$signal
    gait_time <- obj_j$time
    event_anno <- obj_j$event
    
    is_event <- event_anno[,grepl(".event",colnames(event_anno))]
    n_event0 <- apply(is_event,2,sum)
    
    #if(n_event0[1] != n_event0[2]){
    if(length(unique(n_event0))!=1){
      warning("Number of event is different. Automatically select smaller one.\n")
      n_event <- min(n_event0)
      cat("Trimming event annotation adjusting to smaller number of event ... ")
      trim_idx <- apply(is_event,2,function(x)cumsum(x) <= n_event)
      is_event <- is_event[trim_idx,]
      event_anno <- event_anno[trim_idx,]
      cat("done\n")
    }else{
      cat("Number of event is",n_event0,"\n")
    }
    
    event_min_indice <- apply(is_event,2,function(xx)min(which(xx==1)))
    event_min_idx <- which.min(event_min_indice)
    
    cat("Specified number of gait cycle per epoch is: ",number_of_cycle,"\n")
    
    # 一つのイベントを想定 --> イベントそれぞれに要変更
    epoch_li <- vector("list",ncol(is_event))
    
    
    for(evi in 1:ncol(is_event)){
      epoch_idx <- which(is_event[,evi] == 1)
      
      if(number_of_cycle > 1){
        sti <- seq(1,length(epoch_idx),by = number_of_cycle - 1)
        edi <- seq(number_of_cycle,length(epoch_idx),by = number_of_cycle - 1)
      }else{
        sti <- seq(1,length(epoch_idx),by = number_of_cycle)
        edi <- seq(1 + number_of_cycle,length(epoch_idx),by = number_of_cycle)
      }
      
      if(length(sti)!=length(edi)){sti <- sti[-length(sti)]}
      st <- epoch_idx[sti]
      ed <- epoch_idx[edi]
      
      n_epoch <- length(st)
      
      cat("Number of epoch is",n_epoch,"\n")
      
      cat("Segmenting ... ")
      
      n1 <- epoch_size
      n2 <- ncol(signal) + ncol(event_anno) + 1
      n3 <- length(st)
      
      signal_array <- array(NA,c(n1,n2,n3),
                            dimnames = list(1:n1 ,
                                            c("time",colnames(signal),
                                              colnames(event_anno)),
                                            paste("cycle",1:n3,sep="_")))
      
      
      for(i in 1:n3){
        
        epoch_i <- signal[st[i]:ed[i],]
        time_i <- gait_time[st[i]:ed[i],]
        event_i <- event_anno[st[i]:ed[i],]
        
        time_i$time_step <- 1:nrow(time_i)
        time_i$time_stamp <- time_i$time_stamp - time_i$time_stamp[1] + 1/freq
        time_i$time_stamp = seq(1 / freq, nrow(time_i) / freq,by = 1 / freq)
        
        norm_obj <- time_normalize(signal = epoch_i,
                                   time = time_i$time_stamp,
                                   event = event_i,
                                   n = epoch_size)
        
        epoch_i <- norm_obj$signal_approx
        time_i <- norm_obj$time_approx
        event_i <- norm_obj$event_approx
        
        
        signal_array[,,i] <- as.matrix(cbind(time_i$time_stamp_raw,epoch_i,event_i[,!colnames(event_i)=="time_stamp_raw"]))
        
      }
      
      
      epoch_li[[evi]] <- signal_array
      
    }
    
    names(epoch_li) <- gsub(".event_","",colnames(is_event))
    obj_j$epoch <- epoch_li
    
    obj[[j]] <- obj_j
    cat(j,"th","subject","done\n")
    
  }
  
  obj
  
}


create_event_epoch <- function(obj = NULL,
                               signal = NULL,
                               frequency = NULL,
                               event_anno = NULL,
                               epoch_event = NULL,
                               epoch_window_lower = 600,
                               epoch_window_upper = 1600){
  
  if(is.null(obj) & (is.null(signal) & is.null(signal) & is.null(event_anno))){
    stop("Specifying any of arguments; 'obj' or ('signal' / 'frequency' / 'event_anno')")  
  }
  
  for(j in seq_along(obj)){
    cat("Processing",j,"th file \n\n")
    obj_j <- obj[[j]]
    
    if(!is.null(obj)){
      signal <- obj_j$signal
      frequency <- obj_j$frequency
      event_anno <- obj_j$event
    }
    
    epoch_window_len <- epoch_window_lower + epoch_window_upper
    
    cat("Specified epoch window size is ... [",
        epoch_window_lower,",",epoch_window_upper,"]\n")
    cat("Epoch window size is",epoch_window_len,"\n")
    
    cat("Check number of event ... ")
    #is_event <- event_anno[,grepl(".event",colnames(event_anno))]
    select_colum_event_anno <- match(paste0(".event_",epoch_event),colnames(event_anno))
    is_event <- event_anno[,select_colum_event_anno]
    n_event0 <- apply(is_event,2,sum)
    
    if(n_event0[1] != n_event0[2]){
      warning("Number of event is different. Automatically select smaller one.\n")
      replace_event <- which.max(n_event0)
      replace_idx <- (min(n_event0) + 1):max(n_event0)
      cat("Trimming event annotation adjusting to smaller number of event ... ")
      trim_idx <- which(is_event[,replace_event] == 1)[replace_idx]
      is_event[trim_idx,replace_event] <- 0
      event_anno[,select_colum_event_anno][trim_idx,replace_event] <- 0
      cat("done\n")
    }else{
      cat("Number of event is",n_event0,"\n")
    }
    
    is_event_idx <- apply(is_event,2,function(x)which(x==1))
    
    
    cat("Calibrating data to build epochs ... ")
    
    # Left side calibration
    calib_check1 <- apply(is_event_idx,2,function(x)x - epoch_window_lower > 0)
    calib_num1 <- sum(apply(calib_check1,1,function(x)!all(x)))
    calib_idx1 <- which(apply(calib_check1,1,function(x)!all(x)))
    
    # Right side calibration 
    calib_check2 <- apply(is_event_idx,2,function(x)x + epoch_window_upper < nrow(signal))
    calib_num2 <- sum(apply(calib_check2,1,function(x)!all(x)))
    calib_idx2 <- which(apply(calib_check2,1,function(x)!all(x)))
    
    
    calib_idx <- c(calib_idx1,calib_idx2)
    
    if(length(calib_idx)!=0){
      is_event_idx_calib <- is_event_idx[-calib_idx,,drop = FALSE]
    }else{
      cat("No calibration is needed ... \n\n")
      is_event_idx_calib <- is_event_idx
    }
    
    cat("done\n")
    cat(calib_num1 + calib_num2,"event(s) is calibrated from data\n")
    
    
    cat("Check the specified epoch window size is enough\n")
    
    cat("Counting number of time points within \n[epoch_window_lower, epoch_window_upper] ... \n")
    
    event_lag_idx <- apply(is_event_idx_calib,2,diff)
    epoch_window_len0 <- apply(event_lag_idx,2,range)
    rownames(epoch_window_len0) <- c("Minimum of interval","Maximum of interval")
    
    
    cat("Check number of time points is enough to construct specified epoch window size ... \n")
    is_enough_epoch_len <- apply(event_lag_idx,2,
                                 function(x)x > epoch_window_len)
    
    remove_epoch <- apply(is_enough_epoch_len,2,function(x)!any(x))
    
    if(any(remove_epoch)){
      warning("Specified several epoch window length is not enough ... ")
      remove_idx <- which(remove_epoch)
      
      for(i in seq_along(remove_idx)){
        idx <- remove_idx[i]
        target_column <- gsub(".event_","",colnames(is_enough_epoch_len)[idx])
        msg <- paste("In",target_column,"following indices are not enough time points:")
        msg <- paste(which(is_enough_epoch_len[,idx]),sep=" ")
        stop(msg)
      }
      
    }else{
      cat("Enough time points to construct the epoch! \n")
    }
    
    event_name <- gsub(".event_","",colnames(is_enough_epoch_len))
    cat("Now creating epochs for",event_name,"... \n\n")
    
    epoch_li <- vector("list",ncol(is_event))
    names(epoch_li) <- event_name
    
    n1 <- epoch_window_len + 1
    
    time_stamp <- seq((-1)*epoch_window_lower * 1000/frequency,
                      epoch_window_upper * 1000/frequency,by = 1000/frequency)
    
    for(i in 1:ncol(is_event_idx_calib)){
      
      cat("Processing event of",event_name[i],"...")
      st <- is_event_idx_calib[,i] - epoch_window_lower
      ed <- is_event_idx_calib[,i] + epoch_window_upper
      epoch_i <- lapply(seq_along(st),function(ii)signal[st[ii]:ed[ii],])
      event_i <- lapply(seq_along(st),function(ii)event_anno[st[ii]:ed[ii],])
      
      n2 <- ncol(signal) + ncol(event_i[[1]])
      n3 <- length(epoch_i)
      
      signal_array <- array(NA,c(n1,n2 + 1,n3),
                            dimnames = list(1:n1 ,c("time",colnames(signal),colnames(event_i[[1]])),
                                            paste("cycle",1:n3,sep="_")))
      
      for(k in seq_along(epoch_i)){
        signal_array[,,k] <- as.matrix(cbind(time_epoch = time_stamp,epoch_i[[k]],event_i[[k]]))
      }
      
      epoch_li[[i]] <- signal_array
      cat("done\n\n")
    }
    
    obj[[j]]$epoch <- epoch_li
    
    cat("Processing",j,"th file finished \n\n")
    
  }
  obj
}

frequency_check <- function(ch_annotation){
  cat("Frequency check ... \n\n")
  
  ch_freq <- unique(ch_annotation$ch_frequency)
  if(length(ch_freq)!=1){
    cat("There are different sampling rate channels ... \n")
    
    for(i in seq_along(ch_freq))cat(i,":",unique(ch_freq[i]),"Hz\n")
    
    freq_sel <- as.numeric(readline("Which frequency do you want to tune to ? (Select no.):"))
    target_freq <- ch_freq[freq_sel]
    src_freq <- ch_freq[-freq_sel]
    
    idx <- ch_annotation$ch_frequency != target_freq
    
    return(list(index = idx, target_freq = target_freq, src_freq = src_freq))
    
  }else{
    cat("All the sampling rate is the same. No adjusting will be required ... \n")
    return(NULL)
  }
  
}

plot_epoch_wt <- function(obj,
                          subject_index = 1,
                          event_type = 1,
                          region_index = 1,
                          epoch_index = 1,
                          interactive = TRUE,
                          return_val = FALSE,
                          save = NULL){
  
  if(interactive){
    cat("subject_index:\n",seq_along(obj))
    subject_index <- as.numeric(readline("Select subject index:"))
  }
  
  i <- subject_index
  
  # select i-th object (i-th file/subject)
  obj_i <- obj[[i]]
  
  if(interactive){
    cat("event_type:\n",paste0(seq_along(obj_i$epoch),":",names(obj_i$epoch),"\n"))
    event_type <- as.numeric(readline("Select event_type:"))
    
    cat("region_index:\n",paste0(1:ncol(obj_i$signal),":",colnames(obj_i$signal),"\n"))
    region_index <- as.numeric(readline("Select region index:"))
  }
  
  j <- event_type
  k <- region_index
  
  if(interactive){
    cat("epoch_index:\n","from",1,"to",dim(obj_i$epoch[[j]])[3],"\n")
    epoch_index <- as.numeric(readline("Select epoch index:"))
  }
  
  l <- epoch_index
  
  raw_array <- obj_i$epoch[[j]]
  epoch_time <- raw_array[,1,1]
  event_array <- raw_array[,grep(".event",dimnames(raw_array)[[2]]),]
  signal_array <- raw_array[,colnames(obj_i$signal),]
  
  # select k-th region
  signal_array_k <- signal_array[,k,]
  
  xt <- cbind(epoch_time/1000,signal_array_k[,l])
  
  wtx <- biwavelet::wt(xt)
  
  if(!is.null(save)){
    pdf(paste(paste("plot","obj",i,j,k,l,sep="_"),"pdf",sep="."))
  }
  
  layout_mat <- matrix(c(1,2,1,3,1,4),ncol=2,byrow = TRUE)
  #  layout_mat <- matrix(c(1,1,1,2,2,2,3,4,5),ncol=3)
  layout(layout_mat)
  
  biwavelet::plot.biwavelet(wtx,xaxt = "n",yaxt = "n",xlab = "Epoch time",ylab = "Frequency(Hz)")
  title("Wavelet time-frequcency")
  
  per2freq <- cbind(index = seq_along(wtx$period), period = wtx$period,freq = 1/wtx$period)
  
  axis.locs <- log2(1/c(seq(1000,1,by=-10),seq(9,0,by = -1)))
  yticklab <- format(c(seq(1000,1,by=-10),seq(9,0,by = -1)),digits = 1)
  axis(side = 2, at = axis.locs, labels = yticklab)
  
  axis.locs <- axTicks(side = 1)
  xticklab <- format(axis.locs*1000,digits = 1)
  axis(side = 1, at = axis.locs, labels = xticklab)
  
  event_idx <- matrix(c(1,2,3,4),byrow = TRUE,ncol = 2)
  
  abline(v=epoch_time[event_array[,event_idx[j,1],l]==1][1]/1000,lwd=3,col="grey",lty=3)
  abline(v=epoch_time[event_array[,event_idx[j,2],l]==1][1]/1000,lwd=3,col="grey",lty=3)
  
  title1 <- paste("original singnals for each epoch\n","subject:",subject_index,"region:",region_name,"event type",event_type_name)
  title2 <- paste("wavelet filtered signals for each epoch\n","subject:",subject_index,"region:",region_name,"event type",event_type_name)
  title3 <- paste("wavelet coherence for each epoch\n","subject:",subject_index,"region:",region_name,"event type",event_type_name)
  
  matplot(epoch_time,xt[,-1],type="l",ylab="origin",xlab="Epoch time",main=title1)
  grid()
  abline(v=epoch_time[event_array[,event_idx[j,1],l]==1][1]/1000,lwd=3,col="grey",lty=3)
  abline(v=epoch_time[event_array[,event_idx[j,2],l]==1][1]/1000,lwd=3,col="grey",lty=3)
  
  matplot(epoch_time,t(Re(wtx$wave)),type="l",ylab="wave",xlab="Epoch time",main=title2)
  grid()
  abline(v=epoch_time[event_array[,event_idx[j,1],l]==1][1]/1000,lwd=3,col="grey",lty=3)
  abline(v=epoch_time[event_array[,event_idx[j,2],l]==1][1]/1000,lwd=3,col="grey",lty=3)
  
  matplot(epoch_time,t(wtx$power),type="l",ylab="power",xlab="Epoch time",main=title3)
  grid()
  abline(v=epoch_time[event_array[,event_idx[j,1],l]==1][1]/1000,lwd=3,col="grey",lty=3)
  abline(v=epoch_time[event_array[,event_idx[j,2],l]==1][1]/1000,lwd=3,col="grey",lty=3)
  
  if(!is.null(save)){
    dev.off()
  }
  if(return_val){
    return(list(wt = wtx, per2freq = per2freq))
  }
}


diag_epoch <- function(obj,
                       subject_index = 1,
                       event_type = 1,
                       region_index_1 = 1,
                       region_index_2 = 2,
                       epoch_index = 1,
                       interactive = TRUE,
                       return_val = FALSE,
                       save = NULL){
  
  if(interactive){
    cat("subject_index:\n",seq_along(obj))
    subject_index <- as.numeric(readline("Select subject index:"))
  }
  
  i <- subject_index
  
  # select i-th object (i-th file/subject)
  obj_i <- obj[[i]]
  
  if(interactive){
    cat("event_type:\n",paste0(seq_along(obj_i$epoch),":",names(obj_i$epoch),"\n"))
    event_type <- as.numeric(readline("Select event_type:"))
    
    cat("region_index:\n",paste0(1:ncol(obj_i$signal),":",colnames(obj_i$signal),"\n"))
    region_index_1 <- as.numeric(readline("Select region index_1:"))
    region_index_2 <- as.numeric(readline("Select region index_2:"))
  }
  
  j <- event_type
  k_1 <- region_index_1
  k_2 <- region_index_2
  
  if(interactive){
    cat("epoch_index:\n","from",1,"to",dim(obj_i$epoch[[j]])[3],"\n")
    epoch_index <- as.numeric(readline("Select epoch index:"))
  }
  
  l <- epoch_index
  
  raw_array <- obj_i$epoch[[j]]
  epoch_time <- raw_array[,1,1]
  event_array <- raw_array[,grep(".event",dimnames(raw_array)[[2]]),]
  signal_array <- raw_array[,colnames(obj_i$signal),]
  
  # select k-th region
  signal_array_k_1 <- signal_array[,k_1,]
  signal_array_k_2 <- signal_array[,k_2,]
  
  region_name_1 <- names(obj_i$signal)[k_1]
  region_name_2 <- names(obj_i$signal)[k_2]
  
  xt_1 <- cbind(epoch_time/1000,signal_array_k_1[,l])
  xt_2 <- cbind(epoch_time/1000,signal_array_k_2[,l])
  
  wtx1 <- biwavelet::wt(xt_1)
  wtx2 <- biwavelet::wt(xt_2)
  xwtx <- biwavelet::xwt(xt_1,xt_2)
  wtcx <- biwavelet::wtc(xt_1,xt_2,nrands = 0)
  per2freq <- cbind(index = seq_along(xwtx$period), period = xwtx$period,freq = 1/xwtx$period)
  
  if(!is.null(save)){
    pdf(paste(paste("plot","obj",i,j,k,l,sep="_"),"pdf",sep="."))
  }
  
  event_idx <- matrix(c(1,2,3,4),byrow = TRUE,ncol = 2)
  p <- par(cex = 0.5)
  split.screen(c(3,3))
  
  screen(1)
  biwavelet::plot.biwavelet(wtx1,xaxt = "n",yaxt = "n",xlab = "Epoch time",ylab = "Frequency(Hz)")
  title(paste("Scalogram of wavelet of",region_name_1,"\nin subject",i,"during epoch",epoch_index))
  
  axis.locs <- log2(1/c(seq(1000,1,by=-10),seq(9,0,by = -1)))
  yticklab <- format(c(seq(1000,1,by=-10),seq(9,0,by = -1)),digits = 1)
  axis(side = 2, at = axis.locs, labels = yticklab)
  
  axis.locs <- axTicks(side = 1)
  xticklab <- format(axis.locs*1000,digits = 1)
  axis(side = 1, at = axis.locs, labels = xticklab)
  
  abline(v=epoch_time[event_array[,event_idx[j,1],l]==1][1]/1000,lwd=3,col="grey",lty=3)
  abline(v=epoch_time[event_array[,event_idx[j,2],l]==1][1]/1000,lwd=3,col="grey",lty=3)
  
  
  screen(2)
  biwavelet::plot.biwavelet(wtx2,xaxt = "n",yaxt = "n",xlab = "Epoch time",ylab = "Frequency(Hz)")
  title(paste("Scalogram of wavelet of",region_name_2,"\nin subject",subject_index,"during epoch",epoch_index))
  
  axis.locs <- log2(1/c(seq(1000,1,by=-10),seq(9,0,by = -1)))
  yticklab <- format(c(seq(1000,1,by=-10),seq(9,0,by = -1)),digits = 1)
  axis(side = 2, at = axis.locs, labels = yticklab)
  
  axis.locs <- axTicks(side = 1)
  xticklab <- format(axis.locs*1000,digits = 1)
  axis(side = 1, at = axis.locs, labels = xticklab)
  
  abline(v=epoch_time[event_array[,event_idx[j,1],l]==1][1]/1000,lwd=3,col="grey",lty=3)
  abline(v=epoch_time[event_array[,event_idx[j,2],l]==1][1]/1000,lwd=3,col="grey",lty=3)
  
  
  
  logxwtx <- xwtx
  logxwtx$power <- t(apply(xwtx$power.corr,1,function(x)scale(x)^2))
  screen(4)
  biwavelet::plot.biwavelet(logxwtx,xaxt = "n",type="power",plot.phase = TRUE,
                            yaxt = "n",
                            xlab = "Epoch time",
                            ylab = "Frequency(Hz)",zlim=c(-5,5))
  title(paste("Scalogram of cross-wavelet between\n",region_name_1,"and", region_name_2,"in subject",subject_index,"during epoch",epoch_index))
  
  
  axis.locs <- log2(1/c(seq(1000,1,by=-10),seq(9,0,by = -1)))
  yticklab <- format(c(seq(1000,1,by=-10),seq(9,0,by = -1)),digits = 1)
  axis(side = 2, at = axis.locs, labels = yticklab)
  
  axis.locs <- axTicks(side = 1)
  xticklab <- format(axis.locs*1000,digits = 1)
  axis(side = 1, at = axis.locs, labels = xticklab)
  
  abline(v=epoch_time[event_array[,event_idx[j,1],l]==1][1]/1000,lwd=3,col="grey",lty=3)
  abline(v=epoch_time[event_array[,event_idx[j,2],l]==1][1]/1000,lwd=3,col="grey",lty=3)
  
  
  
  screen(5)
  biwavelet::plot.biwavelet(wtcx,xaxt = "n",plot.phase = TRUE,
                            yaxt = "n",
                            xlab = "Epoch time",
                            ylab = "Frequency(Hz)")
  
  title(paste("Scalogram of wavelet coherence between\n",region_name_1,"and", region_name_2,"in subject",subject_index,"during epoch",epoch_index))
  
  axis.locs <- log2(1/c(seq(1000,1,by=-10),seq(9,0,by = -1)))
  yticklab <- format(c(seq(1000,1,by=-10),seq(9,0,by = -1)),digits = 1)
  axis(side = 2, at = axis.locs, labels = yticklab)
  
  axis.locs <- axTicks(side = 1)
  xticklab <- format(axis.locs*1000,digits = 1)
  axis(side = 1, at = axis.locs, labels = xticklab)
  
  abline(v=epoch_time[event_array[,event_idx[j,1],l]==1][1]/1000,lwd=3,col="grey",lty=3)
  abline(v=epoch_time[event_array[,event_idx[j,2],l]==1][1]/1000,lwd=3,col="grey",lty=3)
  
  
  title1 <- paste("Cross-wavelet power between",region_name_1,"and",region_name_2,"\nin subject:",subject_index,"during epoch",epoch_index)
  title2 <- paste("Wavelet coherence between",region_name_1,"and",region_name_2,"\nin subject:",subject_index,"during epoch",epoch_index)
  
  title3 <- paste("Original signal of\n",region_name_1,"and",region_name_2,"in subject:",subject_index,"during epoch",epoch_index)
  title4 <- paste("Wavelet spectrum of\n",region_name_1,"in subject:",subject_index,"during epoch",epoch_index)
  title5 <- paste("Wavelet spectrum of\n",region_name_2,"in subject:",subject_index,"during epoch",epoch_index)
  
  par(p)
  
  p <- par(cex=0.5,oma=c(0,0,3,0))
  screen(7)
  matplot(epoch_time,t(xwtx$power.corr),type="l",ylab="wave",xlab="Epoch time",main=title1)
  grid()
  abline(v=epoch_time[event_array[,event_idx[j,1],l]==1][1],lwd=3,col="grey",lty=3)
  abline(v=epoch_time[event_array[,event_idx[j,2],l]==1][1],lwd=3,col="grey",lty=3)
  
  screen(8)
  matplot(epoch_time,t(wtcx$rsq),type="l",ylab="power",xlab="Epoch time",main=title2)
  grid()
  abline(v=epoch_time[event_array[,event_idx[j,1],l]==1][1],lwd=3,col="grey",lty=3)
  abline(v=epoch_time[event_array[,event_idx[j,2],l]==1][1],lwd=3,col="grey",lty=3)
  
  
  screen(3)
  matplot(epoch_time,xt_1[,-1],type="l",ylab="vol",xlab="Epoch time",main=title3)
  matplot(epoch_time,xt_2[,-1],type="l",ylab="vol",xlab="Epoch time",col = 3,lty = 1,main=title1,add = TRUE)
  legend("topright",legend = c(region_name_1,region_name_2),lty=c(1,1),col=c(1,3))
  grid()
  abline(v=epoch_time[event_array[,event_idx[j,1],l]==1][1],lwd=3,col="grey",lty=3)
  abline(v=epoch_time[event_array[,event_idx[j,2],l]==1][1],lwd=3,col="grey",lty=3)
  
  screen(6)
  matplot(epoch_time,t(wtx1$power.corr),type="l",ylab="power",xlab="Epoch time",main=title4)
  grid()
  abline(v=epoch_time[event_array[,event_idx[j,1],l]==1][1],lwd=3,col="grey",lty=3)
  abline(v=epoch_time[event_array[,event_idx[j,2],l]==1][1],lwd=3,col="grey",lty=3)
  
  screen(9)
  matplot(epoch_time,t(wtx2$power.corr),type="l",ylab="power",xlab="Epoch time",main=title5)
  grid()
  abline(v=epoch_time[event_array[,event_idx[j,1],l]==1][1],lwd=3,col="grey",lty=3)
  abline(v=epoch_time[event_array[,event_idx[j,2],l]==1][1],lwd=3,col="grey",lty=3)
  par(p)
  
  close.screen(all.screens = TRUE)
  
  
  
  if(!is.null(save)){
    dev.off()
  }
  if(return_val){
    return(list(wt = wtx, per2freq = per2freq))
  }
}


diag_simple <- function(obj,
                        subject_index = 1,
                        event_type = 1,
                        region_index_1 = 1,
                        region_index_2 = 2,
                        freq = c(5,10,20,30,50),
                        interactive = TRUE,
                        return_val = FALSE,
                        save = NULL){
  
  if(interactive){
    cat("subject_index:\n",seq_along(obj))
    subject_index <- as.numeric(readline("Select subject index:"))
  }
  
  i <- subject_index
  
  # select i-th object (i-th file/subject)
  obj_i <- obj[[i]]
  
  if(interactive){
    cat("event_type:\n",paste0(seq_along(obj_i$epoch),":",names(obj_i$epoch),"\n"))
    event_type <- as.numeric(readline("Select event_type:"))
    
    cat("region_index:\n",paste0(1:ncol(obj_i$signal),":",colnames(obj_i$signal),"\n"))
    region_index_1 <- as.numeric(readline("Select region index_1:"))
    region_index_2 <- as.numeric(readline("Select region index_2:"))
  }
  
  j <- event_type
  k_1 <- region_index_1
  k_2 <- region_index_2
  
  # if(interactive){
  #   cat("epoch_index:\n","from",1,"to",dim(obj_i$epoch[[j]])[3],"\n")
  #   epoch_index <- as.numeric(readline("Select epoch index:"))
  # }
  # 
  # l <- epoch_index
  # 
  
  raw_array <- obj_i$epoch[[j]]
  epoch_time <- raw_array[,1]
  event_array <- raw_array[,grep(".event",dimnames(raw_array)[[2]])]
  signal_array <- raw_array[,colnames(obj_i$signal)]
  
  # select k-th region
  signal_array_k_1 <- signal_array[,k_1]
  signal_array_k_2 <- signal_array[,k_2]
  
  region_name_1 <- names(obj_i$signal)[k_1]
  region_name_2 <- names(obj_i$signal)[k_2]
  
  xt_1 <- cbind(epoch_time/1000,signal_array_k_1)
  xt_2 <- cbind(epoch_time/1000,signal_array_k_2)
  
  wtx1 <- biwavelet::wt(xt_1)
  wtx2 <- biwavelet::wt(xt_2)
  xwtx <- biwavelet::xwt(xt_1,xt_2)
  wtcx <- biwavelet::wtc(xt_1,xt_2,nrands = 0)
  per2freq <- cbind(index = seq_along(xwtx$period), period = xwtx$period,freq = 1/xwtx$period)
  
  cat("Approximating specified frequency to obtain wavelet frequency ... \n")
  
  idx <- per2freq[sapply(freq,function(xx)which.min((per2freq[,3] - xx)^2)),1]
  
  cat("The closest wavelet frequency was: ",paste(round(per2freq[idx,3],2),"Hz (in terms of period:",round(per2freq[idx,2],2),")"),sep="\n")
  
  if(!is.null(save)){
    pdf(paste(paste("plot","obj",i,j,k,l,sep="_"),"pdf",sep="."))
  }
  
  event_idx <- matrix(c(1,2,3,4),byrow = TRUE,ncol = 2)
  
  p <- par(cex = 0.5)
  split.screen(c(3,3))
  
  screen(1)
  biwavelet::plot.biwavelet(wtx1,xaxt = "n",yaxt = "n",xlab = "Epoch time",ylab = "Frequency(Hz)")
  title(paste("Scalogram of wavelet of",region_name_1,"\nin subject",subject_index))
  
  axis.locs <- log2(1/c(seq(1000,1,by=-10),seq(9,0,by = -1)))
  yticklab <- format(c(seq(1000,1,by=-10),seq(9,0,by = -1)),digits = 1)
  axis(side = 2, at = axis.locs, labels = yticklab)
  
  axis.locs <- axTicks(side = 1)
  xticklab <- format(axis.locs*1000,digits = 1)
  axis(side = 1, at = axis.locs, labels = xticklab)
  
  abline(v=epoch_time[event_array[,event_idx[j,1]]==1][1]/1000,lwd=3,col="grey",lty=3)
  abline(v=median(epoch_time[event_array[,event_idx[j,2]]!=0])/1000,lwd=3,col="grey",lty=3)
  
  
  screen(2)
  biwavelet::plot.biwavelet(wtx2,xaxt = "n",yaxt = "n",xlab = "Epoch time",ylab = "Frequency(Hz)")
  title(paste("Scalogram of wavelet of",region_name_2,"\nin subject",subject_index))
  
  axis.locs <- log2(1/c(seq(1000,1,by=-10),seq(9,0,by = -1)))
  yticklab <- format(c(seq(1000,1,by=-10),seq(9,0,by = -1)),digits = 1)
  axis(side = 2, at = axis.locs, labels = yticklab)
  
  axis.locs <- axTicks(side = 1)
  xticklab <- format(axis.locs*1000,digits = 1)
  axis(side = 1, at = axis.locs, labels = xticklab)
  
  abline(v=epoch_time[event_array[,event_idx[j,1]]==1][1]/1000,lwd=3,col="grey",lty=3)
  abline(v=median(epoch_time[event_array[,event_idx[j,2]]!=0])/1000,lwd=3,col="grey",lty=3)
  
  
  logxwtx <- xwtx
  logxwtx$power <- t(apply(xwtx$power.corr,1,function(x)scale(x)^2))
  screen(4)
  biwavelet::plot.biwavelet(logxwtx,xaxt = "n",type="power",plot.phase = TRUE,
                            yaxt = "n",
                            xlab = "Epoch time",
                            ylab = "Frequency(Hz)",zlim=c(-5,5))
  title(paste("Scalogram of z-transformed cross-wavelet between\n",region_name_1,"and", region_name_2,"in subject",subject_index))
  
  
  axis.locs <- log2(1/c(seq(1000,1,by=-10),seq(9,0,by = -1)))
  yticklab <- format(c(seq(1000,1,by=-10),seq(9,0,by = -1)),digits = 1)
  axis(side = 2, at = axis.locs, labels = yticklab)
  
  axis.locs <- axTicks(side = 1)
  xticklab <- format(axis.locs*1000,digits = 1)
  axis(side = 1, at = axis.locs, labels = xticklab)
  
  abline(v=epoch_time[event_array[,event_idx[j,1]]==1][1]/1000,lwd=3,col="grey",lty=3)
  abline(v=median(epoch_time[event_array[,event_idx[j,2]]!=0])/1000,lwd=3,col="grey",lty=3)
  
  
  screen(5)
  biwavelet::plot.biwavelet(wtcx,xaxt = "n",plot.phase = TRUE,
                            yaxt = "n",
                            xlab = "Epoch time",
                            ylab = "Frequency(Hz)")
  
  title(paste("Scalogram of wavelet coherence between\n",region_name_1,"and", region_name_2,"in subject",subject_index))
  
  axis.locs <- log2(1/c(seq(1000,1,by=-10),seq(9,0,by = -1)))
  yticklab <- format(c(seq(1000,1,by=-10),seq(9,0,by = -1)),digits = 1)
  axis(side = 2, at = axis.locs, labels = yticklab)
  
  axis.locs <- axTicks(side = 1)
  xticklab <- format(axis.locs*1000,digits = 1)
  axis(side = 1, at = axis.locs, labels = xticklab)
  
  abline(v=epoch_time[event_array[,event_idx[j,1]]==1][1]/1000,lwd=3,col="grey",lty=3)
  abline(v=median(epoch_time[event_array[,event_idx[j,2]]!=0][1])/1000,lwd=3,col="grey",lty=3)
  
  title1 <- paste("Cross-wavelet power between",region_name_1,"and",region_name_2,"\nin subject:",subject_index)
  title2 <- paste("Wavelet coherence between",region_name_1,"and",region_name_2,"\nin subject:",subject_index)
  
  title3 <- paste("Original signal of\n",region_name_1,"and",region_name_2,"in subject:",subject_index)
  title4 <- paste("Wavelet spectrum of\n",region_name_1,"in subject:",subject_index)
  title5 <- paste("Wavelet spectrum of\n",region_name_2,"in subject:",subject_index)
  
  par(p)
  
  p <- par(cex=0.5,oma=c(0,0,3,0))
  screen(7)
  matplot(epoch_time,t(xwtx$power.corr)[,idx],type="l",ylab="wave",xlab="Epoch time",main=title1)
  legend("topright",legend=paste(round(per2freq[idx,3],2),"Hz"),col=seq_along(idx),lty=seq_along(idx),bty="n")
  grid()
  abline(v=epoch_time[event_array[,event_idx[j,1]]==1][1],lwd=3,col="grey",lty=3)
  abline(v=median(epoch_time[event_array[,event_idx[j,2]]!=0]),lwd=3,col="grey",lty=3)
  
  
  screen(8)
  matplot(epoch_time,t(wtcx$rsq)[,idx],type="l",ylab="power",xlab="Epoch time",main=title2)
  legend("topright",legend=paste(round(per2freq[idx,3],2),"Hz"),col=seq_along(idx),lty=seq_along(idx),bty="n")
  grid()
  abline(v=epoch_time[event_array[,event_idx[j,1]]==1][1],lwd=3,col="grey",lty=3)
  abline(v=median(epoch_time[event_array[,event_idx[j,2]]!=0]),lwd=3,col="grey",lty=3)
  
  
  screen(3)
  
  xt_1s <- scale(xt_1,scale = FALSE)
  xt_2s <- scale(xt_2,scale = FALSE)
  matplot(epoch_time,xt_1s[,-1],type="l",ylab="vol",xlab="Epoch time",main=title3)
  matplot(epoch_time,xt_2s[,-1],type="l",ylab="vol",xlab="Epoch time",col = 3,lty = 1,main=title1,add = TRUE)
  legend("topright",legend=paste(round(per2freq[idx,3],2),"Hz"),col=seq_along(idx),lty=seq_along(idx),bty="n")
  grid()
  abline(v=epoch_time[event_array[,event_idx[j,1]]==1][1],lwd=3,col="grey",lty=3)
  abline(v=median(epoch_time[event_array[,event_idx[j,2]]!=0]),lwd=3,col="grey",lty=3)
  
  screen(6)
  matplot(epoch_time,t(wtx1$power.corr)[,idx],type="l",ylab="power",xlab="Epoch time",main=title4)
  legend("topright",legend=paste(round(per2freq[idx,3],2),"Hz"),col=seq_along(idx),lty=seq_along(idx),bty="n")
  grid()
  abline(v=epoch_time[event_array[,event_idx[j,1]]==1][1],lwd=3,col="grey",lty=3)
  abline(v=median(epoch_time[event_array[,event_idx[j,2]]!=0]),lwd=3,col="grey",lty=3)
  
  screen(9)
  matplot(epoch_time,t(wtx2$power.corr)[,idx],type="l",ylab="power",xlab="Epoch time",main=title5)
  legend("topright",legend=paste(round(per2freq[idx,3],2),"Hz"),col=seq_along(idx),lty=seq_along(idx),bty="n")
  grid()
  abline(v=epoch_time[event_array[,event_idx[j,1]]==1][1],lwd=3,col="grey",lty=3)
  abline(v=median(epoch_time[event_array[,event_idx[j,2]]!=0]),lwd=3,col="grey",lty=3)
  par(p)
  
  close.screen(all.screens = TRUE)
  
  
  
  if(!is.null(save)){
    dev.off()
  }
  if(return_val){
    return(list(wt = wtx, per2freq = per2freq))
  }
}


aggr_epoch_wt <- function(obj,
                          subject_index = 1,
                          event_type = 1,
                          region_index = 1,
                          interactive = TRUE,
                          parallel = FALSE){
  
  if(interactive){
    cat("subject_index:\n",seq_along(obj))
    subject_index <- as.numeric(readline("Select subject index:"))
  }
  
  i <- subject_index
  
  # select i-th object (i-th file/subject)
  obj_i <- obj[[i]]
  
  
  if(interactive){
    cat("event_type:\n",paste0(seq_along(obj_i$epoch),":",names(obj_i$epoch),"\n"))
    event_type <- as.numeric(readline("Select event_type:"))
    
    cat("region_index:\n",paste0(1:ncol(obj_i$signal),":",colnames(obj_i$signal),"\n"))
    region_index <- as.numeric(readline("Select region index:"))
  }
  
  j <- event_type
  k <- region_index
  
  
  #  l <- epoch_index
  
  region_name <- names(obj_i$signal)[k]
  
  # select j-th event (e.g. right/left)
  
  raw_array <- obj_i$epoch[[j]]
  epoch_time <- raw_array[,1,1]
  event_array <- raw_array[,grep(".event",dimnames(raw_array)[[2]]),]
  signal_array <- raw_array[,colnames(obj_i$signal),]
  
  event_type_name <- names(obj_i$epoch)[j]
  
  # select k-th region
  signal_array_k <- signal_array[,k,]
  signal_array_k_sub <- lapply(1:ncol(signal_array_k),function(ii)signal_array_k[,ii])
  
  if(parallel){
    cat("Wavelet transform over the all epochs ... \n")
    cl <- parallel::makeCluster(detectCores())
    
    cat("Parallel computing is applied. Number of available cores is",detectCores(),"\n")
    
    cat("Exporting objects into cores ... \n")
    
    parallel::clusterExport(cl=cl, varlist="epoch_time",envir = environment())
    cat("done\n")
    
    cat("Now parallel wavelet transforming ... ")
    
    
    tf_worker <- function(arr_sub){
      xt <- cbind(epoch_time/1000,arr_sub)
      wtx <- biwavelet::wt(xt)
      per2freq <- data.frame(index = paste("periodID",seq_along(wtx$period),sep="_"), period = wtx$period,freq = 1/wtx$period)
      wave <- t(Re(wtx$wave))
      power <- t(wtx$power)
      phase <- t(wtx$phase)
      
      colnames(wave) <- colnames(phase) <- colnames(power) <- per2freq[,1]
      list(per2freq=per2freq,wave = wave,power = power,phase = phase)
    }
    
    environment(tf_worker) <- .GlobalEnv
    wt_obj <- parallel::parLapply(cl, X = signal_array_k_sub,fun = tf_worker)
    
    parallel::stopCluster(cl)
    cat("done\n")
    
    
  }else{
    cat("Wavelet transform over the all epochs ... \n")
    wt_obj <- vector("list",ncol(signal_array_k))
    for(l in seq_along(wt_obj)){
      xt <- cbind(epoch_time/1000,signal_array_k[,l])
      wtx <- biwavelet::wt(xt)
      per2freq <- data.frame(index = paste("periodID",seq_along(wtx$period),sep="_"), period = wtx$period,freq = 1/wtx$period)
      wave <- t(Re(wtx$wave))
      power <- t(wtx$power)
      phase <- t(wtx$phase)
      
      colnames(wave) <- colnames(phase) <- colnames(power) <- per2freq[,1]
      wt_obj[[l]] <- list(per2freq=per2freq,wave = wave,power = power,phase = phase)
      cat(".")
    }
    cat("\n")
  }
  
  per2freq <- wt_obj[[1]]$per2freq
  
  cat("Creating array of wave, power, phase obtained from wavelet transform ... \n")
  wave_array <- power_array <- phase_array <- array(NA,dim = c(nrow(signal_array_k),
                                                               ncol(signal_array_k),
                                                               nrow(per2freq)),
                                                    dimnames = list(seq_along(epoch_time),colnames(signal_array_k),wt_obj[[1]]$per2freq[,1]))
  
  for(l in 1:ncol(signal_array_k)){
    wave_array[,l,] <- wt_obj[[l]]$wave
    power_array[,l,] <- wt_obj[[l]]$power
    phase_array[,l,] <- wt_obj[[l]]$phase
    cat(".")
  }
  cat("\n")
  cat("Aggregating array of wave, power, phase ... \n")
  
  aggr_origin_1 <- rowMeans(signal_array_k,1)
  aggr_wave <- apply(wave_array,c(1,3),median)
  aggr_power <- apply(power_array,c(1,3),median)
  
  
  cat("done! \n")
  
  list(subject = subject_index,
       region_index = region_index,
       region_name = region_name,
       event_type = event_type,
       event_type_name = event_type_name,
       origin_1 = signal_array_k,
       wave = wave_array,
       power = power_array,
       phase = phase_array,
       aggr_origin_1 = aggr_origin_1,
       aggr_wave = aggr_wave,
       aggr_power = aggr_power,
       epoch_time = epoch_time,
       per2freq = per2freq,
       Method = "wt")
}


aggr_epoch_xwtc <- function(obj,subject_index = 1,Method="wcoh",
                            event_type = 1,
                            region_index_1 = 1,
                            region_index_2 = 2,
                            interactive = TRUE,
                            parallel = FALSE){
  
  
  if(interactive){
    cat("subject_index:\n",seq_along(obj))
    subject_index <- as.numeric(readline("Select subject index:"))
  }
  
  
  i <- subject_index
  
  # select i-th object (i-th file/subject)
  obj_i <- obj[[i]]
  
  if(interactive){
    cat("event_type:\n",paste0(seq_along(obj_i$epoch),":",names(obj_i$epoch),"\n"))
    event_type <- as.numeric(readline("Select event_type:"))
    
    cat("region_index 1:\n",paste0(1:ncol(obj_i$signal),":",colnames(obj_i$signal),"\n"))
    region_index_1 <- as.numeric(readline("Select region index 1:"))
    
    cat("region_index 2:\n",paste0(1:ncol(obj_i$signal),":",colnames(obj_i$signal),"\n"))
    region_index_2 <- as.numeric(readline("Select region index 2:"))
    
  }
  
  j <- event_type
  k_1 <- region_index_1
  k_2 <- region_index_2
  
  region_name_1 <- names(obj_i$signal)[k_1]
  region_name_2 <- names(obj_i$signal)[k_2]
  
  # select j-th event (e.g. right/left)
  raw_array <- obj_i$epoch[[j]]
  epoch_time <- raw_array[,1,1]
  event_array <- raw_array[,grep(".event",dimnames(raw_array)[[2]]),]
  signal_array <- raw_array[,colnames(obj_i$signal),]
  
  event_type_name <- names(obj_i$epoch)[j]
  
  # select k-th region
  signal_array_k_1 <- signal_array[,k_1,]
  signal_array_k_2 <- signal_array[,k_2,]
  
  signal_array_k_sub <- lapply(1:ncol(signal_array_k_1),function(ii)cbind(signal_array_k_1[,ii],signal_array_k_2[,ii]))
  
  
  
  if(parallel){
    cat("Wavelet transform over the all epochs ... \n")
    cl <- parallel::makeCluster(detectCores())
    
    cat("Parallel computing is applied. Number of available cores is",detectCores(),"\n")
    
    cat("Objects are exported to each core. ... \n")
    parallel::clusterExport(cl=cl, varlist=c("signal_array_k_sub", "epoch_time","Method"),envir = environment())
    cat("done\n")
    
    if(Method == "xwt"){
      cat("Now parallel calculation of cross-wavelet over all epochs... \n")
    }else if(Method =="wcoh"){
      cat("Now parallel calculation of wavelet coherence over all epochs... \n")
      message("Since we are calculating with nrands = 0, \nif you want to estimate the signif value correctly, recalculate with larger nrands!")
    }else{
      stop("There is no corresponding method, please select either 'xwt' or 'wcoh' as the method argument.")
    }
    
    
    tf_worker <- function(arr_sub){
      
      
      xt_1 <- cbind(epoch_time/1000,arr_sub[,1])
      xt_2 <- cbind(epoch_time/1000,arr_sub[,2])
      
      if(Method == "xwt"){
        wtx <- biwavelet::xwt(xt_1,xt_2)
      }else if(Method == "wcoh"){
        wtx <- biwavelet::wtc(xt_1,xt_2,nrands = 0,quiet = TRUE)
      }
      
      per2freq <- data.frame(index = paste("periodID",seq_along(wtx$period),sep="_"), period = wtx$period,freq = 1/wtx$period)
      
      wave <- t(Re(wtx$wave))
      power <- t(wtx$power)
      phase <- t(wtx$phase)
      
      if(Method == "wcoh"){
        rsq <- t(wtx$rsq)
        colnames(rsq) <-  per2freq[,1]
      }
      
      colnames(wave) <- colnames(wave) <- colnames(phase) <- per2freq[,1]
      if(Method == "xwt"){
        
        list(per2freq=per2freq,
             wave = wave,
             power = power,
             phase = phase)
        
      }else if(Method == "wcoh"){
        
        list(per2freq=per2freq,
             wcoh = rsq,
             wave = wave,
             power = power,
             phase = phase)
        
      }
    }
    
    environment(tf_worker) <- .GlobalEnv
    wt_obj <- parallel::parLapply(cl, X = signal_array_k_sub,fun = tf_worker)
    
    
    parallel::stopCluster(cl)
    cat("done\n")
    
    
  }else{
    
    cat("Now calculating cross-wavelet and wavelet coherence over all epochs ... \n")
    message("Since we are calculating with nrands = 0,\n if you want to estimate the signif value correctly, recalculate with larger nrands!")
    
    
    wt_obj <- vector("list",ncol(signal_array_k_1))
    for(l in seq_along(wt_obj)){
      
      xt_1 <- cbind(epoch_time/1000,signal_array_k_1[,l])
      xt_2 <- cbind(epoch_time/1000,signal_array_k_2[,l])
      
      if(Method == "xwt"){
        wtx <- biwavelet::xwt(xt_1,xt_2)
      }else if (Method == "wcoh"){
        wtx <- biwavelet::wtc(xt_1,xt_2,nrands = 0,quiet = TRUE)
      }
      
      per2freq <- data.frame(index = paste("periodID",seq_along(wtx$period),sep="_"), period = wtx$period,freq = 1/wtx$period)
      
      wave <- t(Re(wtx$wave))
      power <- t(wtx$power)
      phase <- t(wtx$phase)
      
      if(Method == "wcoh"){
        rsq <- t(wtx$rsq)
        colnames(rsq) <-  per2freq[,1]
      }
      
      colnames(wave) <- colnames(wave) <- colnames(phase) <- per2freq[,1]
      
      if(Method == "xwt"){
        wt_obj[[l]] <- list(per2freq=per2freq,
                            wave = wave,
                            power = power,
                            phase = phase)
      }else if(Method == "wcoh"){
        wt_obj[[l]] <- list(per2freq=per2freq,
                            wcoh = rsq,
                            wave = wave,
                            power = power,
                            phase = phase)
      }
      
      cat(".")
    }
    cat("\n")
  }
  
  
  
  wave_array <- power_array <- phase_array <- 
    array(NA,dim = c(nrow(signal_array_k_1),ncol(signal_array_k_1),nrow(per2freq)),
          dimnames = list(seq_along(epoch_time),colnames(signal_array_k_1),wt_obj[[1]]$per2freq[,1]))
  
  per2freq <- wt_obj[[1]]$per2freq
  
  if(Method == "wcoh"){
    cat("Creating array of wave, power, phase, coherence obtained from wavelet transform ... \n")
    rsq_array <- array(NA,dim = c(nrow(signal_array_k_1),ncol(signal_array_k_1),nrow(per2freq)),
                       dimnames = list(seq_along(epoch_time),colnames(signal_array_k_1),per2freq[,1]))
    for(l in 1:ncol(signal_array_k_1)){
      rsq_array[,l,] <- wt_obj[[l]]$wcoh
      wave_array[,l,] <- wt_obj[[l]]$wave
      power_array[,l,] <- wt_obj[[l]]$power
      phase_array[,l,] <- wt_obj[[l]]$phase
      cat(".")
    }
    
    cat("\nAggregating array of wave, power, coherence ... \n")
    
    aggr_rsq <- apply(rsq_array,c(1,3),median)
    aggr_wave <- apply(wave_array,c(1,3),median)
    aggr_power <- apply(power_array,c(1,3),median)
    aggr_origin_1 <- rowMeans(signal_array_k_1)
    aggr_origin_2 <- rowMeans(signal_array_k_2)
    
  }else if(Method == "xwt"){
    cat("Creating array of wave, power, phase obtained from wavelet transform ... \n")
    
    for(l in 1:ncol(signal_array_k_1)){
      wave_array[,l,] <- wt_obj[[l]]$wave
      power_array[,l,] <- wt_obj[[l]]$power
      phase_array[,l,] <- wt_obj[[l]]$phase
      cat(".")
    }
    cat("Aggregating array of wave, power ... \n")
    
    aggr_wave <- apply(wave_array,c(1,3),median)
    aggr_power <- apply(power_array,c(1,3),median)
    aggr_origin_1 <- rowMeans(signal_array_k_1)
    aggr_origin_2 <- rowMeans(signal_array_k_2)
    
  }
  cat("\n")
  cat("done! \n")
  
  
  if(Method == "xwt"){
    
    list(subject = subject_index,
         region_index = c(region_index_1,region_index_2),
         region_name_1 = region_name_1,
         region_name_2 = region_name_2, 
         event_type = event_type,
         event_type_name = event_type_name,
         origin_1 = signal_array_k_1,
         origin_2 = signal_array_k_2,
         wave = wave_array,
         power = power_array,
         phase = phase_array,
         aggr_origin_1 = aggr_origin_1,
         aggr_origin_2 = aggr_origin_2,
         aggr_wave = aggr_wave,
         aggr_power = aggr_power,
         epoch_time = epoch_time,
         per2freq = per2freq,
         Method = Method)
    
  }else if(Method == "wcoh"){
    
    list(subject = subject_index,
         region_index = c(region_index_1,region_index_2),
         region_name_1 = region_name_1,
         region_name_2 = region_name_2, 
         event_type = event_type,
         event_type_name = event_type_name,
         origin_1 = signal_array_k_1,
         origin_2 = signal_array_k_2,
         rsq = rsq_array,
         wave = wave_array,
         power = power_array,
         phase = phase_array,
         aggr_origin_1 = aggr_origin_1,
         aggr_origin_2 = aggr_origin_2,
         aggr_rsq = aggr_rsq,
         aggr_wave = aggr_wave,
         aggr_power = aggr_power,
         epoch_time = epoch_time,
         per2freq = per2freq,
         Method = Method)
  }
}

aggr_epoch <- function(obj,subject_index = 1,Method="wcoh",
                       event_type = 1,
                       region_index_1 = 1,
                       region_index_2 = 2,
                       interactive = TRUE,
                       parallel = FALSE){
  
  
  if(Method == "wt"){
    aggr_epoch_wt(obj,subject_index = subject_index,
                  event_type = event_type,
                  region_index = region_index_1,
                  interactive = interactive,
                  parallel = parallel)
    
  }else if(Method %in% c("wcoh","xwt")){
    
    aggr_epoch_xwtc(obj,subject_index = subject_index,Method = Method,
                    event_type = event_type,
                    region_index_1 = region_index_1,
                    region_index_2 = region_index_2,
                    interactive = interactive,
                    parallel = parallel)
    
  }
  
}


simple_wavelet <- function(obj,subject_index = 1,Method="wcoh",
                           event_type = 1,
                           region_index_1 = 1,
                           region_index_2 = 2,
                           interactive = TRUE){
  
  
  if(Method == "wt"){
    simple_wt(obj,subject_index = subject_index,
              event_type = event_type,
              region_index_1 = region_index_1,
              interactive = interactive)
    
  }else if(Method %in% c("wcoh","xwt")){
    
    simple_xwtc(obj,subject_index = subject_index,Method = Method,
                event_type = event_type,
                region_index_1 = region_index_1,
                region_index_2 = region_index_2,
                interactive = interactive)
    
  }
  
}

simple_wt <- function(obj,subject_index = 1,
                      event_type = 1,
                      region_index_1 = 1,
                      interactive = TRUE){
  
  if(interactive){
    cat("subject_index:\n",seq_along(obj))
    subject_index <- as.numeric(readline("Select subject index:"))
  }
  
  i <- subject_index
  
  # select i-th object (i-th file/subject)
  obj_i <- obj[[i]]
  
  
  if(interactive){
    cat("event_type:\n",paste0(seq_along(obj_i$epoch),":",names(obj_i$epoch),"\n"))
    event_type <- as.numeric(readline("Select event_type:"))
    
    cat("region_index:\n",paste0(1:ncol(obj_i$signal),":",colnames(obj_i$signal),"\n"))
    region_index <- as.numeric(readline("Select region index:"))
  }
  
  j <- event_type
  k <- region_index
  
  
  #  l <- epoch_index
  
  region_name <- names(obj_i$signal)[k]
  
  # select j-th event (e.g. right/left)
  
  raw_array <- obj_i$epoch[[j]]
  epoch_time <- raw_array[,1]
  event_array <- raw_array[,grep(".event",dimnames(raw_array)[[2]])]
  signal_array <- raw_array[,colnames(obj_i$signal)]
  
  event_type_name <- names(obj_i$epoch)[j]
  
  # select k-th region
  signal_array_k <- signal_array[,k]
  
  xt <- cbind(epoch_time/1000,signal_array_k)
  wtx <- biwavelet::wt(xt)
  per2freq <- data.frame(index = paste("periodID",seq_along(wtx$period),sep="_"), period = wtx$period,freq = 1/wtx$period)
  wave <- t(Re(wtx$wave))
  power <- t(wtx$power)
  phase <- t(wtx$phase)
  
  colnames(wave) <- colnames(phase) <- colnames(power) <- per2freq[,1]
  wt_obj <- list(per2freq=per2freq,wave = wave,power = power,phase = phase)
  
  list(subject = subject_index,
       region_index = region_index,
       region_name = region_name,
       event_type = event_type,
       event_type_name = event_type_name,
       origin_1 = signal_array_k,
       wave = wave,
       power = power,
       phase = phase,
       epoch_time = epoch_time,
       per2freq = per2freq,
       Method = "wt")
}

simple_xwtc <- function(obj,subject_index = 1,Method="wcoh",
                        event_type = 1,
                        region_index_1 = 1,
                        region_index_2 = 2,
                        interactive = TRUE){
  
  if(interactive){
    cat("subject_index:\n",seq_along(obj))
    subject_index <- as.numeric(readline("Select subject index:"))
  }
  
  i <- subject_index
  
  # select i-th object (i-th file/subject)
  obj_i <- obj[[i]]
  
  if(interactive){
    cat("event_type:\n",paste0(seq_along(obj_i$epoch),":",names(obj_i$epoch),"\n"))
    event_type <- as.numeric(readline("Select event_type:"))
    
    cat("region_index 1:\n",paste0(1:ncol(obj_i$signal),":",colnames(obj_i$signal),"\n"))
    region_index_1 <- as.numeric(readline("Select region index 1:"))
    
    cat("region_index 2:\n",paste0(1:ncol(obj_i$signal),":",colnames(obj_i$signal),"\n"))
    region_index_2 <- as.numeric(readline("Select region index 2:"))
    
  }
  
  j <- event_type
  k_1 <- region_index_1
  k_2 <- region_index_2
  
  region_name_1 <- names(obj_i$signal)[k_1]
  region_name_2 <- names(obj_i$signal)[k_2]
  
  # select j-th event (e.g. right/left)
  raw_array <- obj_i$epoch[[j]]
  epoch_time <- raw_array[,1]
  event_array <- raw_array[,grep(".event",dimnames(raw_array)[[2]])]
  signal_array <- raw_array[,colnames(obj_i$signal)]
  
  event_type_name <- names(obj_i$epoch)[j]
  
  # select k-th region
  signal_array_k_1 <- signal_array[,k_1]
  signal_array_k_2 <- signal_array[,k_2]
  
  signal_array_k_sub <- cbind(signal_array_k_1,signal_array_k_2)
  
  # wavelet analysis
  
  xt_1 <- cbind(epoch_time/1000,signal_array_k_1)
  xt_2 <- cbind(epoch_time/1000,signal_array_k_2)
  
  if(Method == "xwt"){
    wtx <- biwavelet::xwt(xt_1,xt_2)
  }else if (Method == "wcoh"){
    wtx <- biwavelet::wtc(xt_1,xt_2,nrands = 0,quiet = TRUE)
  }
  
  per2freq <- data.frame(index = paste("periodID",seq_along(wtx$period),sep="_"), period = wtx$period,freq = 1/wtx$period)
  
  wave <- t(Re(wtx$wave))
  power <- t(wtx$power)
  phase <- t(wtx$phase)
  
  if(Method == "wcoh"){
    rsq <- t(wtx$rsq)
    colnames(rsq) <-  per2freq[,1]
  }
  
  colnames(wave) <- colnames(wave) <- colnames(phase) <- per2freq[,1]
  
  if(Method == "xwt"){
    
    wt_obj <- list(per2freq=per2freq,
                   wave = wave,
                   power = power,
                   phase = phase)
    
  }else if(Method == "wcoh"){
    
    wt_obj <- list(per2freq=per2freq,
                   wcoh = rsq,
                   wave = wave,
                   power = power,
                   phase = phase)
    
  }
  
  
  
  if(Method == "xwt"){
    
    list(subject = subject_index,
         region_index = c(region_index_1,region_index_2),
         region_name_1 = region_name_1,
         region_name_2 = region_name_2, 
         event_type = event_type,
         event_type_name = event_type_name,
         origin_1 = signal_array_k_1,
         origin_2 = signal_array_k_2,
         wave = wave,
         power = power,
         phase = phase,
         epoch_time = epoch_time,
         per2freq = per2freq,
         Method = Method)
    
  }else if(Method == "wcoh"){
    
    list(subject = subject_index,
         region_index = c(region_index_1,region_index_2),
         region_name_1 = region_name_1,
         region_name_2 = region_name_2, 
         event_type = event_type,
         event_type_name = event_type_name,
         origin_1 = signal_array_k_1,
         origin_2 = signal_array_k_2,
         rsq = rsq,
         wave = wave,
         power = power,
         phase = phase,
         epoch_time = epoch_time,
         per2freq = per2freq,
         Method = Method)
  }
  
}

plot_wavelet <- function(wt_obj,freq = 20){
  
  wave <- wt_obj$wave
  power <- wt_obj$power
  aggr_wave <- wt_obj$aggr_wave
  aggr_power <- wt_obj$aggr_power
  
  if(wt_obj$Method == "wcoh"){
    rsq <- wt_obj$rsq 
    aggr_rsq <- wt_obj$aggr_rsq
  }
  
  epoch_time <- wt_obj$epoch_time
  
  subject <- wt_obj$subject
  region_name <- wt_obj$region_name
  event_type_name <- wt_obj$event_type_name
  
  cat("Approximating specified frequency to obtain wavelet frequency ... \n")
  per2freq <- wt_obj$per2freq
  idx <- which.min(abs(per2freq[,3] - freq))
  
  cat("The closest wavelet frequency was",round(per2freq[idx,3],2),"Hz (in terms of period:",round(per2freq[idx,2],2),")\n")
  
  cat("Creating graphics ... ")
  p <- par(cex=0.5)
  if(wt_obj$Method == "wcoh"){
    
    origin_1 <- wt_obj$origin_1
    origin_2 <- wt_obj$origin_2
    aggr_origin_1 <- wt_obj$aggr_origin_1
    aggr_origin_2 <- wt_obj$aggr_origin_2
    
    split.screen(c(3,2))
    title1 <- paste("Aggregated singnal 1 over epochs\n","subject:",subject,"region:",region_name,"event type",event_type_name)
    title2 <- paste("Aggregated singnal 2 over epochs\n","subject:",subject,"region:",region_name,"event type",event_type_name)
    title3 <- paste("Powerband by",freq,"Hz for each epoch\n","subject:",subject,"region:",region_name,"event type",event_type_name)
    title4 <- paste("Aggregated powerband over epochs by",freq,"Hz for each epoch\n","subject:",subject,"region:",region_name,"event type",event_type_name)
    title5 <- paste("Wavelet coherence by",freq,"Hz for each epoch\n","subject:",subject,"region:",region_name,"event type",event_type_name)
    title6 <- paste("Aggregated wavelet coherence by",freq,"Hz for each epoch\n","subject:",subject,"region:",region_name,"event type",event_type_name)
    
    screen(1)
    matplot(epoch_time,aggr_origin_1,type="l",main = title1,xlab="epoch time")
    
    screen(2)
    matplot(epoch_time,aggr_origin_2,type="l",main = title2,xlab="epoch time")
    
    screen(3)
    matplot(epoch_time,power[,,idx],type="l",main = title3,xlab="epoch time")
    
    screen(4)
    aggr_power_ci <- t(apply(power[,,idx],1,function(xx)quantile(xx,probs = c(.05,.95))))
    matplot(epoch_time,aggr_power[,idx],ylim = range(aggr_power_ci),type="l",main = title4,xlab="epoch time")
    matplot(epoch_time,aggr_power_ci,type="l",col = 4,lty=2,xlab="epoch time",add = TRUE)
    
    screen(5)
    matplot(epoch_time,rsq[,,idx],type="l",main = title5,xlab="epoch time",ylim=c(0,1))
    
    screen(6)
    aggr_rsq_ci <- t(apply(rsq[,,idx],1,function(xx)quantile(xx,probs = c(.05,.95))))
    matplot(epoch_time,aggr_rsq[,idx],ylim=c(0,1),type="l",main = title6,xlab="epoch time")
    matplot(epoch_time,aggr_rsq_ci,type="l",col = 4,lty=2,xlab="epoch time",add = TRUE)
    
    close.screen(all.screens = TRUE)
    
    
  }else if(wt_obj$Method == "xwt"){
    aggr_origin_1 <- wt_obj$aggr_origin_1
    aggr_origin_2 <- wt_obj$aggr_origin_2
    aggr_power <- wt_obj$aggr_power
    
    split.screen(c(2,2))
    title1 <- paste("Aggregated original singnal 1 over epochs\n","subject:",subject,"region:",region_name,"event type",event_type_name)
    title2 <- paste("Aggregated original singnal 2 over epochs\n","subject:",subject,"region:",region_name,"event type",event_type_name)
    title3 <- paste("Cross-wavelet powerband by",freq,"Hz for each epoch\n","subject:",subject,"region:",region_name,"event type",event_type_name)
    title4 <- paste("Aggregated corss-wavelet powerband over epochs by",freq,"Hz for each epoch\n","subject:",subject,"region:",region_name,"event type",event_type_name)
    
    screen(1)
    matplot(epoch_time,aggr_origin_1,type="l",main = title1,xlab="epoch time")
    
    screen(2)
    matplot(epoch_time,aggr_origin_2,type="l",main = title2,xlab="epoch time")
    
    screen(3)
    matplot(epoch_time,power[,,idx],type="l",main = title3,xlab="epoch time")
    
    screen(4)
    aggr_power_ci <- t(apply(power[,,idx],1,function(xx)quantile(xx,probs = c(.05,.95))))
    matplot(epoch_time,aggr_power[,idx],type="l",main = title4,xlab="epoch time")
    matplot(epoch_time,aggr_power_ci,type="l",col = 4,lty=2,xlab="epoch time",add = TRUE)
    
    
    close.screen(all.screens = TRUE)
    
  }else if(wt_obj$Method == "wt"){
    aggr_origin_1 <- wt_obj$aggr_origin_1
    aggr_power <- wt_obj$aggr_power
    
    split.screen(c(3,1))
    title1 <- paste("Aggregated original singnal over epochs\n","subject:",subject,"region:",region_name,"event type",event_type_name)
    title2 <- paste("Wavelet powerband by",freq,"Hz for each epoch\n","subject:",subject,"region:",region_name,"event type",event_type_name)
    title3 <- paste("Aggregated wavelet powerband over epochs by",freq,"Hz for each epoch\n","subject:",subject,"region:",region_name,"event type",event_type_name)
    
    screen(1)
    matplot(epoch_time,aggr_origin_1,type="l",main = title1,xlab="epoch time")
    
    screen(2)
    matplot(epoch_time,power[,,idx],type="l",main = title2,xlab="epoch time")
    
    screen(3)
    matplot(epoch_time,aggr_power[,idx],type="l",main = title3,xlab="epoch time")
    
    close.screen(all.screens = TRUE)
    
  }
  cat("done\n")
  par(p)
}

plot_wavelet_simple <- function(wt_obj,freq = 20){
  
  wave <- wt_obj$wave
  power <- wt_obj$power
  phase <- wt_obj$phase
  
  if(wt_obj$Method == "wcoh"){
    rsq <- wt_obj$rsq 
    aggr_rsq <- wt_obj$aggr_rsq
  }
  
  epoch_time <- wt_obj$epoch_time
  
  subject <- wt_obj$subject
  region_name <- wt_obj$region_name
  event_type_name <- wt_obj$event_type_name
  
  cat("Approximating specified frequency to obtain wavelet frequency ... \n")
  per2freq <- wt_obj$per2freq
  idx <- which.min(abs(per2freq[,3] - freq))
  
  cat("The closest wavelet frequency was",round(per2freq[idx,3],2),"Hz (in terms of period:",round(per2freq[idx,2],2),")\n")
  
  cat("Creating graphics ... ")
  p <- par(cex = 0.5)
  
  if(wt_obj$Method == "wcoh"){
    
    origin_1 <- wt_obj$origin_1
    origin_2 <- wt_obj$origin_2
    region_name_1 <- wt_obj$region_name_1
    region_name_2 <- wt_obj$region_name_2
    
    split.screen(c(2,2))
    title1 <- paste("Singnal 1\n","subject:",subject,"region:",region_name_1,"event type",event_type_name)
    title2 <- paste("Singnal 2\n","subject:",subject,"region:",region_name_2,"event type",event_type_name)
    title3 <- paste("Phase by",freq,"Hz for each epoch\n","subject:",subject,"region:",region_name,"event type",event_type_name)
    title4 <- paste("Wavelet coherence by",freq,"Hz for each epoch\n","subject:",subject,"region:",region_name_1,"and",region_name_2,"event type",event_type_name)
    
    screen(1)
    matplot(epoch_time,origin_1,type="l",main = title1,xlab="epoch time")
    grid()
    
    screen(2)
    matplot(epoch_time,phase[,idx],type="l",main = title3,xlab="epoch time")
    grid()
    
    screen(3)
    matplot(epoch_time,origin_2,type="l",main = title2,xlab="epoch time")
    grid()
    
    
    screen(4)
    matplot(epoch_time,rsq[,idx],type="l",main = title4,xlab="epoch time",ylim=c(0,1))
    grid()
    
    
    close.screen(all.screens = TRUE)
    
    
  }else if(wt_obj$Method == "xwt"){
    
    origin_1 <- wt_obj$origin_1
    origin_2 <- wt_obj$origin_2
    region_name_1 <- wt_obj$region_name_1
    region_name_2 <- wt_obj$region_name_2
    
    split.screen(c(2,2))
    title1 <- paste("Singnal 1\n","subject:",subject,"region:",region_name_1,"event type",event_type_name)
    title2 <- paste("Singnal 2\n","subject:",subject,"region:",region_name_2,"event type",event_type_name)
    title3 <- paste("Phase by",freq,"Hz for each epoch\n","subject:",subject,"region:",region_name,"event type",event_type_name)
    title4 <- paste("Cross-wavelet spectrum by",freq,"Hz for each epoch\n","subject:",subject,"region:",region_name_1,"and",region_name_2,"event type",event_type_name)
    
    screen(1)
    matplot(epoch_time,origin_1,type="l",main = title1,xlab="epoch time")
    grid()
    
    screen(2)
    matplot(epoch_time,phase[,idx],type="l",main = title3,xlab="epoch time")
    grid()
    
    screen(3)
    matplot(epoch_time,origin_2,type="l",main = title2,xlab="epoch time")
    grid()
    
    screen(4)
    matplot(epoch_time,power[,idx],type="l",main = title4,xlab="epoch time")
    grid()
    
    close.screen(all.screens = TRUE)
    
  }else if(wt_obj$Method == "wt"){
    origin_1 <- wt_obj$origin_1
    
    split.screen(c(2,1))
    title1 <- paste("Singnal over epochs\n","subject:",subject,"region:",region_name,"event type",event_type_name)
    title2 <- paste("Wavelet powerband by",freq,"Hz for each epoch\n","subject:",subject,"region:",region_name,"event type",event_type_name)
    
    screen(1)
    matplot(epoch_time,origin_1,type="l",main = title1,xlab="epoch time")
    grid()
    
    screen(2)
    matplot(epoch_time,power[,idx],type="l",main = title2,xlab="epoch time")
    grid()
    
    
    close.screen(all.screens = TRUE)
    
  }
  cat("done\n")
  par(p)
}

select_epoch <- function(obj,epoch_key = 1:10,resume = FALSE){
  
  for(i in seq_along(obj)){
    obj_i <- obj[[i]]
    for(j in seq_along(obj_i$epoch)){
      epoch <- obj_i$epoch[[j]]
      obj_i$epoch[[j]] <- epoch[,,epoch_key]
    }
    obj[[i]] <- obj_i
  }
  obj
}

simplify <- function(obj,aggr_method = "median"){
  names(obj) <- paste("Subject",seq_along(obj),sep="_")
  for(i in seq_along(obj)){
    obj_i <- obj[[i]]
    for(j in seq_along(obj_i$epoch)){
      epoch <- obj_i$epoch[[j]]
      epoch_aggr <- apply(epoch,c(1,2),aggr_method)
      obj_i$epoch[[j]] <- epoch_aggr
    }
    obj[[i]] <- obj_i
  }
  #  obj$simplify <- TRUE
  obj
}

xwtc_simple <- function(obj,
                        subject_index = 1,
                        region_indices = c(1,2,3,4,5),
                        freq = c(5,10,15,20),
                        event_type = 1,
                        interactive = TRUE,
                        Method = "wcoh"){
  
  if(interactive){
    cat("subject_index:\n",seq_along(obj))
    subject_index <- as.numeric(readline("Select subject index:"))
  }
  
  i <- subject_index
  
  # select i-th object (i-th file/subject)
  obj_i <- obj[[i]]
  
  if(interactive){
    cat("event_type:\n",paste0(seq_along(obj_i$epoch),":",names(obj_i$epoch),"\n"))
    event_type <- as.numeric(readline("Select event_type:"))
  }
  
  j <- event_type
  
  
  # select j-th event (e.g. right/left)
  raw_array <- obj_i$epoch[[j]]
  epoch_time <- raw_array[,1]
  event_array <- raw_array[,grep(".event",dimnames(raw_array)[[2]])]
  signal_array <- raw_array[,colnames(obj_i$signal)]
  event_type_name <- names(obj_i$epoch)[j]
  
  tmp_signal_array_k_1 <- signal_array[,1]
  tmp_xt_1 <- cbind(epoch_time / 1000, tmp_signal_array_k_1)
  tmp_wtx <- wt(tmp_xt_1)
  
  
  cat("Approximating specified frequency to obtain wavelet frequency ... \n")
  per2freq <- data.frame(index = paste("periodID",seq_along(tmp_wtx$period),sep="_"), period = tmp_wtx$period,freq = 1/tmp_wtx$period)
  idx <- sapply(freq,function(xx)which.min((xx - per2freq[,3])^2))
  cat("The closest wavelet frequency was: ",paste(round(per2freq[idx,3],2),"Hz (in terms of period:",round(per2freq[idx,2],2),")"),sep="\n")
  
  
  if(Method == "xwxt"){
    cat("Creating cross-spectrum array ... \n")
  }else if(Method == "wcoh"){
    cat("Creating coherence array ... \n")
  }else if(Method == "wt"){
    cat("Creating wavelet array ... \n")
  }
  
  if(Method == "wt"){
    
    xwtc_array <- inv_xwtc_array <- array(NA,
                                          dim = c(length(region_indices),
                                                  nrow(obj_i$epoch[[j]]),
                                                  length(idx)),
                                          dimnames = list(region_name_1 = names(obj_i$signal)[region_indices],
                                                          epoch_time = epoch_time,
                                                          frequency = freq))
    
    xwtc_li <- vector("list",length(region_indices))
    
    for(k in seq_along(region_indices)){
      
      k_1 <- region_indices[k]
      region_name_1 <- names(obj_i$signal)[k_1]
      signal_array_k_1 <- signal_array[,k_1]
      xt_1 <- cbind(epoch_time/1000,signal_array_k_1)
      
      wtx <- biwavelet::wt(xt_1)
      power <- t(wtx$power)
      wave <- t(Re(wtx$wave))
      colnames(power) <- per2freq[,1]
      colnames(wave) <- per2freq[,1]
      xwtc_array[k,,] <- power[,idx]
      inv_xwtc_array[k,,] <- wave[,idx]
      
      xwtc_li[[k]]$wtx <- wtx
      xwtc_li[[k]]$k_1 <- k_1
      cat(".")
    }
    cat("done! \n")
    region_names <- names(obj_i$signal)[region_indices]
    
    return(list(xwtc_array = xwtc_array,
                inv_xwtc_array = inv_xwtc_array,
                xwtc = xwtc_li,
                subject_index = subject_index,
                region = region_indices,
                region_name = region_names,
                freq = freq,
                event_array = event_array,
                event_type = event_type,
                event_type_name = event_type_name,
                epoch_time = epoch_time,
                Method = Method))
    
  }else{
    xwtc_array <- array(NA,
                        dim = c(length(region_indices),
                                length(region_indices),
                                nrow(obj_i$epoch[[j]]),
                                length(idx)),
                        dimnames = list(region_name_1 = names(obj_i$signal)[region_indices],
                                        region_name_2 =names(obj_i$signal)[region_indices],
                                        epoch_time = epoch_time,
                                        frequency = freq))
    
    region_xindex <- combn(region_indices,2)
    xwtc_li <- vector("list",ncol(region_xindex))
    
    for(k in 1:ncol(region_xindex)){
      k_1 <- region_xindex[1,k]
      k_2 <- region_xindex[2,k]
      
      region_name_1 <- names(obj_i$signal)[k_1]
      region_name_2 <- names(obj_i$signal)[k_2]
      
      
      # select k-th region
      signal_array_k_1 <- signal_array[,k_1]
      signal_array_k_2 <- signal_array[,k_2]
      
      #signal_array_k_sub <- cbind(signal_array_k_1,signal_array_k_2)
      
      # wavelet analysis
      
      xt_1 <- cbind(epoch_time/1000,signal_array_k_1)
      xt_2 <- cbind(epoch_time/1000,signal_array_k_2)
      
      if(Method == "xwt"){
        wtx <- biwavelet::xwt(xt_1,xt_2)
        power <- t(wtx$power)
        colnames(power) <- per2freq[,1]
        xwtc_array[k_1,k_2,] <- power[,idx]
        xwtc_array[k_2,k_1,] <- power[,idx]
        
      }else if(Method == "wcoh"){
        wtx <- biwavelet::wtc(xt_1,xt_2,nrands = 0,quiet = TRUE)
        rsq <- t(wtx$rsq)
        colnames(rsq) <-  per2freq[,1]
        xwtc_array[k_1,k_2,,] <- rsq[,idx]
        xwtc_array[k_2,k_1,,] <- rsq[,idx]
      }
      
      xwtc_li[[k]]$wtx <- wtx
      xwtc_li[[k]]$k_1 <- k_1
      xwtc_li[[k]]$k_2 <- k_2
      
      cat(".")
    }
    if(Method == "xwt"){
      xwtc_array[is.na(xwtc_array)] <- Inf
    }else if(Method == "wcoh"){
      xwtc_array[is.na(xwtc_array)] <- 1
    }
    
    cat("done! \n")
    
    xregion_name <-  rbind(names(obj_i$signal)[region_xindex[1,]],names(obj_i$signal)[region_xindex[2,]])
    
    return(list(xwtc_array = xwtc_array,
                xwtc = xwtc_li,
                subject_index = subject_index,
                xregion = region_xindex,
                xregion_name = xregion_name,
                freq = freq,
                event_array = event_array,
                event_type = event_type,
                event_type_name = event_type_name,
                epoch_time = epoch_time,
                Method = Method))
  }
  
}

xwtc_simple_all <- function(obj,
                            region_indices = c(1,2,3,4,5),
                            freq = c(5,10,15,20),
                            event_type = 1,
                            interactive = TRUE,
                            Method = "wcoh"){
  
  results <- vector("list",length(obj))
  names(results) <- paste("subject",seq_along(results),sep="_")
  
  for(i in seq_along(results)){
    results[[i]] <- xwtc_simple(obj,
                                subject_index = i,
                                region_indices = region_indices,
                                freq = freq,
                                event_type = event_type,
                                interactive = FALSE,
                                Method = Method)
    cat("\n")
    cat("subjct",i,"finished!\n\n")
  }
  
  return(results)
  
}


plot_wt_simple <- function(wt_obj,log=FALSE,norm = TRUE,file_path = "plot.png"){
  
  power_array <- wt_obj$xwtc_array
  region_names <- wt_obj$region_name
  freq <- wt_obj$freq
  epoch_time <- wt_obj$epoch_time
  event_array <- wt_obj$event_array
  
  
  file_path_2 <- paste0(paste(gsub(".png","",basename(file_path)),sep="_"),".png")
  png(file_path_2,width = 2000,height = 2000,type="cairo")
  
  split.screen(c(4,4))
  for(f in seq_along(freq)){
    screen(f)
    p <- par(cex = 0.8,mar = c(3,3,3,1))
    xx <- t(power_array[,,f])
    
    if(log){
      xx <- log10(xx)
    }
    
    if(norm){
      #xx <- apply(xx,2,function(x)x/max(x))
      xx <- scale(xx)
    }
    
    event_idx <- matrix(c(1,2,3,4),byrow = TRUE,ncol = 2)
    
    matplot(epoch_time,xx,xlab = "epoch time",ylab="wavelet power",type = "l",lty = 1:ncol(xx),col = 1:ncol(xx))
    legend("topleft",
           legend = region_names,
           lty = seq_along(region_names),
           col = seq_along(region_names),
           bty = "n")
    title(paste(freq[f],"Hz","wavelet power among muscles"))
    
    abline(v=epoch_time[event_array[,event_idx[j,1]]==1][1],lwd=3,col="grey",lty=3)
    abline(v=median(epoch_time[event_array[,event_idx[j,2]]!=0]),lwd=3,col="grey",lty=3)
    
    
    grid()
    par(p)
  }
  close.screen(all.screens = TRUE)
  dev.off()
}

plot_xwtc_simple <- function(xwtc_obj,plot_type = "sg",file_path = "plot.png"){
  
  xregion <- xwtc_obj$xregion
  xregion_name <- xwtc_obj$xregion_name
  event_array <- xwtc_obj$event_array
  event_type <- xwtc_obj$event_type
  Method <- xwtc_obj$Method
  subject_index <- xwtc_obj$subject
  freq <- xwtc_obj$freq
  
  event_idx <- matrix(c(1,2,3,4),byrow = TRUE,ncol = 2)
  
  if(plot_type == "sg"){
    
    xwtc_li <- xwtc_obj$xwtc
    xwtc_array <- xwtc_obj$xwtc_array
    epoch_time <- as.numeric(dimnames(xwtc_array)[[3]])
    
    
    p <- par(cex = 0.7)
    
    file_path_2 <- paste0(paste(gsub(".png","",basename(file_path)),plot_type,sep="_"),".png")
    png(file_path_2,width = 2000,height = 2000,type="cairo")
    
    idx <- split.screen(dim(xwtc_array)[1:2])
    layout_matrix <- matrix(idx,dim(xwtc_array)[1],byrow = TRUE)
    upper_idx <- layout_matrix[upper.tri(layout_matrix)]
    
    for(dp in seq_along(upper_idx)){
      region_name_1 <- xregion_name[1,dp]
      region_name_2 <- xregion_name[2,dp]
      
      screen(upper_idx[dp])
      wtx <- xwtc_li[[dp]]$wtx
      biwavelet::plot.biwavelet(wtx,xaxt = "n",yaxt = "n",xlab = "Epoch time",ylab = "Frequency(Hz)")
      
      if(Method == "xwt"){
        title(paste("Cross-wavelet","in subject",subject_index))
      }else if(Method == "wcoh"){
        title(paste("Wavelet coherence","in subject",subject_index))
      }
      
      axis.locs <- log2(1/c(seq(1000,1,by=-10),seq(9,0,by = -1)))
      yticklab <- format(c(seq(1000,1,by=-10),seq(9,0,by = -1)),digits = 1)
      axis(side = 2, at = axis.locs, labels = yticklab)
      
      axis.locs <- axTicks(side = 1)
      xticklab <- format(axis.locs*1000,digits = 1)
      axis(side = 1, at = axis.locs, labels = xticklab)
      
      # abline(v=epoch_time[event_array[,1]==1][1]/1000,lwd=3,col="grey",lty=3)
      # abline(v=median(epoch_time[which(event_array[,2]!=0)])/1000,lwd=3,col="grey",lty=3) 
      
      abline(v=epoch_time[event_array[,event_idx[j,1]]==1][1]/1000,lwd=3,col="grey",lty=3)
      abline(v=median(epoch_time[event_array[,event_idx[j,2]]!=0])/1000,lwd=3,col="grey",lty=3)
      
    }
    
    lower_idx <- layout_matrix[lower.tri(layout_matrix)]
    
    for(dp in seq_along(lower_idx)){
      region_name_1 <- xregion_name[1,dp]
      region_name_2 <- xregion_name[2,dp]
      
      title1 <- paste("Wavelet coherence of","subject",subject_index,"\n")
      ridx <- which(layout_matrix==lower_idx[dp],arr.ind = TRUE)[1]
      cidx <- which(layout_matrix==lower_idx[dp],arr.ind = TRUE)[2]
      
      screen(lower_idx[dp])
      
      if(Method == "xwt"){
        matplot(epoch_time,xwtc_array[ridx,cidx,,],type="l",main = title1,xlab="epoch time",ylim=c(0,1),ylab="cross-wavelet spectrum")
        legend("topright",legend=paste(round(freq),"Hz"),col=seq_along(freq),lty=seq_along(freq),bty="n")
      }else if(Method == "wcoh"){
        matplot(epoch_time,xwtc_array[ridx,cidx,,],type="l",main = title1,xlab="epoch time",ylim=c(0,1),ylab="coherence")
        legend("topright",legend=paste(round(freq),"Hz"),col=seq_along(freq),lty=seq_along(freq),bty="n")
      }
      
      abline(v=epoch_time[event_array[,event_idx[j,1]]==1][1],lwd=3,col="grey",lty=3)
      abline(v=median(epoch_time[event_array[,event_idx[j,2]]!=0]),lwd=3,col="grey",lty=3)
      
      
      grid()
    }
    
    middle_idx <- layout_matrix[diag(layout_matrix)]
    for(dp in seq_along(middle_idx)){
      screen(middle_idx[dp])
      plot(x = 0:1,                   # Create empty plot
           y = 0:1,
           ann = F,
           bty = "n",
           type = "n",
           xaxt = "n",
           yaxt = "n")
      text(x = 0.5,                   # Add text to empty plot
           y = 0.5,
           dimnames(xwtc_array)[[1]][dp], 
           cex = 1.8)
      
    }
  }
  
  if(plot_type == "ht"){
    xwtc_li <- xwtc_obj$xwtc
    xwtc_array <- xwtc_obj$xwtc_array
    epoch_time <- as.numeric(dimnames(xwtc_array)[[3]])
    
    file_path_2 <- paste0(paste(gsub(".png","",basename(file_path)),plot_type,sep="_"),".png")
    png(file_path_2,width = 2000,height = 2000,type="cairo")
    
    p <- par(cex=0.4)
    
    tseq <- seq(head(epoch_time,1)/2,tail(epoch_time,1) / 2, by=100)
    st <- tseq[1:(length(tseq)-1)]
    ed <- tseq[2:(length(tseq))]
    
    n1 <- length(freq) + 1
    n2 <- length(st) + 1
    
    idx <- split.screen(c(n1,n2))
    layout_matrix <- matrix(idx,n1,n2,byrow = TRUE)
    
    for(f in 2:nrow(layout_matrix)){
      for(dp in 2:ncol(layout_matrix)){
        
        screen(layout_matrix[f,dp])
        
        tseq_idx <- st[dp - 1] <= epoch_time & epoch_time < ed[dp - 1]
        aggr_xwtc <- apply(xwtc_array[,,tseq_idx,f - 1],c(1,2),mean)
        cols <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdYlBu")))(100)
        
        
        graphics::image(aggr_xwtc[nrow(aggr_xwtc):1,],
                        xlab="",
                        ylab="",
                        axes = FALSE,
                        col = cols)
        
        title(paste("epoch time",st[dp - 1],":",ed[dp - 1],"(msec)","@",freq[f - 1],"Hz"))
        
        axis.locs <- seq(0,1,length.out = 5)
        xticklab <- format(rownames(aggr_xwtc))
        axis(side = 1,at = axis.locs,labels = xticklab,tick = FALSE)
        
        axis.locs <- seq(1,0,length.out = 5)
        yticklab <- format(rownames(aggr_xwtc))
        axis(side = 2,at = axis.locs,labels = yticklab,tick = FALSE,las=2)
        
        # if(dp==1){
        #   legend.col(col = cols, lev = c(0,1))
        # }
        # legend.col(col = cols, lev = c(0,1))
      }
    }
    
    label_ridx <- layout_matrix[,1]
    for(dp in 2:length(label_ridx)){
      screen(label_ridx[dp])
      plot(x = 0:1,                   # Create empty plot
           y = 0:1,
           ann = F,
           bty = "n",
           type = "n",
           xaxt = "n",
           yaxt = "n")
      text(x = 0.5,                   # Add text to empty plot
           y = 0.5,
           paste(freq[dp - 1],"Hz"),
           cex = 10)
      
    }
    
    label_cidx <- layout_matrix[1,]
    for(dp in 2:length(label_cidx)){
      screen(label_cidx[dp])
      plot(x = 0:1,                   # Create empty plot
           y = 0:1,
           ann = F,
           bty = "n",
           type = "n",
           xaxt = "n",
           yaxt = "n")
      text(x = 0.5,                   # Add text to empty plot
           y = 0.5,
           paste("t =",st[dp - 1],":",ed[dp - 1]),
           cex = 10)
      
    }
    
  }
  
  close.screen(all.screens = TRUE)
  par(p)
  dev.off()
  
}

plot_wt_simple_all <- function(wt_obj,log=FALSE,norm = TRUE,file_path = "plot.png"){
  
  cat("Aggregating results over subjects ... ")
  nobj <- length(wt_obj)
  tmp_obj <- wt_obj[[1]]
  tmp_power_array <- tmp_obj$xwtc_array
  tmp_event_array <- tmp_obj$event_array
  
  power_array_all <- array(NA,dim = c(dim(tmp_power_array),nobj))
  event_array_all <- array(NA,dim = c(dim(tmp_event_array),nobj))
  
  region_indices <- tmp_obj$region
  region_names <- tmp_obj$region_name
  freq <- tmp_obj$freq
  epoch_time <- tmp_obj$epoch_time
  event_type <- tmp_obj$event_type
  
  for(i in seq_along(wt_obj)){
    power_array_all[,,,i] <- wt_obj[[i]]$xwtc_array
    event_array_all[,,i] <- wt_obj[[i]]$event_array
  }
  
  file_path_2 <- paste0(paste(gsub(".png","",basename(file_path)),sep="_"),".png")
  png(file_path_2,width = 2000,height = 2000,type="cairo")
  
  layout_idx <- split.screen(c(length(freq),length(region_indices)))
  layout_mat <- matrix(layout_idx,nrow = length(freq),ncol = length(region_indices),byrow = TRUE)
  
  for(f in seq_along(freq)){
    for(m in seq_along(region_indices)){
      
      screen(layout_mat[f,m])
      p <- par(cex = 0.5,mar = c(3,3,3,1))
      xx <- power_array_all[m,,f,]
      
      if(log){
        xx <- log10(xx)
      }
      
      if(norm){
        #xx <- apply(xx,2,function(x)x/max(x))
        xx <- scale(xx)
      }
      
      event_idx <- matrix(c(1,2,3,4),byrow = TRUE,ncol = 2)
      
      matplot(epoch_time,xx,xlab = "epoch time",ylab="wavelet power",type = "l",lty = 1:ncol(xx),col = 1:ncol(xx))
      legend("topleft",
             legend = paste("Subject",1:nobj),
             lty = 1:nobj,
             col = 1:nobj,
             bty = "n")
      
      title(paste(freq[f],"Hz","of",region_names[m],"wavelet power among subjects"))
      
      abline(v=epoch_time[event_array_all[,event_idx[event_type,1],i]==1][1],lwd=3,col="grey",lty=3)
      abline(v=median(epoch_time[event_array_all[,event_idx[event_type,2],i]!=0]),lwd=3,col="grey",lty=3)
      
      grid()
      par(p)
    }
  }
  
  close.screen(all.screens = TRUE)
  dev.off()
}

plot_xwtc_simple_all <- function(xwtc_obj,plot_type = "ts",file_path = "plot.png"){
  
  
  nobj <- length(xwtc_obj)
  tmp_obj <- xwtc_obj[[1]]
  tmp_power_array <- tmp_obj$xwtc_array
  tmp_event_array <- tmp_obj$event_array
  
  power_array_all <- array(NA,dim = c(dim(tmp_power_array),nobj))
  event_array_all <- array(NA,dim = c(dim(tmp_event_array),nobj))
  
  xregion_indices <- tmp_obj$xregion
  xregion_names <- tmp_obj$xregion_name
  freq <- tmp_obj$freq
  epoch_time <- tmp_obj$epoch_time
  event_type <- tmp_obj$event_type
  Method <- tmp_obj$Method
  
  
  for(i in seq_along(xwtc_obj)){
    power_array_all[,,,,i] <- xwtc_obj[[i]]$xwtc_array
    event_array_all[,,i] <- xwtc_obj[[i]]$event_array
  }
  
  event_idx <- matrix(c(1,2,3,4),byrow = TRUE,ncol = 2)
  
  if(plot_type == "ts"){
    n_pair <- ncol(xregion_indices)
    n_freq <- length(freq)
    
    for(f in 1:n_freq){
      
      file_path_2 <- file.path(dirname(file_path),paste0(paste(gsub(".png","",basename(file_path)),"freq",freq[f],sep="_"),".png"))
      png(file_path_2,width = 2000,height = 2000,type="cairo")
      
      layout_idx <- split.screen(dim(power_array_all)[1:2])
      layout_matrix <- matrix(layout_idx,dim(power_array_all)[1],byrow = TRUE)
      upper_idx <- layout_matrix[upper.tri(layout_matrix)]
      
      for(i in seq_along(upper_idx)){
        screen(upper_idx[i])
        p <- par(mar=c(3,3,1,1),cex = 1)
        
        idx1 <- xregion_indices[1,i]
        idx2 <- xregion_indices[2,i]
        
        heel <- epoch_time[event_array_all[,event_idx[event_type,1],1] == 1]
        foemax <- apply(event_array_all[,event_idx[event_type,2],],2,max)
        foe <- sapply(1:nobj,function(ii)median(epoch_time[event_array_all[,event_idx[event_type,2],ii] == foemax[ii]]))
        foe <- median(foe)
        
        matplot(epoch_time,power_array_all[idx1,idx2,,f,],type="l",lwd=2,ylim=c(0,1),col=1:nobj,lty=1:nobj)
        abline(v=heel,lwd=3,col="grey",lty=3)
        abline(v=foe,lwd=3,col="grey",lty=3)
        
        par(p)
      }
      
      middle_idx <- layout_matrix[diag(layout_matrix)]
      
      for(i in seq_along(middle_idx)){
        screen(middle_idx[i])
        p <- par(mar=c(3,3,1,1),cex = 1.2)
        
        plot(x = 0:1,                   # Create empty plot
             y = 0:1,
             ann = F,
             bty = "n",
             type = "n",
             xaxt = "n",
             yaxt = "n")
        
        if(i == 1){
          legend(x = 0,
                 y = 1,
                 legend = paste("sbj",1:nobj,sep="_"),
                 col = 1:nobj, lty = 1:nobj,lwd=2,bty = "n",cex=0.7)
          text(x = 0.7,                   # Add text to empty plot
               y = 0.5,
               paste0(dimnames(tmp_power_array)[[1]][i],"\n",freq[f],"Hz"))
        }else{
          text(x = 0.5,                   # Add text to empty plot
               y = 0.5,
               paste0(dimnames(tmp_power_array)[[1]][i],"\n",freq[f],"Hz"))
          
        }
      }
      par(p)
      close.screen(all.screens = TRUE)
      
      dev.off()
      
      #par(p)
    }
  }
  #  }
  
  if(plot_type == "ht"){
    xwtc_li <- xwtc_obj$xwtc
    xwtc_array <- xwtc_obj$xwtc_array
    epoch_time <- as.numeric(dimnames(xwtc_array)[[3]])
    
    file_path_2 <- paste0(paste(gsub(".png","",basename(file_path)),plot_type,sep="_"),".png")
    png(file_path_2,width = 2000,height = 2000,type="cairo")
    
    p <- par(cex=0.4)
    
    tseq <- seq(head(epoch_time,1)/2,tail(epoch_time,1) / 2, by=100)
    st <- tseq[1:(length(tseq)-1)]
    ed <- tseq[2:(length(tseq))]
    
    n1 <- length(freq) + 1
    n2 <- length(st) + 1
    
    idx <- split.screen(c(n1,n2))
    layout_matrix <- matrix(idx,n1,n2,byrow = TRUE)
    
    for(f in 2:nrow(layout_matrix)){
      for(dp in 2:ncol(layout_matrix)){
        
        screen(layout_matrix[f,dp])
        
        tseq_idx <- st[dp - 1] <= epoch_time & epoch_time < ed[dp - 1]
        aggr_xwtc <- apply(xwtc_array[,,tseq_idx,f - 1],c(1,2),mean)
        cols <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdYlBu")))(100)
        
        
        graphics::image(aggr_xwtc[nrow(aggr_xwtc):1,],
                        xlab="",
                        ylab="",
                        axes = FALSE,
                        col = cols)
        
        title(paste("epoch time",st[dp - 1],":",ed[dp - 1],"(msec)","@",freq[f - 1],"Hz"))
        
        axis.locs <- seq(0,1,length.out = 5)
        xticklab <- format(rownames(aggr_xwtc))
        axis(side = 1,at = axis.locs,labels = xticklab,tick = FALSE)
        
        axis.locs <- seq(1,0,length.out = 5)
        yticklab <- format(rownames(aggr_xwtc))
        axis(side = 2,at = axis.locs,labels = yticklab,tick = FALSE,las=2)
        
        # if(dp==1){
        #   legend.col(col = cols, lev = c(0,1))
        # }
        # legend.col(col = cols, lev = c(0,1))
      }
    }
    
    label_ridx <- layout_matrix[,1]
    for(dp in 2:length(label_ridx)){
      screen(label_ridx[dp])
      plot(x = 0:1,                   # Create empty plot
           y = 0:1,
           ann = F,
           bty = "n",
           type = "n",
           xaxt = "n",
           yaxt = "n")
      text(x = 0.5,                   # Add text to empty plot
           y = 0.5,
           paste(freq[dp - 1],"Hz"),
           cex = 10)
      
    }
    
    label_cidx <- layout_matrix[1,]
    for(dp in 2:length(label_cidx)){
      screen(label_cidx[dp])
      plot(x = 0:1,                   # Create empty plot
           y = 0:1,
           ann = F,
           bty = "n",
           type = "n",
           xaxt = "n",
           yaxt = "n")
      text(x = 0.5,                   # Add text to empty plot
           y = 0.5,
           paste("t =",st[dp - 1],":",ed[dp - 1]),
           cex = 10)
      
    }
    
  }
  
  close.screen(all.screens = TRUE)
  par(p)
  dev.off()
  
}


synergy <- function(obj,
                    region_indices = c(1,2,3,4,5),
                    subject_index = 1,
                    event_type = 1,
                    freq = c(5,10,15,20),
                    freq_range = NULL,
                    interactive = TRUE,
                    filter_method = c("btw","wt"),
                    syn_method = c("pca","fact","nmf"),
                    n_syn = 5,
                    plotit = TRUE,
                    file_path = NULL){
  
  if(interactive){
    cat("subject_index:\n",seq_along(obj))
    subject_index <- as.numeric(readline("Select subject index:"))
  }
  
  i <- subject_index
  
  # select i-th object (i-th file/subject)
  obj_i <- obj[[i]]
  ff <- obj_i$frequency
  if(interactive){
    cat("event_type:\n",paste0(seq_along(obj_i$epoch),":",names(obj_i$epoch),"\n"))
    event_type <- as.numeric(readline("Select event_type:"))
  }
  
  j <- event_type
  
  # select j-th event (e.g. right/left)
  raw_array <- obj_i$epoch[[j]]
  epoch_time <- raw_array[,1]
  event_array <- raw_array[,grep(".event",dimnames(raw_array)[[2]])]
  signal_array <- raw_array[,colnames(obj_i$signal)]
  event_type_name <- names(obj_i$epoch)[j]
  
  tmp_signal_array_k_1 <- signal_array[,1]
  tmp_xt_1 <- cbind(epoch_time / 1000, tmp_signal_array_k_1)
  tmp_wtx <- wt(tmp_xt_1)
  
  
  cat("Approximating specified frequency to obtain wavelet frequency ... \n")
  per2freq <- data.frame(index = paste("periodID",seq_along(tmp_wtx$period),sep="_"), period = tmp_wtx$period,freq = 1/tmp_wtx$period)
  if(is.null(freq_range)){
    idx <- sapply(freq,function(xx)which.min((xx - per2freq[,3])^2))
    cat("The closest wavelet frequency was: ",paste(round(per2freq[idx,3],2),"Hz (in terms of period:",round(per2freq[idx,2],2),")"),sep="\n")
  }else{
    idx <- which(freq_range[1] < per2freq[,3]  & per2freq[,3] <= freq_range[2])
    freq <- round(per2freq[idx,3],1)
    cat("The closest wavelet frequency was: ",paste(round(per2freq[idx,3],2),"Hz (in terms of period:",round(per2freq[idx,2],2),")"),sep="\n")
  }
  
  
  
  if(filter_method == "wt"){
    
    xwtc_array <- select_wave <- array(NA,
                                       dim = c(length(region_indices),
                                               nrow(obj_i$epoch[[j]]),
                                               length(idx)),
                                       dimnames = list(region_name_1 = names(obj_i$signal)[region_indices],
                                                       epoch_time = epoch_time,
                                                       frequency = freq))
    
    for(k in seq_along(region_indices)){
      
      k_1 <- region_indices[k]
      region_name_1 <- names(obj_i$signal)[k_1]
      signal_array_k_1 <- signal_array[,k_1]
      xt_1 <- cbind(epoch_time/1000,signal_array_k_1)
      
      wtx <- biwavelet::wt(xt_1)
      power <- t(wtx$power)
      wave <- t(Re(wtx$wave))
      xwtc_array[k,,] <- power[,idx]
      
      select_wave[k,,] <- wave[,idx]
      signal_recov_each <- apply(select_wave,c(1,3),envelope)
      signal_recov_all <-  t(apply(select_wave,c(1,2),sum))
      signal_recov_all <-  apply(signal_recov_all,2,envelope)
      cat(".")
      
    }
    
    region_names <- names(obj_i$signal)[region_indices]
    
    signal_array_filt <- signal_recov_all
    
    if(syn_method == "pca"){
      res <- pcaMethods::pca(signal_array_filt,method = "svd",nPcs =  n_syn)
      
      if(plotit){
        if(!is.null(file_path)){
          png(file_path,width = 2000,height = 2000,type="cairo")
        }
        
        screen_idx <- split.screen(c(1,2))
        screen_idx1 <- split.screen(c(n_syn,1),1)
        screen_idx2 <- split.screen(c(2,1),2)
        
        
        screen(screen_idx2[1])
        cols1 <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))(ncol(signal_array_filt))
        
        matplot(epoch_time,
                signal_array_filt,
                type="l",
                main="Filtered signals",
                ylab="vol",
                lwd=2,
                col = cols1,
                lty = 1:ncol(signal_array_filt))
        grid()
        
        screen(screen_idx2[2])
        cols2 <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(7, "Set1"))(ncol(res@scores))
        
        matplot(epoch_time,
                res@scores,
                type="l",
                main="Synergy basis",
                ylab = "Basis signal",
                lwd=2,
                col = cols2,
                lty = 1:ncol(res@scores))
        
        legend("bottomright",legend = paste("Syn",1:ncol(res@scores)),
               col = cols2,lty = 1:ncol(res@scores),cex = 0.7,lwd = 2,bty = "n")
        
        grid()
        for(b in seq_along(screen_idx1)){
          screen(screen_idx1[b])
          p <- par(mar=c(3,2,2,0),cex = 0.7)
          
          barplot(res@loadings[,b],
                  las=2,
                  col=cols1,
                  main = paste("Synergy",b,paste0("(",round(res@R2[b]*100,1),"%)")))
          
          par(p)
          
        }
        
        close.screen(all.screens = TRUE)
        
        if(!is.null(file_path)){
          dev.off() 
        }
      }
      
      synergy_signal <- res@scores
      synergy_weight <- res@loadings
      misc <- list(R2 = res@R2, R2cum = res@R2cum)
      
      output <- list(synergy_signal = synergy_signal, synergy_weight = synergy_weight, misc = misc)
      return(output)
      
    }else if(syn_method == "fact"){
      
      res <- stats::factanal(signal_array_filt,factors = n_syn,rotation = "varimax",scores = "Bartlett")
      ld <- loadings(res)
      SS <- apply(ld,2,function(x)sum(x^2))
      propvar <- SS / nrow(ld)
      cumprop <- cumsum(propvar)
      res$SS <- SS
      res$R2 <- propvar
      res$cumR2 <- cumprop
      
      if(plotit){
        
        if(!is.null(file_path)){
          png(file_path,width = 2000,height = 2000,type="cairo")
        }
        
        screen_idx <- split.screen(c(1,2))
        screen_idx1 <- split.screen(c(n_syn,1),1)
        screen_idx2 <- split.screen(c(2,1),2)
        
        
        screen(screen_idx2[1])
        cols1 <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))(ncol(signal_array_filt))
        
        matplot(epoch_time,
                signal_array_filt,
                type="l",
                main="Filtered signals",
                ylab="vol",
                lwd=2,
                col = cols1,
                lty = 1:ncol(signal_array_filt))
        grid()
        
        screen(screen_idx2[2])
        cols2 <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(7, "Set1"))(ncol(res$scores))
        
        matplot(epoch_time,
                res$scores,
                type="l",
                main="Synergy basis",
                ylab = "Basis signal",
                lwd=2,
                col = cols2,
                lty = 1:ncol(res$scores))
        
        legend("bottomright",legend = paste("Syn",1:ncol(res$scores)),
               col = cols2,lty = 1:ncol(res$scores),cex = 0.7,lwd = 2,bty = "n")
        
        grid()
        for(b in seq_along(screen_idx1)){
          screen(screen_idx1[b])
          p <- par(mar=c(3,2,2,0),cex = 0.7)
          
          barplot(res$loadings[,b],
                  las=2,
                  col=cols1,
                  main = paste("Synergy",b,paste0("(",round(propvar[b]*100,1),"%)")))
          
          par(p)
          
        }
        close.screen(all.screens = TRUE)
        
        if(!is.null(file_path)){
          dev.off() 
        }
        
      }
      
      synergy_signal <- res$scores
      synergy_weight <- res$loadings
      misc <- list(SS = res$SS, R2 = res$R2, R2cum = res$cumR2)
      
      output <- list(synergy_signal = synergy_signal, synergy_weight = synergy_weight, misc = misc)
      return(output)
      
      
    }else if(syn_method == "nmf"){
      
      res <- NMF::nmf(signal_array_filt,rank= n_syn)
      
      if(plotit){
        if(!is.null(file_path)){
          png(file_path,width = 2000,height = 2000,type="cairo")
        }
        screen_idx <- split.screen(c(1,2))
        screen_idx1 <- split.screen(c(n_syn,1),1)
        screen_idx2 <- split.screen(c(2,1),2)
        
        
        screen(screen_idx2[1])
        cols1 <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))(ncol(signal_array_filt))
        
        matplot(epoch_time,
                signal_array_filt,
                type="l",
                main="Filtered signals",
                ylab="vol",
                lwd=2,
                col = cols1,
                lty = 1:ncol(signal_array_filt))
        abline(v = 0, col="grey",lwd=2,lty = 2)
        
        grid()
        
        screen(screen_idx2[2])
        cols2 <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(7, "Set1"))(ncol(res@fit@W))
        
        matplot(epoch_time,
                res@fit@W,
                type="l",
                main="Synergy basis",
                ylab = "Basis signal",
                lwd=2,
                col = cols2,
                lty = 1:ncol(res@fit@W))
        abline(v = 0, col="grey",lwd=2,lty = 2)
        legend("topright",legend = paste("Syn",1:ncol(res@fit@W)),
               col = cols2,lty = 1:ncol(res@fit@W),cex = 0.7,lwd = 2,bty = "n")
        
        grid()
        for(b in seq_along(screen_idx1)){
          screen(screen_idx1[b])
          p <- par(mar=c(3,2,2,0),cex = 0.7)
          
          barplot(res@fit@H[b,],
                  las=2,
                  col=cols1,
                  main = paste("Synergy",b))
          
          par(p)
          
        }
        close.screen(all.screens = TRUE)
        
        if(!is.null(file_path)){
          dev.off() 
        }
        
      }
      
      synergy_signal <- res@fit@W
      synergy_weight <- res@fit@H
      output <- list(synergy_signal = synergy_signal, synergy_weight = synergy_weight, misc = NULL)
      return(output)
      
    }
    
    
    
    
  }else if(filter_method == "btw"){
    
    # signal_array_filt <- apply(signal_array[,region_indices],2,
    #                            function(x)btw(x,freq = ff,order = 4,bandpass = freq_range))
    # signal_array_filt <- apply(signal_array_filt,2,function(x) 100 * x / max(x))
    
    signal_array_filt <- apply(signal_array,2,function(x) 100 * x / max(x))
    region_names <- names(obj_i$signal)[region_indices]
    
    if(syn_method == "pca"){
      res <- pcaMethods::pca(signal_array_filt,method = "svd",nPcs =  n_syn)
      
      if(plotit){
        if(!is.null(file_path)){
          png(file_path,width = 2000,height = 2000,type="cairo")
        }
        
        screen_idx <- split.screen(c(1,2))
        screen_idx1 <- split.screen(c(n_syn,1),1)
        screen_idx2 <- split.screen(c(2,1),2)
        
        
        screen(screen_idx2[1])
        cols1 <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))(ncol(signal_array_filt))
        
        matplot(epoch_time,
                signal_array_filt,
                type="l",
                main="Filtered signals",
                ylab="vol",
                lwd=2,
                col = cols1,
                lty = 1:ncol(signal_array_filt))
        grid()
        
        screen(screen_idx2[2])
        cols2 <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(7, "Set1"))(ncol(res@scores))
        
        matplot(epoch_time,
                res@scores,
                type="l",
                main="Synergy basis",
                ylab = "Basis signal",
                lwd=2,
                col = cols2,
                lty = 1:ncol(res@scores))
        
        legend("bottomright",legend = paste("Syn",1:ncol(res@scores)),
               col = cols2,lty = 1:ncol(res@scores),cex = 0.7,lwd = 2,bty = "n")
        
        grid()
        for(b in seq_along(screen_idx1)){
          screen(screen_idx1[b])
          p <- par(mar=c(3,2,2,0),cex = 0.7)
          
          barplot(res@loadings[,b],
                  las=2,
                  col=cols1,
                  main = paste("Synergy",b,paste0("(",round(res@R2[b]*100,1),"%)")))
          
          par(p)
          
        }
        
        close.screen(all.screens = TRUE)
        if(!is.null(file_path)){
          dev.off() 
        }
      }
      
      synergy_signal <- res@scores
      synergy_weight <- res@loadings
      misc <- list(R2 = res@R2, R2cum = res@R2cum)
      
      output <- list(synergy_signal = synergy_signal, synergy_weight = synergy_weight, misc = misc)
      return(output)
      
    }else if(syn_method == "fact"){
      
      res <- stats::factanal(signal_array_filt,factors = n_syn,rotation = "varimax",scores = "Bartlett")
      ld <- loadings(res)
      SS <- apply(ld,2,function(x)sum(x^2))
      propvar <- SS / nrow(ld)
      cumprop <- cumsum(propvar)
      res$SS <- SS
      res$R2 <- propvar
      res$cumR2 <- cumprop
      
      if(plotit){
        
        if(!is.null(file_path)){
          png(file_path,width = 2000,height = 2000,type="cairo")
        }
        
        screen_idx <- split.screen(c(1,2))
        screen_idx1 <- split.screen(c(n_syn,1),1)
        screen_idx2 <- split.screen(c(2,1),2)
        
        
        screen(screen_idx2[1])
        cols1 <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))(ncol(signal_array_filt))
        
        matplot(epoch_time,
                signal_array_filt,
                type="l",
                main="Filtered signals",
                ylab="vol",
                lwd=2,
                col = cols1,
                lty = 1:ncol(signal_array_filt))
        grid()
        
        screen(screen_idx2[2])
        cols2 <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(7, "Set1"))(ncol(res$scores))
        
        matplot(epoch_time,
                res$scores,
                type="l",
                main="Synergy basis",
                ylab = "Basis signal",
                lwd=2,
                col = cols2,
                lty = 1:ncol(res$scores))
        
        legend("bottomright",legend = paste("Syn",1:ncol(res$scores)),
               col = cols2,lty = 1:ncol(res$scores),cex = 0.7,lwd = 2,bty = "n")
        
        grid()
        for(b in seq_along(screen_idx1)){
          screen(screen_idx1[b])
          p <- par(mar=c(3,2,2,0),cex = 0.7)
          
          barplot(res$loadings[,b],
                  las=2,
                  col=cols1,
                  main = paste("Synergy",b,paste0("(",round(propvar[b]*100,1),"%)")))
          
          par(p)
          
        }
        close.screen(all.screens = TRUE)
        
        if(!is.null(file_path)){
          dev.off() 
        }
        
      }
      
      synergy_signal <- res$scores
      synergy_weight <- res$loadings
      misc <- list(SS = res$SS, R2 = res$R2, R2cum = res$cumR2)
      
      output <- list(synergy_signal = synergy_signal, synergy_weight = synergy_weight, misc = misc)
      return(output)
      
      
      
    }else if(syn_method == "nmf"){
      
      res <- NMF::nmf(signal_array_filt,rank= n_syn)
      
      if(plotit){
        if(!is.null(file_path)){
          png(file_path,width = 2000,height = 2000,type="cairo")
        }
        screen_idx <- split.screen(c(1,2))
        screen_idx1 <- split.screen(c(n_syn,1),1)
        screen_idx2 <- split.screen(c(2,1),2)
        
        
        screen(screen_idx2[1])
        cols1 <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(8, "Accent"))(ncol(signal_array_filt))
        
        matplot(epoch_time,
                signal_array_filt,
                type="l",
                main="Filtered signals",
                ylab="vol",
                lwd=2,
                col = cols1,
                lty = 1:ncol(signal_array_filt))
        abline(v = 0, col="grey",lwd=2,lty = 2)
        
        grid()
        
        screen(screen_idx2[2])
        cols2 <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(7, "Set1"))(ncol(res@fit@W))
        
        matplot(epoch_time,
                res@fit@W,
                type="l",
                main="Synergy basis",
                ylab = "Basis signal",
                lwd=2,
                col = cols2,
                lty = 1:ncol(res@fit@W))
        abline(v = 0, col="grey",lwd=2,lty = 2)
        legend("topright",legend = paste("Syn",1:ncol(res@fit@W)),
               col = cols2,lty = 1:ncol(res@fit@W),cex = 0.7,lwd = 2,bty = "n")
        
        grid()
        for(b in seq_along(screen_idx1)){
          screen(screen_idx1[b])
          p <- par(mar=c(3,2,2,0),cex = 0.7)
          
          barplot(res@fit@H[b,],
                  las=2,
                  col=cols1,
                  main = paste("Synergy",b))
          
          par(p)
          
        }
        close.screen(all.screens = TRUE)
        
        if(!is.null(file_path)){
          dev.off() 
        }
        
      }
      
      synergy_signal <- res@fit@W
      synergy_weight <- res@fit@H
      
      output <- list(synergy_signal = synergy_signal, synergy_weight = synergy_weight, misc = NULL)
      return(output)
      
    }
    
    
  }
  
}


synergy_all <- function(obj,
                        region_indices = c(1,2,3,4,5),
                        freq = c(5,10,15,20),
                        freq_range = c(5,20),
                        event_type = 1,
                        filter_method = c("btw","wt"),
                        syn_method = c("pca","fact","nmf"),
                        n_syn = 5){
  
  obj_i <- obj[[1]]
  raw_array <- obj_i$epoch[[event_type]]
  epoch_time <- raw_array[,1]
  event_array <- raw_array[,grep(".event",dimnames(raw_array)[[2]])]
  signal_array <- raw_array[,colnames(obj_i$signal)]
  event_type_name <- names(obj_i$epoch)[event_type]
  region_names <- colnames(signal_array)[region_indices]
  
  
  
  synergy_obj <- vector("list",length(obj))
  cat("Filtering signals")
  
  for(i in seq_along(obj)){
    synergy_obj[[i]] <- synergy(obj = obj,
                                region_indices = region_indices,
                                subject_index = i,
                                event_type = event_type,
                                freq = freq,
                                freq_range = freq_range,
                                interactive = FALSE,
                                filter_method = filter_method,
                                syn_method = syn_method,
                                n_syn = n_syn,
                                plotit = FALSE,
                                file_path = FALSE)
  }
  
  synergy_signal_array <- array(NA, dim = c(length(epoch_time),
                                            length(obj),
                                            n_syn),
                                dimnames = list(epoch_time,
                                                names(obj),
                                                paste("Synergy",1:n_syn)
                                ))
  
  synergy_weight_array <- array(NA, dim = c(length(region_indices),
                                            n_syn,
                                            length(obj)),
                                dimnames = list(region_names,
                                                paste("Synergy",1:n_syn),
                                                names(obj)
                                ))
  
  if(syn_method %in% c("pca","fact")){
    R2_array <- array(NA, dim = c(n_syn,
                                  length(obj)),
                      dimnames = list(paste("Synergy",1:n_syn),
                                      names(obj)))
    
    for(i in seq_along(synergy_obj)){
      R2_array[,i] <- round(synergy_obj[[i]]$misc$R2*100,2)
    }
    R2_array_mean <- apply(R2_array,1,mean)
    R2_array_sd <- apply(R2_array,1,median)
  }
  
  
  for(i in seq_along(synergy_obj)){
    synergy_signal_array[,i,] <- synergy_obj[[i]]$synergy_signal
    synergy_weight_array[,,i] <- synergy_obj[[i]]$synergy_weight
  }
  
  synergy_signal_array_mean <- apply(synergy_signal_array,c(1,3),median)
  synergy_weight_array_mean <- apply(synergy_weight_array,c(1,2),median)
  synergy_weight_array_sd <- apply(synergy_weight_array,c(1,2),sd)
  
  layout_idx <- split.screen(c(n_syn,2))
  layout_matrix <- matrix(layout_idx,nrow = n_syn,ncol = 2,byrow = TRUE)
  
  for(i in 1:dim(synergy_signal_array)[3]){
    
    screen(layout_matrix[i,1])
    p <- par(mar=c(3,2,1,1),cex=0.6)
    matplot(epoch_time,synergy_signal_array[,,i],type = "l",ylab="Synergy basis",main=paste(toupper(syn_method),"Synergy",i,"base signal"))
    matplot(epoch_time,synergy_signal_array_mean[,i],type = "l",ylab="Synergy basis",add = TRUE, lwd = 3)
    grid()
    par(p)
    
    screen(layout_matrix[i,2])
    p <- par(mar=c(3,2,1,1),cex = 0.6)
    bar <- barplot(synergy_weight_array_mean[,i],main=paste(toupper(syn_method),"Synergy",i,"weight of muscles"),las = 2)
    par(p)
    
    # arrows(bar, synergy_weight_array_mean[,i] - synergy_weight_array_sd[,i], 
    #        bar, synergy_weight_array_mean[,i] + synergy_weight_array_sd[,i], code = 3, lwd = 1, angle = 90, length = 0.1)
  }
  
  close.screen(all.screens = TRUE) 
  
}


hilbert <- function(x,freq=freq,thres = thres) {
  
  # 信号長が奇数の時には末尾にゼロを加えて計算
  if (length(x) %% 2 == 1) {
    is.odd.flag <- TRUE
    x <- c(x, 0)
  } else {
    is.odd.flag <- FALSE
  }
  
  x <- lowpass(x,freq = freq,order = 2,thres = thres)
  
  # ヒルベルト変換はフーリエ変換→負領域をゼロに→逆変換で計算
  X <- fft(x, inverse=FALSE)
  H <- c(1, rep(2, length.out=length(X)/2-1),
         1, rep(0, length.out=length(X)/2-1))
  Y <- X * H
  #Y[f >= lp] <- 0+0i
  y <- fft(Y, inverse=TRUE) / length(Y)
  
  # 信号長が奇数の時には末尾につけたゼロの分を削除
  if (is.odd.flag == TRUE) {
    return(y[1:(length(y)-1)])
  } else {
    return(y)
  }
}

envelope <- function(x,freq=2000,thres = 10){
  
  len <- length(x)
  x <- c(rev(x),x,rev(x))
  
  x <- abs(x - mean(x))
  res <- abs(hilbert(x,freq,thres))
  
  res <- res[(len + 1):(2*len)]
  
  return(res)
  
}


btw <- function(x,freq = freq,
                order = 2,
                bandpass = c(5,15)){
  
  len <- length(x)
  x <- c(rev(x),x,rev(x))
  
  norm_freq <- 2 * bandpass / freq
  pf <- signal::butter(order,norm_freq,type = "pass",plane = "z")
  
  filt <- signal::filtfilt(pf,x)
  #  filt <- envelope(filt)
  filt[(len + 1):(2*len)]
}

lowpass <- function(x,freq,
                    order,
                    thres){
  len <- length(x)
  x <- c(rev(x),x,rev(x))
  
  norm_freq <- 2 * thres / freq
  pf <- signal::butter(order,norm_freq,type = "low")
  filt <- signal::filtfilt(pf,x)
  filt[(len + 1):(2*len)]
}

legend_col <- function(col, lev){
  
  opar <- par
  
  n <- length(col)
  
  bx <- par("usr")
  
  box.cx <- c(bx[2] + (bx[2] - bx[1]) / 1000,
              bx[2] + (bx[2] - bx[1]) / 1000 + (bx[2] - bx[1]) / 50)
  box.cy <- c(bx[3], bx[3])
  box.sy <- (bx[4] - bx[3]) / n
  
  xx <- rep(box.cx, each = 2)
  
  par(xpd = TRUE)
  for(i in 1:n){
    
    yy <- c(box.cy[1] + (box.sy * (i - 1)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i)),
            box.cy[1] + (box.sy * (i - 1)))
    polygon(xx, yy, col = col[i], border = col[i])
    
  }
  par(new = TRUE)
  plot(0, 0, type = "n",
       ylim = c(min(lev), max(lev)),
       yaxt = "n", ylab = "",
       xaxt = "n", xlab = "",
       frame.plot = FALSE)
  axis(side = 4, las = 2, tick = FALSE, line = .25)
  par <- opar
}

data_import <- function(file_path,time_column = NULL, ch_annotation = NULL){
  
  is_valid <- data_check(file_path)
  if(!is_valid){stop("Data is invalid")}
  
  freq_check <- frequency_check(ch_annotation)
  freq_idx <- freq_check$index
  target_freq <- freq_check$target_freq
  src_freq <- freq_check$src_freq
  
  cat("Importing dataset and adding annotations...\n\n")
  dataset <- vector("list",length(file_path))
  
  for(i in seq_along(file_path)){
    
    if(!is.null(time_column)){
      
      signal_column <- which(ch_annotation$ch_signal==1) + 1
      event_column <- which(ch_annotation$ch_event==1) + 1
      
      signal_column_name <- ch_annotation$ch_name[ch_annotation$ch_signal==1]
      event_column_name <- ch_annotation$ch_name[ch_annotation$ch_event==1]
      
      signal <- fread(file_path[i],header = FALSE,select = signal_column,col.names = signal_column_name,data.table = FALSE)
      event <- fread(file_path[i],header = FALSE,select = event_column,col.names = event_column_name,data.table = FALSE)
      time <- fread(file_path[i],header = FALSE,select = time_column,col.names = "time_stamp",data.table = FALSE)
      
      time$time_step <- 1:nrow(time)
      #time$time_stamp = seq(1 / target_freq, max(time[,1]),by = 1 / target_freq)
      time$time_stamp = seq(1 / target_freq, nrow(time) / target_freq,by = 1 / target_freq)
      
      cat(basename(file_path[i]),"imported.\n")
      
    }else{
      
      signal_column <- which(ch_annotation$ch_signal==1)
      event_column <- which(ch_annotation$ch_event==1)
      
      signal_column_name <- ch_annotation$ch_name[ch_annotation$ch_signal==1]
      event_column_name <- ch_annotation$ch_name[ch_annotation$ch_event==1]
      
      signal <- fread(file_path[i],header = FALSE,select = signal_column,col.names = signal_column_name,data.table = FALSE)
      event <- fread(file_path[i],header = FALSE,select = event_column,col.names = event_column_name,data.table = FALSE)
      
      #time <- data.frame(time_stamp = seq(1 / frequency, nrow(signal) * (1 / frequency), by= 1 / frequency))
      time <- data.frame(time_stamp = seq(1 / target_freq, nrow(signal) * (1 / target_freq), by= 1 / target_freq))
      time$time_step <- 1:nrow(time)
      #time$time_stamp <- seq(1 / target_freq, max(time[,1]),by= 1 / target_freq)
      
      cat(basename(file_path[i]),"imported.\n")
    }
    
    
    
    if(!is.null(freq_check)){
      
      is_signal <- ch_annotation$ch_signal == 1
      is_event <- ch_annotation$ch_event == 1
      
      is_signal_mod <- freq_idx[is_signal]
      is_event_mod <- freq_idx[is_event]
      
      
      if(any(is_signal_mod)){
        
        signal_target <- signal[,is_signal_mod,,drop=FALSE]
        
        cat("Adjusting frequcency of ...",paste(colnames(signal[,is_signal_mod]),collapse=", "),"\n")
        cat(src_freq[1],"Hz","-->",target_freq[1],"Hz\n")
        
        src_freq <- ch_annotation[freq_idx,,drop = FALSE]$ch_frequency
        target_range <- apply(event_target,2,function(x)range(which(x!=0)))
        
        signal_target_adj　<- sapply(1:ncol(target_range),
                                    function(ii)time_normalize(signal_target[target_range[1,ii]:target_range[2,ii],ii,drop=FALSE],
                                                               time = time$time_stamp[target_range[1,ii]:target_range[2,ii]],
                                                               n = nrow(signal_target))$signal_approx)
        
        replace_idx <- which(is_signal_mod)
        for(j in seq_along(replace_idx)){
          idx <- replace_idx[j]
          signal[,idx] <- signal_target_adj[,j]
        }
        
      }
      
      if(any(is_event_mod)){
        
        event_target <- event[,is_event_mod,drop = FALSE]
        
        cat("Adjusting frequcency of ...",paste(colnames(event[,is_event_mod]),collapse=", "),"\n")
        cat(src_freq[1],"Hz","-->",target_freq[1],"Hz\n")
        
        
        src_freq <- ch_annotation[freq_idx,,drop = FALSE]$ch_frequency
        target_range <- apply(event_target,2,function(x)range(which(x!=0)))
        
        # time_normalize変更
        event_target_adj　<- sapply(1:ncol(target_range),
                                   function(ii)time_normalize(event_target[target_range[1,ii]:target_range[2,ii],ii,drop=FALSE],
                                                              time = time$time_stamp[target_range[1,ii]:target_range[2,ii]],
                                                              n = nrow(event_target))$signal_approx)
        
        replace_idx <- which(is_event_mod)
        for(j in seq_along(replace_idx)){
          idx <- replace_idx[j]
          event[,idx] <- event_target_adj[,j]
        }
        
      }
      
    }
    
    dataset[[i]] <- list(signal = signal,event = event,time = time, frequency = target_freq)
    cat("\n")
  }
  
  cat("\ndata importing done !\n\n")
  cat("Next step is epoching ...\n")
  dataset
  
}

t_col <- function(color, percent = 50, name = NULL) {
  #      color = color name
  #    percent = % transparency
  #       name = an optional name for the color
  
  ## Get RGB values for named color
  rgb.val <- col2rgb(color)
  
  ## Make new color using input color as base and alpha set by transparency
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3],
               maxColorValue = 255,
               alpha = (100 - percent) * 255 / 100,
               names = name)
  
  ## Save the color
  invisible(t.col)
}


