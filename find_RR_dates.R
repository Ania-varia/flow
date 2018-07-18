find_RR_dates <- function(whole_dt, param, by_sex=FALSE, min_ind=70){
  require(data.table)
  #whole_dt  data frame/data.table with columns "Assay_Date", "Gender", param, "Genotype", KO and WT
  #by_sex= whether to split by sex when determining controls 
  #min_ind - how many controls individuals
  #param - on which parameter to base control finding
  if(!is.data.table(whole_dt)){
    whole_dt <-data.table(whole_dt)
  }
  if(by_sex){
    all_res <- rbind(find_RR_dates(whole_dt[Gender=="Female"], param, by_sex=FALSE)[,sex:="female"],
                     find_RR_dates(whole_dt[Gender=="Male"], param, by_sex=FALSE)[,sex:="male"])
    return(all_res)
  }
  
  test_BM <- whole_dt[!is.na(param),c("Assay_Date", "Gender", param, "Genotype"), with=F]
  test_BM_WT <- take_only_WTs(test_BM)
  test_BM_WT[,wt_date:=Assay_Date]
  test_BM_WT  <-test_BM_WT[,.(.N), by=.(Assay_Date,wt_date)]
  test_BM_WT[, cumsum_fromstart:=cumsum(N)]
  test_BM_WT[, cumsum_fromend:=rev(cumsum(rev(N)))]
  setkey(test_BM_WT, Assay_Date)
  #find start indices
  st <- apply(test_BM_WT,1, function(x) {max(which(as.numeric(x["cumsum_fromstart"])- as.numeric(test_BM_WT$cumsum_fromstart) - as.numeric(x["N"])>=ceiling((min_ind-as.numeric(x["N"]))/2)) )})
  st[st<0] <-1
  #find end indices
  end <-apply(test_BM_WT,1, function(x) {min(which(as.numeric(test_BM_WT$cumsum_fromstart) - as.numeric(x["cumsum_fromstart"])>=ceiling((min_ind-as.numeric(x["N"]))/2))) })
  end[end>length(end)] <-length(end)
  test_BM_WT <-cbind(test_BM_WT, st,end)
  #fill start and end dates for days with enough samples before/after
  test_BM_WT[cumsum_fromstart-N >= ceiling((min_ind-N)/2) & cumsum_fromend-N >= ceiling((min_ind-N)/2)
             ,st_date:=test_BM_WT$Assay_Date[st]]
  test_BM_WT[cumsum_fromstart-N>=ceiling((min_ind-N)/2) & cumsum_fromend-N >= ceiling((min_ind-N)/2)
             ,end_date:=test_BM_WT$Assay_Date[end]]
  
  #fill start and end dates for days at the beginning/end
  test_BM_WT[cumsum_fromstart-N < ceiling((min_ind-N)/2),st_date:=test_BM_WT$Assay_Date[1]]
  test_BM_WT[cumsum_fromstart-N < ceiling((min_ind-N)/2), end_date:=test_BM_WT[min(which(cumsum_fromstart>=min_ind)),Assay_Date]]
  
  test_BM_WT[cumsum_fromend-N < ceiling((min_ind-N)/2), st_date:=test_BM_WT[max(which(cumsum_fromend>=min_ind)),Assay_Date]]
  test_BM_WT[cumsum_fromend-N < ceiling((min_ind-N)/2), end_date:=max(Assay_Date)]
  
  #Fix days when there are no WTs
  missing_dates <- as.Date(setdiff(as.character(test_BM$Assay_Date), as.character(test_BM_WT$Assay_Date)))
  misising_BM <- test_BM[Assay_Date%in%missing_dates][, .(miss_date=Assay_Date, .N), by=Assay_Date]
  setkey(misising_BM, Assay_Date)
  all_res <-rbind(test_BM_WT[,c("Assay_Date", "st_date", "end_date"), with=F], test_BM_WT[misising_BM, roll=T][,c("Assay_Date", "st_date", "end_date"), with=F])[order(Assay_Date)]
  all_res[,sex:="both"]
  return(all_res[,param:=param])
}