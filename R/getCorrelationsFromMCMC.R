
#str(mcmcProcessed$sims.list)
# for ( i in 1:length(mcmcProcessed$sims.list) ) {
#
#   if( names(mcmcProcessed$sims.list)[[i]] != "grIndRE" &
#       names(mcmcProcessed$sims.list)[[i]] != "lengthExp"){
#     print(c(i,names(mcmcProcessed$sims.list)[[i]]))
#
#     tmp=cor(data.frame(mcmcProcessed$sims.list[[i]]))
#
#     print("gt 0.8")
#     print(which(tmp > 0.8 & tmp != 1, arr.ind = TRUE))
#
#     print("lt 0.8")
#     print(which(tmp < -0.8 & tmp != 1, arr.ind = TRUE))
#   }
# }



