##############################################################################################
#' @title How tall is a tower tall enough

#' @author
#' Ke Xu \email{ke.xu@wisc.edu}

#' @description 
#' Workflow. Answer the question: how tall is a tower tall enough based on LES results from Matthias Suring

#' @param Currently none

#' @return Currently none

#' @references Xu ms3

#' @keywords eddy-covariance, flux tower height, turbulent flux, LES

#' @examples Currently none

#' @seealso Currently none

# changelog and author contributions / copyrights
#   Ke Xu (2017-03-24)
#     original creation
##############################################################################################



##############################################################################################
##PREREQUISITE
# #clean workspace
rm(list=ls())
#clean console
cat("\014")
#print POSIX date/timestamps with fractional seconds
options(digits.secs=3)

wd <- "~/eddy/"
setwd(wd)
#libraries
#source(paste('./routines/1libraries.r', sep=""))
#functions and constants
#source('./routines/1functions10.r')

library(ff)
library(ffbase)

dirData <- "~/eddy/"
load("~/eddy/in/How_tall_allMthd.RData")

DirOut <- paste0(dirData, "out/allMethod/evenly/")

#name match
mtch <- data.frame(
  "H" = c("w", "pt"),
  "LE" = c("w", "q"),
  stringsAsFactors=FALSE
)


#Just use the mid hour:
tt <- 2

title_list <- c( expression(paste("H [W ", 'm'^-2 * ']', sep="")),
                 expression(paste("LE [W ",'m'^-2 * ']', sep="")),
                 expression(paste('CO'[2], ' flux ' * "\n" * '['  ,mu, 'mol m'^-2 * 's'^-1 * ']', sep="")))

fluxes <- c("H", "LE", "HnLE")
meanFluxSfc <- c(100, 200, 300)
storFlux <- c(5.6, 6.8, 12.4)



#------Read in data in nc files--------------------------------------------------------------------

readNcFiles <- FALSE
if(readNcFiles){
   
  #read in surface forcing data-----------------------
  
  fileFluxSfc <- "ERF_L075km_A50_u5_parallel_xy.nc"
  
  fmtFluxSfc <- open.ncdf(paste0(dirData, "in/", fileFluxSfc))
  print(fmtFluxSfc)
  # [1] "file /sap/kxu/kxu_NEON_140912/LES/ERF_L075km_A50_u5_parallel_xy.nc has 8 dimensions:"
  # [1] "time   Size: 3"
  # [1] "zu_xy   Size: 1"
  # [1] "zw_xy   Size: 1"
  # [1] "zu1_xy   Size: 1"
  # [1] "x   Size: 1729"
  # [1] "xu   Size: 1729"
  # [1] "y   Size: 2305"
  # [1] "yv   Size: 2305"
  # [1] "------------------------"
  # [1] "file /sap/kxu/kxu_NEON_140912/LES/ERF_L075km_A50_u5_parallel_xy.nc has 7 variables:"
  # [1] "double ind_z_xy[zu_xy]  Longname:ind_z_xy Missval:1e+30"
  # [1] "float shf*_xy[x,y,zu1_xy,time]  Longname:shf*_xy Missval:1e+30" [k*m/s]
  # [1] "float qsws*_xy[x,y,zu1_xy,time]  Longname:qsws*_xy Missval:1e+30" [m/s*kg/kg]
  # [1] "float zim_xy[x,y,zu_xy,time]  Longname:zim_xy Missval:1e+30"
  # [1] "float gr_val_xy[x,y,zu_xy,time]  Longname:gr_val_xy Missval:1e+30"
  # [1] "float ziu_xy[x,y,zu_xy,time]  Longname:ziu_xy Missval:1e+30"
  # [1] "float zithr_xy[x,y,zu_xy,time]  Longname:zithr_xy Missval:1e+30"
  
  
  dataSfc <- list()
  dataSfc$t <- get.var.ncdf(fmtFluxSfc, "time")
  dataSfc$H <- get.var.ncdf(fmtFluxSfc, "shf*_xy")
  dataSfc$LE <- get.var.ncdf(fmtFluxSfc, "qsws*_xy")
  
}



#read in tower data---------------------------------
readTowrData <- FALSE
if(readTowrData){
  filesTower <- c("q", "pt", "w") 
  fmtTower <- list()
  dataTower <- list()
  
  for(file in filesTower){
    #file <- "w"
    fmtTower[[file]] <- open.ncdf(paste0(dirData, "in/", "tower_", file, ".nc"))
    dataTower[[file]] <- get.var.ncdf(fmtTower[[file]], paste0(file, "t"))
    if(length(which(dataTower[[file]] == 1e30)) != 0){
      dataTower[[file]][which(dataTower[[file]] == 1e30)] <- NaN
    }
  }
  
  
  print(fmtTower$q)
  # [1] "file /sap/kxu/kxu_NEON_140912/LES/tower_q.nc has 7 dimensions:"
  # [1] "time   Size: 36000"
  # [1] "zu_3d   Size: 18"
  # [1] "zw_3d   Size: 18"
  # [1] "x   Size: 24"
  # [1] "xu   Size: 24"
  # [1] "y   Size: 16"
  # [1] "yv   Size: 16"
  # [1] "------------------------"
  # [1] "file /sap/kxu/kxu_NEON_140912/LES/tower_q.nc has 1 variables:"
  # [1] "float qt[x,y,zw_3d,time]  Longname:qt Missval:1e+30"
  
  
  #read in measurement levels
  fmtTower[[file]] <- open.ncdf(paste0(dirData, "in/", "tower_", file, ".nc"))
  
  dataTower$hgt <- get.var.ncdf(fmtTower[[file]], "zw_3d")
  
  
  
  dataTower$u <- list()
  dataTower$v <- list()
  
}



#------Compare flux from one random tower at different heights with domain mean "truth"-----------
difHgt <- FALSE
if(difHgt){
  load("/sap/kxu/kxu_NEON_140912/LES/out/loadInForLes_alltower.RData")
  
  
  # calculate storage flux
  
  
  dataTower$stor <- list()
  
  
  for (flux in c("H", "LE")){
    #flux <- "H"
    
    dataTower$stor[[flux]] <- array(NaN, dim=c(nrow(dataTower$w), ncol(dataTower$w), length(dataTower$hgt)), dimnames = c("xx", "yy", "hgt"))
    
    for (zz in 1:12){
      #zz <- 2
      
      for(xx in 1:length(dataTower$x)){
        
        for(yy in 1:length(dataTower$y)){
          
          if(zz == 1){
            dataTower$stor[[flux]][xx, yy, zz] <- (mean(dataTower[[mtch[[flux]][2]]][xx, yy, zz, 23001:24000]) - mean(dataTower[[mtch[[flux]][2]]][xx, yy, zz, 12001:13000]))/55/60*dataTower$hgt[zz]
          } else {
            
            dataTower$stor[[flux]][xx, yy, zz] <- dataTower$stor[[flux]][xx, yy, (zz - 1)] + (dataTower$hgt[zz] - dataTower$hgt[zz - 1]) * mean (c((mean(dataTower[[mtch[[flux]][2]]][xx, yy, (zz-1), 23001:24000]) - mean(dataTower[[mtch[[flux]][2]]][xx, yy, (zz-1), 12001:13000]))/55/60,
                                                                                                                                                   (mean(dataTower[[mtch[[flux]][2]]][xx, yy, zz, 23001:24000]) - mean(dataTower[[mtch[[flux]][2]]][xx, yy, zz, 12001:13000]))/55/60))
            
          } # end of if
          
        }# end of loop around yy
        
      }# end of loop around xx
      
      #[mean(last five minute) - mean(first five minute)]/55 min * 49m, only for z=49
      
    }# end of loop around zz
    
  }# end of loop around flux
  
  
  # dataTower$stor$H <- dataTower$stor$H * 100 / mean(dataSfc$H[,,1])
  # 
  # dataTower$stor$LE <- dataTower$stor$LE * 200 / mean(dataSfc$LE[,,1])
  
  
  
  
  
  
  #calculate tower measured NSAE:
  
  
  dataTower$H <- array(NaN, dim=c(length(dataTower$x), length(dataTower$y), length(dataTower$hgt)), dimnames = c("xx", "yy", "hgt"))
  dataTower$LE <- array(NaN, dim=c(length(dataTower$x), length(dataTower$y), length(dataTower$hgt)), dimnames = c("xx", "yy", "hgt"))
  
  
  #calculate turbulent + storage flux for one tower at different levels
  for (hh in 1:12
       #length(hgtTower)
  ){
    #hh <- 1
    
    ##random tower number:
    #nTowr <- 100
    #as.integer(runif(nTowr, min = 1, max = (ncol(dataTower$w) * nrow(dataTower$w))))
    
    for (rr in 1:length(dataTower$x)){
      
      for(cc in 1:length(dataTower$y)){
        
        for (flux in c("H", "LE")){
          #flux <- "H"
          #[w(x,y,t)*T(x,y,t)] - w_bar(x,y) over t dimension * T_bar(x,y) averaged over t dimension
          dataTower[[flux]][rr, cc, hh] <- mean(as.vector(dataTower[[mtch[[flux]][1]]][rr, cc,
                                                                                       hh, (((tt-1) * 12000 + 1):(tt*12000))]) * 
                                                  as.vector(dataTower[[mtch[[flux]][2]]][rr, cc, 
                                                                                         hh, (((tt-1) * 12000 + 1):(tt*12000))])
          ) - mean(as.vector(dataTower[[mtch[[flux]][1]]][rr, cc, 
                                                          hh, (((tt-1) * 12000 + 1):(tt*12000))])) * 
            mean(as.vector(dataTower[[mtch[[flux]][2]]][rr, cc, 
                                                        hh, (((tt-1) * 12000 + 1):(tt*12000))])) 
          #storage flux
          + dataTower$stor[[flux]][rr, cc, hh]
          
          
          
        }# end of loop around flux
      }#end of loop around yy
      
    }#end of loop around xx
  }# end of loop around hh
  
  
  #Unit conversion
  
  dataTower$H <- unlist(dataTower$H) * 100 / mean(dataSfc$H[,,1])
  dataTower$LE <- unlist(dataTower$LE) * 200 / mean(dataSfc$LE[,,1])
  
  #       ##write out calculated tower measured flux
  #       for (flux in c("H", "LE")){
  #         write.csv(dataTower[[flux]], paste0(dirData, "out/", flux, "_of_one_tower_over_hgt_stor.csv"))
  #       }
  
  whr <- c(1:12)
  
  
  #plot the figure of measured flux at one random tower with measuremnt height
  
  
  
  
  #opar<-par()
  png(filename=paste0(dirData, "out/exceptERF/tower observed flux_over_hgt_stor.png"), width = 800, height = 1000, units = "px", pointsize = 20, bg = "white")
  cexvar=2; par(mfrow=c(2,1), las=0, cex.axis=cexvar*0.7, cex.lab=cexvar*0.7, font.lab=2, mar=c(5,5,2,2),
                mgp=c(3,1,0), family="times", lwd=cexvar, cex.main=cexvar*0.7)
  
  
  #define ylim
  for (flux in c("H", "LE")) {
    # flux <- "H"
    
    ymean <- colMeans(dataTower[[flux]], dims=c(2))[whr]
    ysd <- unlist(sapply(whr, function(zz) sqrt(stats::var(as.vector(dataTower[[flux]][,,zz])))))
    
    plot(ymean ~ dataTower$hgt[whr], type="l", #main=paste(substr(dt_lt_lt_l, 1, 19), " - ", hour_out_aft_str, sep=""),
         xlab=paste("Tower height [m]", sep=""), 
         ylab=title_list[which(flux==fluxes)],
         ylim=range(c((ymean + ysd), (ymean - ysd)))
    )
    
    #error bars
    plot_XYer(
      x=dataTower$hgt[whr],
      y=ymean,
      xbar=rep(0,length(dataTower$hgt[whr])),
      ybar=ysd,
      lwd=par()$lwd*0.5,
      lty=1,
      length=0.1
    )
    
    #   plot_XYer(x=x_axis, y=WAVEmn$LE_en, xbar=NULL, ybar=I(100 / 100 * WAVEmn$LE_en), 
    #             lwd=par()$lwd * 0.5, lty=1, length=0.05, xT=F, yT=T)
    
    points(ymean ~ dataTower$hgt[whr], col="red", pch=4)
    
    # abline(h=0, lwd=1)
    abline(h=meanFluxSfc[which(flux==fluxes)], col = "red")
    #abline(lm(WAVEmn[[vari]] ~ x_axis), lty=3)
  }
  
  dev.off()
  
  
  save.image("/sap/kxu/kxu_NEON_140912/LES/out/exceptERF/mulHgt.RData")
  
  
  
  
  #------At two measurement levels, 50 m and 100 m, flux with tower no. increase----------------------
  
  load("/sap/kxu/kxu_NEON_140912/LES/out/loadInForLes_allTowerAtz01.RData")
  
  
  #Calculate storage flux
  
  dataTower$stor <- list()
  for (flux in c("H", "LE")){
    #flux <- "LE"
    
    dataTower$stor[[flux]] <- matrix(NaN, nrow=24, ncol=16)
    
    
    
    for (xx in 1:nrow(dataTower$pt)){
      #xx <- 1
      for (yy in 1:ncol(dataTower$pt)){
        #yy <- 1
        
        #[mean(last five minute) - mean(first five minute)]/55 min * 49m, only for z=49
        dataTower$stor[[flux]][xx, yy] <- (mean(dataTower[[mtch[[flux]][2]]][xx, yy, 23001:24000]) - mean(dataTower[[mtch[[flux]][2]]][xx, yy, 12001:13000]))/55/60*dataTower$hgt
      }
    }
    
  }
  
  
}



#------Calculate flux with multiple towers (increasing number of towers)------------------------
#When with all towers, it is 100 and 200?

ntiti <- 50000

dataMtpTower <- list()

methods <- c("Naive", "S07", "M08")

#Select measurement height
hgtTowerSlct <- c(49
                  #, 126
)


evenly <- TRUE
if(evenly) whrXTowr <- c(1:23)[-c(15, 17, 20)]
if(!evenly) whrXTowr <- c(1:24)


averTime <- "halfHour"


if(averTime == "halfHour") {
  #mtrx <- sapply(c("w", "pt", "q"), function(vv) unlist(sapply(12001:24000, function(ttt) as.vector(dataTower[[vv]][,,ttt]))))
  mtrx <- list()
  mtrx$w <- rbind(unlist(sapply(12001:18000, function(ttt) as.vector(dataTower$w[whrXTowr,,#1,
                                                                                 ttt]))), unlist(sapply(18001:24000, function(ttt) as.vector(dataTower$w[whrXTowr,,
                                                                                                                                                         #1,
                                                                                                                                                         ttt]))))
  mtrx$pt <- rbind(unlist(sapply(12001:18000, function(ttt) as.vector(dataTower$pt[whrXTowr,,#1,
                                                                                   ttt]))), unlist(sapply(18001:24000, function(ttt) as.vector(dataTower$pt[whrXTowr,,#1,
                                                                                                                                                            ttt]))))
  mtrx$q <- rbind(unlist(sapply(12001:18000, function(ttt) as.vector(dataTower$q[whrXTowr,,#1,
                                                                                 ttt]))), unlist(sapply(18001:24000, function(ttt) as.vector(dataTower$q[whrXTowr,,#1,
                                                                                                                                                         ttt]))))
  #mtrx$stor$H <- as.vector(dataTower$stor$H)
  #mtrx$stor$LE <- as.vector(dataTower$stor$LE)

}


if(averTime == "oneHour"){
  #mtrx <- sapply(c("w", "pt", "q"), function(vv) unlist(sapply(12001:24000, function(ttt) as.vector(dataTower[[vv]][,,ttt]))))
  mtrx <- list()
  mtrx$w <- unlist(sapply(12001:24000, function(ttt) as.vector(dataTower$w[whrXTowr,,#1,
                                                                           ttt])))
  mtrx$pt <- unlist(sapply(12001:24000, function(ttt) as.vector(dataTower$pt[whrXTowr,,
                                                                             #1, 
                                                                             ttt])))
  mtrx$q <- unlist(sapply(12001:24000, function(ttt) as.vector(dataTower$q[whrXTowr,,
                                                                           #1, 
                                                                           ttt])))
  mtrx$stor$H <- as.vector(dataTower$stor$H)
  mtrx$stor$LE <- as.vector(dataTower$stor$LE)

}



print(Sys.time())



for (method in methods
                 ){
  
  
  #method <- "S07"
  
  dataMtpTower[[method]] <- list()
  
  
  #How many towrs are used to calculate domain mean flux?
  for (flux in c("H"
                 , "LE"
  )){
    #flux <- "LE"
    
    dataMtpTower[[method]][[flux]] <- list()
    
   # print(Sys.time())
    #for (hh in 1:length(hgtTowerSlct)){
    #hh <- 1
    
    
    
    for (nTowr in 1:14){
      #nTowr <- 14
      
      #randomly choose nTowr towers ntiti times:
      for (titi in 1:ntiti){
        #titi <- 1
        
        
        dataMtpTowerTmp <- list()
        
        if(averTime == "halfHour") {
          whrTowr <- as.integer(runif(nTowr, 1, 20*16))
          
          if(titi%%2 == 1)  whrTowr <- whrTowr + 20*16
        }
        
        
        if(averTime == "oneHour") {
          whrTowr <- as.integer(runif(nTowr, 1, 320))
          
        }
       
        
        #Prepare spatial average w(t), pt(t), and q(t)
        
        
        # way 1: [(w(x,y,t)-w(t))*(T(x,y,t)-T(t))]
        if( method == "Naive"){
          
           # ccTower <- (as.integer((whrTowr-1)/nrow(dataTower$w)) + 1)
           # rrTower <- ((whrTowr-1)%%nrow(dataTower$w) + 1)
        
          SpatAvrg <- list() 
          
          for (vv in c("w", "pt", "q")){
            #vv <- "w"
            SpatAvrg[[vv]] <- sapply(1:ncol(mtrx$w), 
                                     function(ttt)
                                       mean(mtrx[[vv]][whrTowr, 
                                                            #hh, 
                                                            ttt
                                                            ]))
          }
          
          
          
          
          
          #                   #This is just traditional eddy covariance observations
          #                   tmp[[flux]][hh, nTowr, titi] <- mean(tmp[[mtch[[flux]][1]]] * tmp[[mtch[[flux]][2]]])
          #                   - mean(tmp[[mtch[[flux]][1]]]) * mean(tmp[[mtch[[flux]][2]]])
          
          covOver <- "Spac"
          
          if(covOver == "Time") dataMtpTowerTmp$turb <- sapply(whrTowr, function(ww) mean((mtrx$w[ww,] - SpatAvrg[[mtch[[flux]][1]]]) * (mtrx[[mtch[[flux]][2]]][ww,] - SpatAvrg[[mtch[[flux]][2]]]))) #+ 
            #storage flux
            #mean(dataTower$stor[[flux]][whrTowr])
          
          if(covOver == "Spac") dataMtpTowerTmp$turb <- mean(sapply(1:ncol(mtrx$w), function(ttt) mean((mtrx$w[whrTowr,ttt] - SpatAvrg[[mtch[[flux]][1]]][ttt]) * (mtrx[[mtch[[flux]][2]]][whrTowr,ttt] - SpatAvrg[[mtch[[flux]][2]]][ttt]))))
          
          
          
        }# end of method 1
        
        
        
        # way 2: [(w(x,y,t)-w)*(T(x,y,t)-T)] 
        if(method == "S07" | method == "M08"){
          
          
          
          
          #                
          TempSpatAvrg <- list()
          TempSpatAvrg <- sapply(c("w", "pt", "q"), function(vv) mean(mtrx[[vv]][whrTowr,]))
          
          
          
          
          #Prepare all w(x,y,t), pt(x,y,t), and q(x,y,t)
          
          
          #                   #This is just traditional eddy covariance observations
          #                   tmp[[flux]][hh, nTowr, titi] <- mean(tmp[[mtch[[flux]][1]]] * tmp[[mtch[[flux]][2]]])
          #                   - mean(tmp[[mtch[[flux]][1]]]) * mean(tmp[[mtch[[flux]][2]]])
          
          
          
          #;Steinfeld approach
          #                 ;estimate space_theta for all locations
          #                 ;for average all towers
          #                 ;(w * space_theta-space_time_theta)_bar + (w * theta-space_theta) 
          #                 ;no w prime here  
          
          if(method == "S07") {
            dataMtpTowerTmp$turb <- mean((mtrx$w[whrTowr,] - TempSpatAvrg[[mtch[[flux]][1]]]) * (mtrx[[mtch[[flux]][2]]][whrTowr,] - TempSpatAvrg[[mtch[[flux]][2]]])) 
            #storage flux
            #dataMtpTowerTmp$stor <- mean(dataTower$stor[[flux]][whrTowr])
          }
          
          
          #w <- ff::as.ffdf.data.frame(data.frame(w = mtrx$w[whrTowr,]))
          
          #dataMtpTowerTmp[[toto]] <- mean((tmp[[mtch[[flux]][1]]]) * (tmp[[mtch[[flux]][2]]] - TempSpatAvrg[[mtch[[flux]][2]]])) + 
          #storage flux
          #                       dataTower$stor[[flux]][rrTower[toto], ccTower[toto]]
          
          #                 ;Mauder approach
          #                 ;estimate space_time_theta for all sites,
          #                 ; then for first tower location
          #                 ; compute (w_prime*(T-t_avg))_bar  + (localw_bar*(localt_avg-t_avg))
          
          if(method == "M08") {
            dataMtpTowerTmp$turb <- mean((mtrx$w[whrTowr[1],] - mean(mtrx$w[whrTowr[1],])) * (mtrx[[mtch[[flux]][2]]][whrTowr[1],] - TempSpatAvrg[[mtch[[flux]][2]]])) + mean(mtrx$w[whrTowr[1],]) * (mean(mtrx[[mtch[[flux]][2]]][whrTowr[1],]) - TempSpatAvrg[[mtch[[flux]][2]]]) 
            #storage flux
            #dataMtpTowerTmp$stor <- mean(dataTower$stor[[flux]][1])
          }
          
          
          #                     dataMtpTowerTmp[[toto]] <- mean((tmp[[mtch[[flux]][1]]] - mean(tmp[[mtch[[flux]][1]]])) * (tmp[[mtch[[flux]][2]]] - TempSpatAvrg[[mtch[[flux]][2]]])) + mean(tmp[[mtch[[flux]][1]]]) * (mean(tmp[[mtch[[flux]][2]]]) - TempSpatAvrg[[mtch[[flux]][2]]]) + 
          #                       #storage flux
          #                       dataTower$stor[[flux]][rrTower[toto], ccTower[toto]]
          
          
          
          
          
        }# end of method S07, M08
        
        
        
        
        
        #               if(method == 2){
        #                     dataMtpTower[[flux]][hh, nTowr, titi] <- mean(sapply(((tt-1) * 12000 + 1):(tt*12000), 
        #                                                             function(xx)
        #                                                               (mean(dataTower[[mtch[[flux]][1]]][,,hh,xx][whrTowr] * dataTower[[mtch[[flux]][2]]][,,hh,xx][whrTowr]) - mean(dataTower[[mtch[[flux]][1]]][,,hh,xx][whrTowr]) * mean(dataTower[[mtch[[flux]][2]]][,,hh,xx][whrTowr]))
        #                                                             ))
        #                      
        #                 }# end of method 2
        
        #             if(titi == 1){
        
        if(titi == 1){
          dataMtpTower[[method]][[flux]][[toString(nTowr)]]$turb <- dataMtpTowerTmp$turb
          #dataMtpTower[[method]][[flux]][[toString(nTowr)]]$stor <- dataMtpTowerTmp$stor
          
        } else {
          
          dataMtpTower[[method]][[flux]][[toString(nTowr)]]$turb <- c(dataMtpTower[[method]][[flux]][[toString(nTowr)]]$turb, dataMtpTowerTmp$turb)
         # dataMtpTower[[method]][[flux]][[toString(nTowr)]]$stor <- c(dataMtpTower[[method]][[flux]][[toString(nTowr)]]$stor, dataMtpTowerTmp$stor)
          
          #dataMtpTower[[method]][[flux]][[toString(nTowr)]] <- rbind(dataMtpTower[[method]][[flux]][[toString(nTowr)]], dataMtpTowerTmp)
          
        }
        
        
        
        
        
        
        #             } else {
        #               
        #               dataMtpTower$H[hh, nTowr, titi, 1] <- c(dataMtpTower$H[hh, nTowr, titi, 1], mean(as.vector(unlist(dataMtpTowerTmp))))
        #               dataMtpTower$H[hh, nTowr, titi, 2] <- c(dataMtpTower$H[hh, nTowr, titi, 1], sqrt(var(as.vector(unlist(dataMtpTowerTmp)))))
        #               
        #             }
        
      }#end of titi    
      
      
    }# end of 1 to 30 towers
    # }#end of different heights
    #print(Sys.time())
  }# end of H and LE
  
  #       #Unit conversion
  #       dataMtpTower[[method]]$H <- dataMtpTower[[method]]$H * 100 / mean(dataSfc$H[,,1])
  #       dataMtpTower[[method]]$LE <- dataMtpTower[[method]]$LE * 200 / mean(dataSfc$LE[,,1])
  #       
  
  
  
} # end of different methods
print(Sys.time())




#output formatting
OUT <- list()

OUT$mean <- lapply(methods, function(method) lapply(fluxes, function(flux) sapply(1:nTowr, function(x) lapply(c("turb"#, "stor"
                                                                                                                ), function(tblsto) mean(dataMtpTower[[method]][[flux]][[toString(x)]][[tblsto]])))))

names(OUT$mean) <- methods
names(OUT$mean$S07) <- fluxes
names(OUT$mean$M08) <- fluxes
names(OUT$mean$Naive) <- fluxes


#output formatting
OUT <- list()

OUT$mean <- lapply(methods, function(method) lapply(fluxes, function(flux) sapply(1:nTowr, function(x) mean(dataMtpTower[[method]][[flux]][[toString(x)]]$turb))))

names(OUT$mean) <- methods
names(OUT$mean$Naive) <- fluxes
names(OUT$mean$S07) <- fluxes

#names(OUT$mean$M08) <- fluxes


OUT$sd <- lapply(methods, function(method) lapply(fluxes, function(flux) sapply(2:nTowr, function(x) sqrt(var(dataMtpTower[[method]][[flux]][[toString(x)]]$turb)))))
names(OUT$sd) <- methods
names(OUT$sd$Naive) <- fluxes
names(OUT$sd$S07) <- fluxes
names(OUT$sd$M08) <- fluxes



names(OUT$mean) <- methods
names(OUT$mean$S07) <- fluxes
#names(OUT$mean$M08) <- fluxes

OUT1 <- list()
for (method in methods){
  for(flux in fluxes){
    
    OUT1$mean[[method]][[flux]] <- (data.frame(t(OUT$mean[[method]][[flux]])))
    # dimnames(OUT1$mean[[method]][[flux]])[[2]] <- c("turb"#, "stor"
    #                                                 )
    
    # for(tblsto in c("turb"#, "stor"
    #                 )){
      OUT1$mean[[method]][[flux]]#[[tblsto]] 
      <- unlist(OUT1$mean[[method]][[flux]]#[[tblsto]]
                )
      
    #}
  }
}


save.image(paste0(dirData, "out/allMethod/", "S07M08", "times", ntiti, "averTime", averTime, "DstrEvenly", ".RData"))



plot(OUT1$mean$M08$H$turb ~ c(1:30))


for (method in methods){
  OUT[[method]] <- list()
  
  ##write out calculated tower measured flux
  for(flux in c("H", "LE")){
    #flux <- "H"
    
    OUT[[method]][[flux]] <- list()
    
    OUT[[method]][[flux]]$mean <- colMeans(t(dataMtpTower[[method]][[flux]][1,,,1]), na.rm=T)
    OUT[[method]][[flux]]$sd <- #colMeans(t(dataMtpTower[[method]][[flux]][1,,,2]), na.rm=T) 
      sqrt(colVars(t(dataMtpTower[[flux]][1, ,]), na.rm=T))
    OUT[[method]][[flux]]$sd[1] <- 0
    
    
    
    write.csv(as.matrix(OUT[[method]][[flux]]), paste0(dirData, "out/exceptERF/", flux, "_of_muliple_tower_at_", hgtTowerSlct[hh], "m_method", method, ".csv")) 
  }
}


ERF$mean <- c(91.51871, 92.60554, 93.44201, 93.44746, 93.35405, 92.68202, 92.72210, 93.27153, 92.60393, 92.68973, 92.71162, 92.91823, 93.0834, 93.1205)

OUT3$ERF$evenly <- list()
OUT3$ERF$evenly$H <- c(91.51871, 92.60554, 93.44201, 93.44746, 93.35405, 92.68202, 92.72210, 93.27153, 92.60393, 92.68973, 92.71162, 92.91823, 93.0834, 93.1205)
OUT3$ERF$evenly$LE <- c(184.9976, 182.2292, 180.4804, 180.9791, 180.5262, 179.9984, 179.5561, 179.7256, 180.0051, 180.9870, 180.4546, 179.8863, 179.5964, 179.9082)


OUT3$S07$evenly$H <- c(85.65356, 88.98768, 90.22197, 90.58918, 91.15085, 91.35647, 91.42410, 91.57815, 91.63688, 91.84413, 91.78602, 91.93949, 92.11595, 92.17764)
OUT3$S07$evenly$LE <- c(180.7756, 184.1810, 185.2498, 185.7879, 186.0764, 186.4670, 186.5037, 186.6130, 186.7695, 186.8094, 186.8307, 187.0035, 187.0329, 187.0522) 

OUT3$M08$evenly$H <- c(85.87031, 89.33827, 90.23091, 90.93837, 90.37570, 91.36617, 91.41543, 92.11704, 91.65294, 91.82404, 91.70973, 92.35198, 92.07578, 92.13430)
OUT3$M08$evenly$LE <- c(180.8669, 184.3289, 185.8722, 185.7256, 186.5172, 186.1350, 186.2746, 187.2013, 186.7149, 186.7319, 186.9431, 186.9996, 187.0010, 187.0524)

#ERF$mean <- c(85.3467, 85.41277, 87.73395, 89.54257, 89.13310, 90.78794, 92.22323, 92.00530, 91.78583, 90.91024, 92.59321, 91.66941, 91.92674, 92.00218)


ERF$sd <- c(58.4909, 40.64224, 33.24436, 29.30793, 25.51577, 23.51947, 20.99679, 19.36826, 18.22808, 17.87916, 16.21049, 15.89132, 14.99801, 14.81729)


OUT3$ERF$oneHour$LE <- c(184.9976, 182.2292, 180.4804, 180.9791, 180.5262, 179.9984, 179.5561, 179.7256, 180.0051, 180.9870, 180.4546, 179.8863, 179.5964, 179.9082)
  #c(187.6862, 185.4468, 182.4117, 181.9913, 181.9684, 182.4843, 182.0843, 181.9117, 181.4913, 180.9684, 180.4843, 180.4117, 179.9901, 179.6146)

OUT3$ERF$oneHour$H <- c(91.51871, 92.60554, 93.44201, 93.44746, 93.35405, 92.68202, 92.72210, 93.27153, 92.60393, 92.68973, 92.71162, 92.91823, 93.0834, 93.1205)



OUT3$S07$halfHour$HnLE

OUT2 <- list()
OUT2$mean$S07$halfHour <- OUT$mean$S07
OUT2$mean$S07$oneHour <- data.frame(
  H = na.omit(unlist(OUT$mean$S07$H)),
  LE = na.omit(unlist(OUT$mean$S07$H))
)

OUT2$mean$M08$halfHour <- OUT$mean$M08
OUT2$mean$Naive$halfHour <- OUT$mean$Naive 

OUT2$sd$S07$oneHour <- OUT$sd$S07
OUT2$sd$M08$halfHour <- OUT$sd$M08
OUT2$sd$S07$halfHour <- OUT$sd$S07
OUT2$sd$Naive$halfHour <- OUT$sd$Naive


names(OUT2$sd$S07$oneHour)

save(OUT2, file=paste0(dirData,"out/allMethod/evenly/OUT2.RData"))

#unit conversion
OUT3 <- list()

for (mm in names(OUT2$mean)){
  #mm <- names(OUT2$mean)[1]
  for (tt in names(OUT2$mean[[mm]])){
    #tt <- names(OUT2$mean[[mm]])[1]
    for (ff in names(OUT2$mean[[mm]][[tt]])){
      if(ff == "H"){
        OUT3[[mm]][[tt]][[ff]] <- OUT2$mean[[mm]][[tt]][[ff]] * 100/mean(dataSfc$H[,,2])
      }
      
      if(ff == "LE"){
        OUT3[[mm]][[tt]][[ff]] <- OUT2$mean[[mm]][[tt]][[ff]] * 200/mean(dataSfc$LE[,,2])
      }
    }
  }
}






for (mm in names(OUT2$sd)){
  #mm <- names(OUT2$sd)[1]
  for (tt in names(OUT2$sd[[mm]])){
    #tt <- names(OUT2$sd[[mm]])[1]
    for (ff in names(OUT2$sd[[mm]][[tt]])){
      if(ff == "H"){
        OUT3[[mm]][[tt]][[ff]] <- OUT2$sd[[mm]][[tt]][[ff]] * 100/mean(dataSfc$H[,,2])
      }
      
      if(ff == "LE"){
        OUT3[[mm]][[tt]][[ff]] <- OUT2$sd[[mm]][[tt]][[ff]] * 200/mean(dataSfc$LE[,,2])
      }
    }
  }
  
}


#H + LE

for (mm in names(OUT3)){
  for (tt in names(OUT3[[mm]])){
    OUT3[[mm]][[tt]]$HnLE <- OUT3[[mm]][[tt]]$H + OUT3[[mm]][[tt]]$LE
  }
}

save(OUT3, file=paste0(DirOut, "/OUT3.RData"))


#plot out different methods-------------------------------------------------
plot((OUT1$mean$S07$H$turb) ~ c(1:30))
plot((OUT1$mean$M08$H$turb)[1:26] ~ c(1:26))


plotDifMethFluxOverNoTowr <- TRUE
if(plotDifMethFluxOverNoTowr){
  
  
  red <- rgb(228,26,28, max=255)
  
  green <- rgb(77,175,74, max=255)
  
  blue <- rgb(55,126,184, max=255)
  
  org <- rgb(255,127,0, max=255)
  
  prp <- rgb(152,78,163, max=255)
  
  cols <- list(
    "S07" = green,
    "M08" = prp,
    "ERF" = red
  )
  
  unitConv <- data.frame(
    H = 1,
    LE = 2,
    HnLE = 3
  )
  
    
    #cc <- 14
    
    
    
  whr <- 1:14
  
  cc <- 14
  
  ftxt <- c(
    "H" , "LE", "H + LE"
    
  )
  
  flagPlotEvenly <- TRUE
  
  #plot flux over no. of towers at two different heights
  #opar<-par()
  for (flux in c("H", "LE"#, "HnLE"
                 )){
    #flux <- "LE"
    
  png(filename=paste0(DirOut, "/MtpTowrFlux_", flux, "_1.png"), width = 800, height = 800, units = "px", pointsize = 20, bg = "white")
  cexvar=2; par(mfrow=c(1,1), las=0, cex.axis=cexvar, cex.lab=cexvar, font.lab=2, mar=c(5,5,2,2),
                mgp=c(3,1,0), family="times", lwd=cexvar*1.3, cex.main=cexvar*1.2)
  
  
  #define ylim
  
      plot(OUT3$S07$halfHour[[flux]][whr]/unitConv[[flux]] ~ c(1:30)[whr], type="l", 
           lwd =6,
           #main=paste(substr(dt_lt_lt_l, 1, 19), " - ", hour_out_aft_str, sep=""),
           xlab=paste("number of towers used", sep=""),
           ylab=expression(bold("Ensemble domain mean energy [%]")), main = ftxt[which(flux == c("H", "LE", "HnLE"))],
           ylim=range(OUT3$Naive$halfHour[[flux]][2:whr[length(whr)]]/unitConv[[flux]],
                      OUT3$S07$halfHour[[flux]]/unitConv[[flux]],
                      #OUT3$S07$oneHour[[flux]]/unitConv[[flux]],
                      OUT3$M08$oneHour[[flux]]/unitConv[[flux]],
                      OUT3$ERF$oneHour[[flux]]/unitConv[[flux]],
                      OUT3$ERF$evenly[[flux]]/unitConv[[flux]],
                      (meanFluxSfc[which(flux==c("H", "LE", "HnLE"))] - storFlux[which(flux==c("H", "LE", "HnLE"))])/unitConv[[flux]], na.rm=TRUE),
           col = #"white"
           green#, lty = 2
      )
      
    # lines(x=c(1:30)[whr], y=OUT3$S07$halfHour[[flux]][whr]
    #       #rowMeans(dataMtpTower[[flux]][2,,], na.rm=T)[whr]
    #       , col = green)
    
    if(!flagPlotEvenly){
       points(OUT3$S07$halfHour[[flux]][whr]/unitConv[[flux]] ~ c(1:30)[whr], col=#"black"
             green
           , pch=
           21
    )
    }
   
    
   
     #Naive
   
     # plot(OUT3$Naive$halfHour[[flux]][whr] ~ c(1:30)[whr], type="l", #main=paste(substr(dt_lt_lt_l, 1, 19), " - ", hour_out_aft_str, sep=""),
     #      xlab=paste("number of towers used", sep=""),
     #      ylab=paste0("ensemble mean of ", flux, " [W m-2]"),
     #      ylim=range(OUT3$Naive$halfHour[[flux]][2:whr[length(whr)]],
     #                 OUT3$S07$halfHour[[flux]],
     #                 #OUT3$S07$oneHour[[flux]],
     #                 OUT3$M08$halfHour[[flux]],
     #                 OUT3$ERF$oneHour[[flux]],
     #                 (meanFluxSfc[which(flux == fluxes)] - 6.4) ,
     #                 na.rm=TRUE),
     #      col = org
     # )
      if(!flagPlotEvenly){
        lines(x=c(1:30)[whr[-1]], y=OUT3$Naive$halfHour[[flux]][whr[-1]]/unitConv[[flux]]
              #rowMeans(dataMtpTower[[flux]][2,,], na.rm=T)[whr]
              , col =  org, lwd =6
              )

        points(OUT3$Naive$halfHour[[flux]][whr[-1]]/unitConv[[flux]] ~ c(1:30)[whr[-1]], col= org#, pch=21
        )
      }
     
      
      #method <- "S07"
      if(TRUE){
      lines(x=c(1:30)[whr], y=OUT3$S07$oneHour[[flux]][whr]/unitConv[[flux]]
            #rowMeans(dataMtpTower[[flux]][2,,], na.rm=T)[whr]
            , col = blue)
      
        if(!flagPlotEvenly){
           points(
        #rowMeans(dataMtpTower[[flux]][2,,], na.rm=T)
        OUT3$S07$oneHour[[flux]][whr]/unitConv[[flux]] ~ c(1:30)[whr], col = blue, pch=21)
      
        }
     
      
      }
      
      
      #method <- "M08"
      
      lines(x=c(1:30)[whr], y=OUT3$M08$halfHour[[flux]][whr]/unitConv[[flux]]
            #rowMeans(dataMtpTower[[flux]][2,,], na.rm=T)[whr]
            , lwd =6,col =  prp#, lty = 2
            )
      
      if(!flagPlotEvenly){
      points(
        #rowMeans(dataMtpTower[[flux]][2,,], na.rm=T)
        OUT3$M08$halfHour[[flux]][whr]/unitConv[[flux]] ~ c(1:30)[whr], col =  prp
        , pch=21)
      }
      
      
      #method <- "ERF"
      
      lines(x=c(1:30)[whr], y=OUT3$ERF$oneHour[[flux]][whr]/unitConv[[flux]]
            #rowMeans(dataMtpTower[[flux]][2,,], na.rm=T)[whr]
            , lwd =6,col =  red#, lty = 2
            )
      
      
      if(!flagPlotEvenly){
      points(
        #rowMeans(dataMtpTower[[flux]][2,,], na.rm=T)
        OUT3$ERF$oneHour[[flux]][whr]/unitConv[[flux]] ~ c(1:30)[whr], col =  red
        , pch=21)
      }
      
      
      if(flagPlotEvenly){
        
         #S07, M08, ERF with evenly sampling towers
        
        for(mm in c("S07", "M08", "ERF")){
          lines(x=c(1:30)[whr], y=OUT3[[mm]]$evenly[[flux]][whr]/unitConv[[flux]]
            #rowMeans(dataMtpTower[[flux]][2,,], na.rm=T)[whr]
            , lwd =6,col = cols[[mm]]
            )
      
          
        }
      
      }#end of flagPlotEvenly
      
   
      # #one single Tower with S07
      # points(OUT3$S07$halfHour[[flux]][1] ~ c(1:30)[1], col="black"
      #        #green
      #        , pch=19
      #        #21
      # )
      
    
    abline(h=(meanFluxSfc[which(flux==c("H", "LE", "HnLE"))] - storFlux[which(flux==c("H", "LE", "HnLE"))])/unitConv[[flux]], lwd =6,col =  blue
           )
    #abline(lm(WAVEmn[[vari]] ~ x_axis), lty=3)
  
  
  #legend
    if(flux == c("H", "LE", "HnLE")[1]){
      if(FALSE){
        legend(x="bottomright", ncol=1, lty=rep(1,
                                                5
                                                #2
        ), lwd =rep(6,5),
        col=c(org, 
              green, 
              #blue,
              prp, red, blue #"red",cols$M08,"green"#, cols$S07
        ), legend=c("Spatial EC",
                    "S07",
                    # "S07 one Hour",
                    "M08",
                    "ERF",
                    "Reference"
        ),
        #methods, #c(
        #paste0("49 m"), 
        #paste0("126 m")),
        #eval(substitute(expression(sigma * "(res)" == x1), list(x1=round(summary(LM)$sigma,4)))),
        #eval(substitute(expression(R^{2} == x1), list(x1=round(summary(LM)$r.squared,4))))),
        bty="n", cex=cexvar )
      }
      
      if(TRUE){
        legend(x="bottomright", ncol=1, lty=c(2,1,2,1,2,1,1
                                                
                                                #2
        ), lwd =rep(6,5),
        col=c(#org, 
              green, 
              green,
              #blue,
              prp, 
              prp,
              red, red,
              blue #"red",cols$M08,"green"#, cols$S07
        ), legend=c(#"Spatial EC",
                    "S07 with uneven sampling",
                    "S07 with even sampling",
                    # "S07 one Hour",
                    "M08 with uneven sampling",
                    "M08 with even sampling",
                    "ERF with uneven sampling",
                    "ERF with even sampling",
                    "Reference"
        ),
        #methods, #c(
        #paste0("49 m"), 
        #paste0("126 m")),
        #eval(substitute(expression(sigma * "(res)" == x1), list(x1=round(summary(LM)$sigma,4)))),
        #eval(substitute(expression(R^{2} == x1), list(x1=round(summary(LM)$r.squared,4))))),
        bty="n", cex=cexvar*0.6 )
      }
      
    }
  
  
  dev.off()
  
  
  }# end of plot different fluxes
  
  
}#end of if plot or not


save.image(paste0(DirOut, "/How_tall_allMthd_OUT3.RData"))


#plot flux over no. of towers at two different heights
plotFluxOverNoTowrAt2Hgt <- FALSE
if(plotFluxOverNoTowrAt2Hgt){
  whr <- 1:30
  #opar<-par()
  png(filename=paste0(dirData, "out/exceptERF/MtpTowrFlux_at_different_hgt_mthd", method, ".png"), width = 800, height = 1600, units = "px", pointsize = 20, bg = "white")
  cexvar=2; par(mfrow=c(2,1), las=0, cex.axis=cexvar*0.7, cex.lab=cexvar*0.7, font.lab=2, mar=c(5,5,2,2),
                mgp=c(3,1,0), family="times", lwd=cexvar, cex.main=cexvar*0.7)
  
  
  #define ylim
  for (flux in c("H", "LE")){
    #flux <- "H"
    
    ymean <- colMeans(t(dataMtpTower[[flux]][1,,,1]), na.rm=T)
    ysd <- colMeans(t(dataMtpTower[[flux]][1,,,2]), na.rm=T) #sqrt(colVars(t(dataMtpTower[[flux]][1, ,]), na.rm=T))
    ysd[1] <- 0
    
    ##write out calculated tower measured flux
    
    write.csv(cbind(ymean, ysd), paste0(dirData, "out/exceptERF/", flux, "_of_muliple_tower_at_", hgtTowerSlct[hh], "m_method", method, ".csv"))  
    
    
    
    plot(ymean[whr] ~ c(1:30)[whr], type="l", #main=paste(substr(dt_lt_lt_l, 1, 19), " - ", hour_out_aft_str, sep=""),
         xlab=paste("no. of Towers", sep=""), 
         ylab=title_list[which(flux==fluxes)],
         ylim=range(c((ymean + ysd)[whr], (ymean - ysd)[whr])),
         col = "red"
    )
    
    #error bars
    plot_XYer(
      x=c(1:30)[whr],
      y=ymean[whr],
      xbar=rep(0,30)[whr],
      ybar=ysd[whr],
      lwd=par()$lwd*0.5,
      lty=1,
      length=0.1
    )
    
    #   plot_XYer(x=x_axis, y=WAVEmn$LE_en, xbar=NULL, ybar=I(100 / 100 * WAVEmn$LE_en), 
    #             lwd=par()$lwd * 0.5, lty=1, length=0.05, xT=F, yT=T)
    
    points(ymean[whr] ~ c(1:30)[whr], col="red", pch=21)
    
    #       #height 2
    #       lines(x=c(1:30)[whr], y=dataMtpTower[[flux]][2,,,1]
    #               #rowMeans(dataMtpTower[[flux]][2,,], na.rm=T)[whr]
    #        , col = "green")
    #       
    #       points(
    #         #rowMeans(dataMtpTower[[flux]][2,,], na.rm=T)
    #         dataMtpTower[[flux]][2, , ,1][whr] ~ c(1:30)[whr], col="green", pch=21)
    
    
    abline(h=meanFluxSfc[which(flux==fluxes)], col = "black")
    #abline(lm(WAVEmn[[vari]] ~ x_axis), lty=3)
  }
  
  legend(x="bottomright", ncol=1, lty=rep(2,
                                          1
                                          #2
  ), col=c("red"#,"green"
  ), legend=c(
    paste0("49 m")), 
  #paste0("126 m")),
  #eval(substitute(expression(sigma * "(res)" == x1), list(x1=round(summary(LM)$sigma,4)))),
  #eval(substitute(expression(R^{2} == x1), list(x1=round(summary(LM)$r.squared,4))))),
  bty="n", cex=cexvar /2.5)
  
  dev.off()
}

save.image(paste0(dirData, "out/allMethod/How_tall_mthd", method, "evenly.RData"))









