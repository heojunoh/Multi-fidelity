library(ggplot2)
library(gridExtra)
# install.packages("ggpubr")
library(ggpubr)
# install.packages("dplyr")
library(dplyr)
# install.packages("hrbrthemes")
library(hrbrthemes)
# install.packages("tidyverse")
library(tidyverse)
library(gtable)
library(grid)
# extrafont::font_import() 
extrafont::loadfonts()



###### Nonlinear ######
resultmatc <- matrix(, nrow=101, ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatc[[j]])-1)){
    r <- c(r, rep(rmsematc[[j]][i], costmatc[[j]][i+1] - costmatc[[j]][i]))
  }
  r <- c(r, rmsematc[[j]][length(rmsematc[[j]])])
  resultmatc[,j] <- r[1:101]
}

resultmatc

resultmeanc <- apply(resultmatc, 1, mean)
resultmaxc <- apply(resultmatc, 1, max)
resultminc <- apply(resultmatc, 1, min)



resultmatc2 <- matrix(, nrow=101, ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatc2[[j]])-1)){
    r <- c(r, rep(rmsematc2[[j]][i], costmatc2[[j]][i+1] - costmatc2[[j]][i]))
  }
  r <- c(r, rmsematc2[[j]][length(rmsematc2[[j]])])
  resultmatc2[,j] <- r[1:101]
}

resultmatc2

resultmeanc2 <- apply(resultmatc2, 1, mean)
resultmaxc2 <- apply(resultmatc2, 1, max)
resultminc2 <- apply(resultmatc2, 1, min)



resultmatc3 <- matrix(, nrow=101, ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatc3[[j]])-1)){
    r <- c(r, rep(rmsematc3[[j]][i], costmatc3[[j]][i+1] - costmatc3[[j]][i]))
  }
  r <- c(r, rmsematc3[[j]][length(rmsematc3[[j]])])
  resultmatc3[,j] <- r[1:101]
}

resultmatc3

resultmeanc3 <- apply(resultmatc3, 1, mean)
resultmaxc3 <- apply(resultmatc3, 1, max)
resultminc3 <- apply(resultmatc3, 1, min)



resultmatco <- matrix(, nrow=101, ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatco[[j]])-1)){
    r <- c(r, rep(rmsematco[[j]][i], costmatco[[j]][i+1] - costmatco[[j]][i]))
  }
  r <- c(r, rmsematco[[j]][length(rmsematco[[j]])])
  resultmatco[,j] <- r[1:101]
}

resultmatco

resultmeanco <- apply(resultmatco, 1, mean)
resultmaxco <- apply(resultmatco, 1, max)
resultminco <- apply(resultmatco, 1, min)



resultmatk <- matrix(, nrow=101, ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatk[[j]])-1)){
    r <- c(r, rep(rmsematk[[j]][i], costmatk[[j]][i+1] - costmatk[[j]][i]))
  }
  r <- c(r, rmsematk[[j]][length(rmsematk[[j]])])
  resultmatk[,j] <- r[1:101]
}

resultmatk

resultmeank <- apply(resultmatk, 1, mean)
resultmaxk <- apply(resultmatk, 1, max)
resultmink <- apply(resultmatk, 1, min)


nonlinearalresult <- data.frame(x=rep(0:100,4), 
                                Strategy=factor(c(rep("Comprehensive searching",101), rep("Optimal point-first",101), rep("Cokriging",101), rep("MR-SUR",101)), 
                                                levels=c("Comprehensive searching","Optimal point-first","Cokriging","MR-SUR")),
                                Mean = c(resultmeanc, resultmeanc2, resultmeanco, resultmeank),
                                Max = c(resultmaxc, resultmaxc2, resultmaxco, resultmaxk),
                                Min = c(resultminc, resultminc2, resultminco, resultmink))

plotalnonlinear <- ggplot(nonlinearalresult, aes(x)) +                                     
  geom_line(aes(y=Mean, group=Strategy, color=Strategy, linetype=Strategy), size = 1) +
  geom_ribbon(aes(ymin=Min, ymax=Max, group=Strategy, color=Strategy, fill=Strategy), alpha=0.1, 
              color = "black", linetype = "blank")+
  theme_ipsum() +  ggtitle("Active learning for nonlinear example funcion") +
  theme(axis.title.x = element_text(size=14, margin = margin(t = 10), hjust=0.5),
        # axis.text.x = element_blank(),
        axis.title.y = element_text(size=14, margin = margin(r = 10), hjust=0.5),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="Costs", y = "RMSE") #+ theme(legend.position = c(0.9,0.9))

plotalnonlinear


nonlinearalresult <- data.frame(x=rep(0:100,5), 
                                Strategy=factor(c(rep("ALC",101), rep("ALC+ALM",101), rep("ALM",101), rep("Cokriging",101), rep("MR-SUR",101)), 
                                                levels=c("ALC","ALC+ALM","ALM","Cokriging","MR-SUR")),
                                Mean = c(resultmeanc, resultmeanc2, resultmeanc3, resultmeanco, resultmeank),
                                Max = c(resultmaxc, resultmaxc2, resultmaxc3, resultmaxco, resultmaxk),
                                Min = c(resultminc, resultminc2, resultminc3, resultminco, resultmink))

plotalnonlinear <- ggplot(nonlinearalresult, aes(x)) +                                     
  geom_line(aes(y=Mean, group=Strategy, color=Strategy, linetype=Strategy), size = 1) +
  geom_ribbon(aes(ymin=Min, ymax=Max, group=Strategy, color=Strategy, fill=Strategy), alpha=0.1, 
              color = "black", linetype = "blank")+
  theme_ipsum() +  ggtitle("Active learning for nonlinear example funcion") +
  theme(axis.title.x = element_text(size=14, margin = margin(t = 10), hjust=0.5),
        # axis.text.x = element_blank(),
        axis.title.y = element_text(size=14, margin = margin(r = 10), hjust=0.5),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="Costs", y = "RMSE") #+ theme(legend.position = c(0.9,0.9))

plotalnonlinear



###### Park ######
resultmatc <- matrix(, nrow=101, ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatc[[j]])-1)){
    r <- c(r, rep(rmsematc[[j]][i], costmatc[[j]][i+1] - costmatc[[j]][i]))
  }
  r <- c(r, rmsematc[[j]][length(rmsematc[[j]])])
  resultmatc[,j] <- r[1:101]
}

resultmatc

resultmeanc <- apply(resultmatc, 1, mean)
resultmaxc <- apply(resultmatc, 1, max)
resultminc <- apply(resultmatc, 1, min)



resultmatc2 <- matrix(, nrow=101, ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatc2[[j]])-1)){
    r <- c(r, rep(rmsematc2[[j]][i], costmatc2[[j]][i+1] - costmatc2[[j]][i]))
  }
  r <- c(r, rmsematc2[[j]][length(rmsematc2[[j]])])
  resultmatc2[,j] <- r[1:101]
}

resultmatc2

resultmeanc2 <- apply(resultmatc2, 1, mean)
resultmaxc2 <- apply(resultmatc2, 1, max)
resultminc2 <- apply(resultmatc2, 1, min)



resultmatco <- matrix(, nrow=101, ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatco[[j]])-1)){
    r <- c(r, rep(rmsematco[[j]][i], costmatco[[j]][i+1] - costmatco[[j]][i]))
  }
  r <- c(r, rmsematco[[j]][length(rmsematco[[j]])])
  resultmatco[,j] <- r[1:101]
}

resultmatco

resultmeanco <- apply(resultmatco, 1, mean)
resultmaxco <- apply(resultmatco, 1, max)
resultminco <- apply(resultmatco, 1, min)



resultmatk <- matrix(, nrow=101, ncol=10)

for(j in 1:10){
  r <- c()
  for(i in 1:(length(costmatk[[j]])-1)){
    r <- c(r, rep(rmsematk[[j]][i], costmatk[[j]][i+1] - costmatk[[j]][i]))
  }
  r <- c(r, rmsematk[[j]][length(rmsematk[[j]])])
  resultmatk[,j] <- r[1:101]
}

resultmatk

resultmeank <- apply(resultmatk, 1, mean)
resultmaxk <- apply(resultmatk, 1, max)
resultmink <- apply(resultmatk, 1, min)


parkalresult <- data.frame(x=rep(0:100,4), 
                                Strategy=factor(c(rep("Comprehensive searching",101), rep("Optimal point-first",101), rep("Cokriging",101), rep("MR-SUR",101)), 
                                                levels=c("Comprehensive searching","Optimal point-first","Cokriging","MR-SUR")),
                                Mean = c(resultmeanc, resultmeanc2, resultmeanco, resultmeank),
                                Max = c(resultmaxc, resultmaxc2, resultmaxco, resultmaxk),
                                Min = c(resultminc, resultminc2, resultminco, resultmink))

plotalpark <- ggplot(parkalresult, aes(x)) +                                     
  geom_line(aes(y=Mean, group=Strategy, color=Strategy, linetype=Strategy), size = 1) +
  geom_ribbon(aes(ymin=Min, ymax=Max, group=Strategy, color=Strategy, fill=Strategy), alpha=0.1, 
              color = "black", linetype = "blank")+
  theme_ipsum() +  ggtitle("Active learning for Park funcion") +
  theme(axis.title.x = element_text(size=14, margin = margin(t = 10), hjust=0.5),
        # axis.text.x = element_blank(),
        axis.title.y = element_text(size=14, margin = margin(r = 10), hjust=0.5),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  scale_y_continuous(limits = c(0, 0.58))+
  scale_x_continuous(limits = c(50, 100))+
  labs(x="Costs", y = "RMSE") #+ theme(legend.position = c(0.9,0.9))
plotalpark



parkalresult <- data.frame(x=rep(0:100,5), 
                           Strategy=factor(c(rep("ALC",101), rep("ALC+ALM",101), rep("ALM",101), rep("Cokriging",101), rep("MR-SUR",101)), 
                                           levels=c("ALC","ALC+ALM","ALM","Cokriging","MR-SUR")),
                           Mean = c(resultmeanc, resultmeanc2, resultmeanc3, resultmeanco, resultmeank),
                           Max = c(resultmaxc, resultmaxc2, resultmaxc3, resultmaxco, resultmaxk),
                           Min = c(resultminc, resultminc2, resultminc3, resultminco, resultmink))

plotalpark <- ggplot(parkalresult, aes(x)) +                                     
  geom_line(aes(y=Mean, group=Strategy, color=Strategy, linetype=Strategy), size = 1) +
  geom_ribbon(aes(ymin=Min, ymax=Max, group=Strategy, color=Strategy, fill=Strategy), alpha=0.1, 
              color = "black", linetype = "blank")+
  theme_ipsum() +  ggtitle("Active learning for Park funcion") +
  theme(axis.title.x = element_text(size=14, margin = margin(t = 10), hjust=0.5),
        # axis.text.x = element_blank(),
        axis.title.y = element_text(size=14, margin = margin(r = 10), hjust=0.5),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  scale_y_continuous(limits = c(0, 0.25))+
  scale_x_continuous(limits = c(50, 100))+
  labs(x="Costs", y = "RMSE") #+ theme(legend.position = c(0.9,0.9))

plotalpark



#######


### Numerical studies ###

model <- c(rep("RNAmf",100),rep("CoKriging",100),rep("NARGP",100))
model <- factor(model, levels=c("RNAmf", "CoKriging", "NARGP"))


nonlinear.l2error <- c(result.nonlinear.rmse[,1], result.nonlinear.rmse[,2], c(0.34022762, 0.34210415, 0.34263473, 0.34509214, 0.34552822,
                                                                               0.34640154, 0.34665543, 0.3470831 , 0.34765107, 0.3488693 ,
                                                                               0.34954502, 0.3521199 , 0.35284381, 0.35298299, 0.35411262,
                                                                               0.35513944, 0.35958212, 0.36015785, 0.36130299, 0.36154642,
                                                                               0.36222212, 0.36337812, 0.36387793, 0.36430938, 0.365791  ,
                                                                               0.3663411 , 0.36829544, 0.36864227, 0.36901338, 0.36922447,
                                                                               0.37038188, 0.37066326, 0.37076633, 0.3707668 , 0.37104885,
                                                                               0.37106391, 0.37142202, 0.37205578, 0.37229922, 0.37258915,
                                                                               0.37288867, 0.37301391, 0.37356443, 0.37848364, 0.37849521,
                                                                               0.37915099, 0.37955322, 0.37983163, 0.38087621, 0.38271045,
                                                                               0.38314046, 0.3835409 , 0.38446642, 0.38478077, 0.38487356,
                                                                               0.38642814, 0.38651514, 0.38747   , 0.38786221, 0.38793908,
                                                                               0.38848576, 0.38956627, 0.39057573, 0.39074492, 0.39074509,
                                                                               0.39206265, 0.39780143, 0.40165672, 0.40224004, 0.40319582,
                                                                               0.40441013, 0.40627892, 0.40637448, 0.40831034, 0.41006016,
                                                                               0.41131778, 0.41437364, 0.41742601, 0.42040924, 0.42262375,
                                                                               0.42381744, 0.42381759, 0.42551418, 0.42596411, 0.42739569,
                                                                               0.42739967, 0.42842417, 0.42889505, 0.43072772, 0.43418908,
                                                                               0.43465899, 0.44231468, 0.44258468, 0.44631899, 0.4543987 ,
                                                                               0.45920467, 0.47018402, 0.49496253, 0.50203073, 0.54418115))

pperd <- ggplot(data.frame(nonlinear.l2error, model), aes(x=model, y=nonlinear.l2error, fill=model)) + 
  geom_boxplot(alpha=0.5)  + 
  theme_bw() +  ggtitle("Perdikaris") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        # axis.title.y = element_text(size=18, margin = margin(r = 10), hjust=0.5),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="", y = "")+
  theme(legend.position="none")


branin.l2error <- c(result.branin.rmse[,1], result.branin.rmse[,2], c(16.51563706,  17.61260983,  19.27095676,  20.03083803,
                                                                      21.17312365,  21.49043772,  21.97356509,  22.07240125,
                                                                      22.54659889,  23.11322589,  24.05164811,  24.24957211,
                                                                      24.38184806,  25.37281651,  26.35344171,  26.50723802,
                                                                      28.36562405,  28.68204586,  28.82430544,  29.0372357 ,
                                                                      29.16104224,  29.48510606,  29.68764215,  29.71943289,
                                                                      29.98071657,  30.18067373,  30.3015576 ,  30.38552358,
                                                                      30.5151775 ,  30.98234873,  31.19926802,  31.28714229,
                                                                      31.59409527,  31.59560226,  31.91475028,  32.21009853,
                                                                      32.64194056,  32.99039166,  34.47904725,  34.65733423,
                                                                      34.77926843,  35.10888961,  36.07956805,  36.19937603,
                                                                      36.45751832,  36.58075663,  36.74966396,  37.3700588 ,
                                                                      38.05483154,  38.08946518,  38.11235063,  38.45652428,
                                                                      39.93109198,  41.07009993,  41.24307847,  42.5570682 ,
                                                                      42.91638998,  43.08970545,  43.33807062,  43.39204447,
                                                                      43.67477756,  44.16319876,  44.42081343,  44.65708265,
                                                                      45.32749562,  46.12288939,  47.11973168,  47.19152415,
                                                                      47.49545135,  47.60377949,  47.87656701,  48.00777565,
                                                                      48.08254957,  48.15646375,  48.53462565,  49.30408324,
                                                                      49.85945489,  50.32123495,  51.69290225,  51.71658022,
                                                                      51.80019439,  52.09135228,  52.57640215,  52.83387251,
                                                                      53.65493215,  53.81735274,  53.95296813,  54.1873426 ,
                                                                      54.53113248,  57.39591913,  57.51211135,  57.82645775,
                                                                      59.24087794,  59.33988191,  59.65029112,  59.76397285,
                                                                      61.93991579,  62.17425702,  62.94331531, 100.68267997))

pbra <- ggplot(data.frame(branin.l2error, model), aes(x=model, y=branin.l2error, fill=model)) + 
  geom_boxplot(alpha=0.5)  + 
  theme_bw() +  ggtitle("Branin") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        # axis.title.y = element_text(margin = margin(r = 10), hjust=0.5),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="", y = "")+
  theme(legend.position="none")


park.l2error <- c(result.park.rmse[,1], result.park.rmse[,2], c(0.01543141, 0.01883533, 0.01972081, 0.0213297 , 0.02240312,
                                                                0.02270992, 0.02487475, 0.02556971, 0.02569795, 0.026073  ,
                                                                0.02608532, 0.02721259, 0.02726363, 0.02741558, 0.02744942,
                                                                0.02746237, 0.02760419, 0.02770499, 0.02801725, 0.02918555,
                                                                0.0300006 , 0.03055055, 0.03081576, 0.03087383, 0.03109054,
                                                                0.03131905, 0.03133921, 0.03221839, 0.03278523, 0.03281587,
                                                                0.03301657, 0.03324136, 0.03391592, 0.03466853, 0.03492661,
                                                                0.03507629, 0.03530289, 0.03556985, 0.03561982, 0.03586804,
                                                                0.03596458, 0.03621814, 0.0363974 , 0.03671399, 0.03684868,
                                                                0.03735715, 0.03737333, 0.03757   , 0.03785317, 0.03818831,
                                                                0.03826123, 0.0386512 , 0.03943591, 0.04010671, 0.04013234,
                                                                0.04107583, 0.0420315 , 0.04234557, 0.04235512, 0.04248894,
                                                                0.04276381, 0.04302345, 0.04372843, 0.04385263, 0.04569816,
                                                                0.0461708 , 0.0462578 , 0.04760623, 0.04815232, 0.04850794,
                                                                0.04945962, 0.04952686, 0.05046353, 0.05080539, 0.05145118,
                                                                0.05168263, 0.05252429, 0.05318832, 0.05348071, 0.05403466,
                                                                0.05678712, 0.0577028 , 0.0577662 , 0.05849376, 0.06050538,
                                                                0.06059728, 0.06067097, 0.06178511, 0.06209458, 0.06247924,
                                                                0.06324842, 0.06385329, 0.06440754, 0.07095728, 0.0732427 ,
                                                                0.08368189, 0.09488221, 0.09797853, 0.10284276, 0.11309177))

ppark <- ggplot(data.frame(park.l2error, model), aes(x=model, y=park.l2error, fill=model)) + 
  geom_boxplot(alpha=0.5)  + 
  theme_bw() +  ggtitle("Park") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        # axis.title.y = element_text(margin = margin(r = 10), hjust=0.5),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="", y = "")+
  theme(legend.position="none")


borehole.l2error <- c(result.borehole.rmse[,1], result.borehole.rmse[,2], c(85.64001612, 85.80437339, 86.73002939, 86.95017722, 86.98712327,
                                                                            87.12769297, 87.45519789, 87.48634341, 87.56742474, 87.78279148,
                                                                            87.85621479, 87.86694683, 87.9233014 , 87.96753205, 88.17682491,
                                                                            88.21380928, 88.26297619, 88.26526763, 88.26864269, 88.29456226,
                                                                            88.38521236, 88.39347625, 88.40050458, 88.43921543, 88.50386939,
                                                                            88.53470627, 88.62088881, 88.68516991, 88.69037   , 88.72536774,
                                                                            88.74113825, 88.74784762, 88.79661657, 88.88463065, 88.88683606,
                                                                            88.93382358, 88.96772845, 89.13823852, 89.24070832, 89.3073373 ,
                                                                            89.31514742, 89.32450153, 89.34169547, 89.3947988 , 89.45621231,
                                                                            89.47074324, 89.47469893, 89.53404981, 89.64991893, 89.7310376 ,
                                                                            89.77329303, 89.77857304, 89.78571346, 89.83626268, 89.92640951,
                                                                            90.01750269, 90.03437723, 90.07276807, 90.21247281, 90.23524965,
                                                                            90.27048604, 90.27910259, 90.3586534 , 90.3712481 , 90.44051929,
                                                                            90.46377021, 90.4640873 , 90.49262182, 90.5057352 , 90.51683443,
                                                                            90.58513476, 90.69702102, 90.71305564, 90.7579166 , 90.79728355,
                                                                            90.80506742, 90.97474317, 91.01390924, 91.0202576 , 91.08809224,
                                                                            91.24277077, 91.24576713, 91.30474663, 91.31146244, 91.39138725,
                                                                            91.47915909, 91.56630331, 91.67621847, 91.69788027, 91.73715048,
                                                                            91.75316665, 91.77531024, 91.83136602, 91.94769331, 92.14544891,
                                                                            92.51136494, 92.76423062, 93.40757245, 94.59379251, 95.43260832))

pborehole <- ggplot(data.frame(borehole.l2error, model), aes(x=model, y=borehole.l2error, fill=model)) + 
  geom_boxplot(alpha=0.5)  + 
  theme_bw() +  ggtitle("Borehole") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        # axis.title.y = element_text(margin = margin(r = 10), hjust=0.5),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  scale_fill_manual(values=c(scales::hue_pal()(3)[1], scales::hue_pal()(3)[2])) +
  labs(x="", y = "")+
  scale_y_continuous(limits = c(0, 5.5))+
  theme(legend.position="none")


currin.l2error <- c(result.currin.rmse[,1], result.currin.rmse[,2], c(0.29442103, 0.3186567 , 0.32653252, 0.35884138, 0.3612266 ,
                                                                      0.3844722 , 0.39507451, 0.39712256, 0.41062758, 0.41489027,
                                                                      0.41807915, 0.4235495 , 0.42748424, 0.42926974, 0.43253144,
                                                                      0.43507348, 0.43554247, 0.43916288, 0.45281416, 0.45451204,
                                                                      0.45868008, 0.45940901, 0.46398707, 0.47300148, 0.47714258,
                                                                      0.47728713, 0.48612575, 0.49195389, 0.49343751, 0.50160201,
                                                                      0.50723251, 0.50786866, 0.5086829 , 0.51449708, 0.52858122,
                                                                      0.53043482, 0.53224675, 0.53875863, 0.54796368, 0.54964504,
                                                                      0.56915951, 0.58196077, 0.59245869, 0.60428343, 0.609985  ,
                                                                      0.63588918, 0.64744829, 0.65831741, 0.66844321, 0.6746748 ,
                                                                      0.69451237, 0.69489135, 0.70941937, 0.71535954, 0.71961617,
                                                                      0.72219468, 0.74525468, 0.74751615, 0.75195187, 0.76328598,
                                                                      0.78247198, 0.78269446, 0.78646322, 0.80000871, 0.81419144,
                                                                      0.82354549, 0.82600036, 0.83651875, 0.8392843 , 0.84059   ,
                                                                      0.87707428, 0.89894201, 0.90958787, 0.93962372, 0.9621871 ,
                                                                      0.96368175, 1.00470257, 1.0155331 , 1.04454381, 1.05469583,
                                                                      1.06000599, 1.06246045, 1.07694135, 1.08235461, 1.08749761,
                                                                      1.10841457, 1.12305717, 1.14630055, 1.15815202, 1.15942962,
                                                                      1.19638412, 1.27647976, 1.40148636, 1.44381612, 1.47485264,
                                                                      1.48467365, 1.48758505, 1.58189863, 1.70460261, 2.00790359))


pcurrin <- ggplot(data.frame(currin.l2error, model), aes(x=model, y=currin.l2error, fill=model)) + 
  geom_boxplot(alpha=0.5)  + 
  theme_bw() +  ggtitle("Currin") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        # axis.title.y = element_text(size=18, margin = margin(r = 10), hjust=0.5),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="", y = "")+
  theme(legend.position="none")


franke.l2error <- c(result.frank.rmse[,1], result.frank.rmse[,2], c(0.5511443885324832, 0.6596750505439009, 0.6849074620635947, 0.6068789900905118, 0.9611439840771514, 
                                                                    0.7710328723117196, 0.6444712743740411, 0.7431999660089059, 0.7845206726961997, 0.5759764789667287, 
                                                                    0.5346417945939718, 0.7071543256872981, 0.6118137770298608, 0.7577922060166076, 0.7878397684517019, 
                                                                    0.7666807626976353, 0.5421965681062701, 0.545833991058973, 0.7638972293390064, 0.7819956555524649, 
                                                                    0.6895415085362098, 0.6316336219290457, 0.7331878573017473, 0.5473299706700595, 0.645152935075807, 
                                                                    0.658652708707487, 0.7269947804708726, 0.7317881940229867, 0.7159739581611791, 0.7330850934213795, 
                                                                    0.7629645399398459, 0.6688970841603887, 0.8547300123428879, 0.6206954813612146, 0.6089246770888753, 
                                                                    0.5704676662722323, 0.6386987795868743, 0.6394497492724655, 0.6203402216163288, 0.86373560213365, 
                                                                    0.6960460619649649, 0.7757263310792275, 0.5804835363241708, 0.6694584145628527, 0.6164537199395331, 
                                                                    0.6760912580346445, 0.6829915842859478, 0.6163662085979954, 0.6472182543008344, 0.7451274508814558, 
                                                                    0.553385608420307, 0.5985391204385572, 0.6725486911652034, 0.6915490689913442, 0.6725604208633503, 
                                                                    0.5500150389853841, 0.7048326743971403, 0.6316724808981606, 0.6366837699955432, 0.5592483178288905, 
                                                                    0.6413031581504399, 0.6413031581504399, 0.7161672309296464, 0.5684778574701264, 0.7972733318813475, 
                                                                    0.6592057639291572, 0.5411544077755798, 0.5616109509795045, 0.5909782876115647, 0.7043338395691602, 
                                                                    0.554034980033292, 0.5692825068446294, 0.6426497978864342, 0.5667769778770291, 0.6797704289124331, 
                                                                    0.5266794608367036, 0.7126113752623242, 0.7564715599685653, 0.630750626680051, 0.5950322307473801, 
                                                                    0.7824573608363246, 0.6438479624136744, 0.8287452873603723, 0.630105545550783, 0.5581817122963557, 
                                                                    0.72566043402063, 0.7098667663519718, 0.6224230110470849, 0.6830183107679724, 0.5528958709325407, 
                                                                    0.49071890796148643, 0.6587663151066209, 0.5562792642491277, 0.7542139774838275, 0.6017206228079495, 
                                                                    0.6669038394139393, 0.7459977160419697, 0.6733708865506871, 0.6501331492665351, 0.7206131698068257))

pfranke <- ggplot(data.frame(franke.l2error, model), aes(x=model, y=franke.l2error, fill=model)) + 
  geom_boxplot(alpha=0.5)  + 
  theme_bw() +  ggtitle("Franke") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        # axis.title.y = element_text(size=18, margin = margin(r = 10), hjust=0.5),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="", y = "")+
  theme(legend.position="none")


nonlinear.meancrps <- c(result.nonlinear.meancrps[,1], result.nonlinear.meancrps[,2], c(0.19601248, 0.19738676, 0.19935524, 0.19936653, 0.20161979,
                                                                                        0.20162493, 0.20333192, 0.20357583, 0.20409573, 0.20424958,
                                                                                        0.20465862, 0.20684365, 0.20722672, 0.20761565, 0.20846041,
                                                                                        0.20983259, 0.21087637, 0.21155828, 0.21189774, 0.21199481,
                                                                                        0.21221467, 0.21250487, 0.21254519, 0.2125626 , 0.21264963,
                                                                                        0.21265009, 0.21268746, 0.21278098, 0.21292538, 0.21302776,
                                                                                        0.21364283, 0.21425327, 0.21429693, 0.21436857, 0.21456071,
                                                                                        0.21468787, 0.21531511, 0.21537192, 0.21556793, 0.21610559,
                                                                                        0.21616138, 0.21639508, 0.21668155, 0.21685781, 0.2176249 ,
                                                                                        0.21920487, 0.21950093, 0.21967   , 0.21967008, 0.22000444,
                                                                                        0.2214316 , 0.22251253, 0.22274926, 0.22527235, 0.22570208,
                                                                                        0.22635634, 0.22637389, 0.22663662, 0.22748233, 0.22788377,
                                                                                        0.22817972, 0.2281916 , 0.22936076, 0.23056528, 0.23141645,
                                                                                        0.2327797 , 0.23307097, 0.23513764, 0.2355944 , 0.23648214,
                                                                                        0.23648222, 0.23787876, 0.23852823, 0.23853046, 0.23913112,
                                                                                        0.23946125, 0.24158279, 0.2440589 , 0.24746223, 0.2483642 ,
                                                                                        0.24880676, 0.24899884, 0.25231504, 0.25331067, 0.25364598,
                                                                                        0.2598985 , 0.2620474 , 0.26451653, 0.27220879, 0.27291315,
                                                                                        0.27703373, 0.30319667, 0.32853959, 0.33597343, 0.37081451,
                                                                                        0.38377157, 0.44192297, 0.55535666, 0.65073497, 0.69197052))

pperd2 <- ggplot(data.frame(nonlinear.meancrps, model), aes(x=model, y=nonlinear.meancrps, fill=model)) + 
  geom_boxplot(alpha=0.5)  + 
  theme_bw() +  ggtitle("Perdikaris") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        # axis.title.y = element_text(size=18, margin = margin(r = 10), hjust=0.5),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="", y = "")+
  theme(legend.position="none")


branin.meancrps <- c(result.branin.meancrps[,1], result.branin.meancrps[,2], c(7.93069595,  9.43343044,  9.54151377,  9.72126493, 10.39037927,
                                                                               10.49107708, 10.51889256, 10.78811531, 10.83320515, 11.01893776,
                                                                               11.04393738, 11.09730405, 11.67525392, 11.96963657, 11.97959774,
                                                                               12.07384909, 12.69938011, 13.23650216, 13.36906042, 13.81575332,
                                                                               14.46785737, 14.56691336, 14.62767005, 14.69378076, 14.77545607,
                                                                               14.79852786, 14.86748539, 15.07782292, 15.20813245, 15.39326357,
                                                                               15.51186074, 15.57391663, 15.68483006, 15.84166467, 15.93522859,
                                                                               16.41943072, 16.47444469, 16.65862628, 16.8445593 , 17.10562516,
                                                                               17.16838018, 17.67307722, 18.03696358, 18.22482711, 18.27345742,
                                                                               18.62117588, 18.67404549, 18.98408725, 19.26929307, 19.52907175,
                                                                               19.7236741 , 19.89496869, 20.74817045, 20.90615823, 20.98521786,
                                                                               21.18327277, 21.24209368, 21.27275237, 21.36112172, 21.36218409,
                                                                               21.92558969, 22.52076681, 22.99033822, 23.01103413, 23.06385702,
                                                                               23.17814015, 23.90019583, 23.9649098 , 24.08005458, 24.35688269,
                                                                               25.2657552 , 25.70082389, 25.8761768 , 25.93351156, 26.4556661 ,
                                                                               26.58978344, 26.74993292, 28.075578  , 28.21124157, 28.30089178,
                                                                               28.52870706, 28.63183036, 28.80295901, 28.98688465, 29.23877078,
                                                                               30.51799923, 30.54782418, 30.99171635, 31.10885636, 31.25312795,
                                                                               31.70788774, 32.2403172 , 32.70030603, 32.95489386, 33.11272083,
                                                                               33.14405864, 36.87733569, 38.90765285, 39.42471151, 70.28433906))

pbra2 <- ggplot(data.frame(branin.meancrps, model), aes(x=model, y=branin.meancrps, fill=model)) + 
  geom_boxplot(alpha=0.5)  + 
  theme_bw() +  ggtitle("Branin") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        # axis.title.y = element_text(size=18, margin = margin(r = 10), hjust=0.5),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="", y = "")+
  theme(legend.position="none")


park.meancrps <- c(result.park.meancrps[,1], result.park.meancrps[,2], c(0.0077333 , 0.00940334, 0.01037537, 0.01071631, 0.01215122,
                                                                         0.01246912, 0.01258697, 0.01325076, 0.01348918, 0.01350012,
                                                                         0.01378779, 0.01380279, 0.01406547, 0.01411974, 0.01423231,
                                                                         0.01428311, 0.01453217, 0.01458184, 0.01484368, 0.01486118,
                                                                         0.01524672, 0.015531  , 0.01562671, 0.01575025, 0.01630846,
                                                                         0.01656625, 0.01669695, 0.01685203, 0.01690171, 0.01693753,
                                                                         0.01722525, 0.01724074, 0.01725906, 0.01732064, 0.01736938,
                                                                         0.01759915, 0.01767862, 0.01793282, 0.01798829, 0.01813647,
                                                                         0.01831486, 0.01885326, 0.01892137, 0.01909389, 0.01918683,
                                                                         0.01923887, 0.01953403, 0.01984809, 0.01990239, 0.0199531 ,
                                                                         0.02001505, 0.02001579, 0.02037022, 0.02047832, 0.02065944,
                                                                         0.02169165, 0.0221693 , 0.02236053, 0.02264802, 0.02279485,
                                                                         0.0231223 , 0.02391821, 0.02394734, 0.02396089, 0.02397694,
                                                                         0.02407729, 0.02412081, 0.02419562, 0.0246353 , 0.02497756,
                                                                         0.02525212, 0.02538408, 0.02634231, 0.02648956, 0.02682843,
                                                                         0.02746085, 0.0278626 , 0.0279504 , 0.02896078, 0.03004806,
                                                                         0.03005748, 0.03046396, 0.03076343, 0.03118792, 0.03157026,
                                                                         0.03158957, 0.03168303, 0.03340951, 0.03344207, 0.0343543 ,
                                                                         0.03574715, 0.03675782, 0.036798  , 0.03707381, 0.03962727,
                                                                         0.04007912, 0.04251469, 0.04393943, 0.04686675, 0.04702597))

ppark2 <- ggplot(data.frame(park.meancrps, model), aes(x=model, y=park.meancrps, fill=model)) + 
  geom_boxplot(alpha=0.5)  + 
  theme_bw() +  ggtitle("Park") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        # axis.title.y = element_text(size=18, margin = margin(r = 10), hjust=0.5),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="", y = "")+
  theme(legend.position="none")


borehole.meancrps <- c(result.borehole.meancrps[,1], result.borehole.meancrps[,2], c(49.29462471, 49.2980267 , 50.09820259, 50.10611749, 50.13905544,
                                                                                     50.2642037 , 50.31809112, 50.45530659, 50.47234833, 50.56785915,
                                                                                     50.63670308, 50.66977955, 50.71428898, 50.73148708, 50.82184035,
                                                                                     50.88616105, 50.88804366, 50.89153971, 50.92732821, 50.93055938,
                                                                                     50.95578799, 50.96262694, 50.9823663 , 51.0253229 , 51.05430568,
                                                                                     51.07780624, 51.0780425 , 51.10396648, 51.1840877 , 51.20991604,
                                                                                     51.28421588, 51.31527123, 51.33451128, 51.34107754, 51.37791952,
                                                                                     51.37958484, 51.40765674, 51.44229555, 51.46384142, 51.59682893,
                                                                                     51.61899525, 51.6746597 , 51.69808218, 51.7027772 , 51.70894838,
                                                                                     51.73378427, 51.76525051, 51.78308338, 51.78315322, 51.80762076,
                                                                                     51.82194404, 51.83817118, 51.84682794, 51.86668421, 51.88695548,
                                                                                     51.89809504, 51.91261588, 51.98959469, 52.00891176, 52.02293739,
                                                                                     52.05633049, 52.0896567 , 52.10920556, 52.12855853, 52.13015747,
                                                                                     52.13356444, 52.2188267 , 52.28982923, 52.2953657 , 52.37094534,
                                                                                     52.38177441, 52.39335955, 52.39580686, 52.47252351, 52.48761522,
                                                                                     52.51523912, 52.52257224, 52.69457141, 52.69949172, 52.77416841,
                                                                                     52.7806585 , 52.7959367 , 52.80766792, 52.89438096, 52.90452621,
                                                                                     52.94153982, 52.94372474, 52.98543036, 52.99177638, 53.05195537,
                                                                                     53.0609831 , 53.08978817, 53.15678613, 53.16424325, 53.38292602,
                                                                                     53.88583814, 53.99809494, 54.41812143, 55.08558337, 56.51043403))

pborehole2 <- ggplot(data.frame(borehole.meancrps, model), aes(x=model, y=borehole.meancrps, fill=model)) + 
  geom_boxplot(alpha=0.5)  + 
  theme_bw() +  ggtitle("Borehole") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        # axis.title.y = element_text(size=18, margin = margin(r = 10), hjust=0.5),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  scale_fill_manual(values=c(scales::hue_pal()(3)[1], scales::hue_pal()(3)[2])) +
  labs(x="", y = "")+
  scale_y_continuous(limits = c(0, 1.5))+
  theme(legend.position="none")



currin.meancrps <- c(result.currin.meancrps[,1], result.currin.meancrps[,2], c(0.14984793, 0.15505572, 0.15752919, 0.1698199 , 0.17900319,
                                                                               0.18479192, 0.18752176, 0.19301432, 0.19358862, 0.19858049,
                                                                               0.19965511, 0.20126463, 0.206678  , 0.2098936 , 0.21043123,
                                                                               0.21064509, 0.21116791, 0.21476699, 0.21666309, 0.21724425,
                                                                               0.21825302, 0.21901363, 0.22262649, 0.2237749 , 0.22923687,
                                                                               0.22969995, 0.23134981, 0.23409576, 0.23670036, 0.23862148,
                                                                               0.2401337 , 0.24342123, 0.24376436, 0.24652922, 0.24822088,
                                                                               0.25209849, 0.25473364, 0.25586932, 0.25760047, 0.25864618,
                                                                               0.26220888, 0.26239111, 0.26497999, 0.26526329, 0.26629332,
                                                                               0.27691318, 0.27783033, 0.28691924, 0.28803394, 0.28989036,
                                                                               0.29701956, 0.29913689, 0.29994336, 0.29995128, 0.30109846,
                                                                               0.30355179, 0.30439197, 0.31807817, 0.31906798, 0.33022416,
                                                                               0.33680343, 0.33803476, 0.34622141, 0.34853128, 0.35212875,
                                                                               0.36136308, 0.37145629, 0.3787886 , 0.38130153, 0.38332794,
                                                                               0.38349665, 0.38485145, 0.38525533, 0.39298484, 0.39986938,
                                                                               0.40135883, 0.40558233, 0.40605587, 0.40960297, 0.41125834,
                                                                               0.42218459, 0.43216168, 0.43272242, 0.43966233, 0.47053592,
                                                                               0.47574082, 0.48443594, 0.48444645, 0.49302779, 0.51296516,
                                                                               0.53790712, 0.55106439, 0.56819775, 0.59170665, 0.61832829,
                                                                               0.64467375, 0.67897264, 0.72814085, 0.76438428, 0.80044238))

pcurrin2 <- ggplot(data.frame(currin.meancrps, model), aes(x=model, y=currin.meancrps, fill=model)) + 
  geom_boxplot(alpha=0.5)  + 
  theme_bw() +  ggtitle("Currin") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        # axis.title.y = element_text(size=18, margin = margin(r = 10), hjust=0.5),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="", y = "")+
  theme(legend.position="none")


franke.meancrps <- c(result.frank.meancrps[,1], result.frank.meancrps[,2], c(0.28897173135762844, 0.38001987647033136, 0.3901871913076911, 0.3357667552290378, 0.6264033853526362, 
                                                                             0.4529880578354316, 0.37249206921637074, 0.4369112078310833, 0.46669746029992437, 0.30873200355816544, 
                                                                             0.271851155047444, 0.4117122740247307, 0.33595121771805747, 0.4643573793720561, 0.47851924121490436, 
                                                                             0.4378370647683306, 0.2930681478049687, 0.2760983591543816, 0.44896803838417826, 0.43619765767098767, 
                                                                             0.40120302561334636, 0.35660728459283403, 0.44093907597107945, 0.2948376108011204, 0.35395133930696315, 
                                                                             0.3671513201565945, 0.4215738457264353, 0.4272297102921997, 0.41812540499164264, 0.42769246404724426, 
                                                                             0.4794921021578388, 0.35065196224582773, 0.5286052847680187, 0.33986295354986046, 0.31850489724131875, 
                                                                             0.3023930642857322, 0.3453335292423643, 0.3494420314138799, 0.34152032294634777, 0.5693658970144243, 
                                                                             0.40025192749570815, 0.4576638706000773, 0.3068501793036992, 0.3816272007445675, 0.3269210695917484, 
                                                                             0.3792022791785919, 0.3945377293233353, 0.3466095127194851, 0.3579534395119872, 0.43948083088108875, 
                                                                             0.28396416946041725, 0.3282977689523598, 0.3793415474288098, 0.40008890731927227, 0.3902492069684566, 
                                                                             0.2859709903516653, 0.4124655131866682, 0.3399022895697458, 0.3479461745272443, 0.28775951779433645, 
                                                                             0.3592086329002883, 0.3592086329002883, 0.4227412534088286, 0.30679821168632987, 0.5090249013225981, 
                                                                             0.36806014901993556, 0.28371676325676054, 0.28360807883858463, 0.3132943576721815, 0.4019755267308435, 
                                                                             0.2841378338294835, 0.3018195349438574, 0.3526153178360991, 0.2971083951586694, 0.3832082712950782, 
                                                                             0.2808789198484714, 0.4282899995226224, 0.4438494715903681, 0.34094704751487814, 0.31135270388383013, 
                                                                             0.45026819296648113, 0.36120073045966944, 0.5036575037944042, 0.3392850145050295, 0.30108261551912213, 
                                                                             0.4153547538434583, 0.4168722046483455, 0.34997293150970393, 0.39011951821543445, 0.29832171560471216, 
                                                                             0.25423372078041456, 0.36628755380871325, 0.2936671677835252, 0.4555573892485559, 0.32913400796557085, 
                                                                             0.3857929357877079, 0.42711877998828196, 0.35034509792775165, 0.34785563483732046, 0.4205077288360665))

pfranke2 <- ggplot(data.frame(franke.meancrps, model), aes(x=model, y=franke.meancrps, fill=model)) + 
  geom_boxplot(alpha=0.5)  + 
  theme_bw() +  ggtitle("Franke") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        # axis.title.y = element_text(size=18, margin = margin(r = 10), hjust=0.5),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_blank(),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="", y = "")+
  theme(legend.position="none")




ggarrange(pperd, pbra, ppark, pborehole, pcurrin, pfranke, ncol=3, nrow=2, common.legend = TRUE, legend="bottom")
# annotate_figure(ggarrange(pperd, pbra, ppark, pborehole, pcurrin, pfranke, ncol=3, nrow=2, common.legend = TRUE, legend="bottom"), 
#                 bottom = text_grob("RMSE", color = "black", face = "bold", size = 20))
ggarrange(pperd2, pbra2, ppark2, pborehole2, pcurrin2, pfranke2, ncol=3, nrow=2, common.legend = TRUE, legend="bottom")
# ggarrange(p10, p11, p12, pp10, pp11, pp12, ncol=3, nrow=2, common.legend = TRUE, legend="bottom")



perd.comptime <- c(result.nonlinear.comptime[,1], result.nonlinear.comptime[,2], c(12.638431072235107, 14.603119134902954, 9.038669109344482, 12.611816167831421, 8.181637287139893, 
                                                                                   8.41358995437622, 11.093848943710327, 12.64332103729248, 10.260518789291382, 8.836445093154907, 
                                                                                   7.823691129684448, 9.36992883682251, 7.161639928817749, 8.045645236968994, 7.5129311084747314, 
                                                                                   10.962071180343628, 9.90575623512268, 13.129024028778076, 9.923938035964966, 10.478762149810791, 
                                                                                   7.7480552196502686, 5.4767560958862305, 5.874589920043945, 7.085026025772095, 9.507341861724854, 
                                                                                   7.090926647186279, 6.58402419090271, 11.1159029006958, 7.510082721710205, 9.444695949554443, 
                                                                                   11.468750715255737, 7.939141035079956, 9.712474822998047, 7.494318008422852, 6.00114893913269, 
                                                                                   7.759785890579224, 9.164674997329712, 8.371340036392212, 12.275887966156006, 8.428921937942505, 
                                                                                   8.46586298942566, 9.01035213470459, 11.580531120300293, 8.277020931243896, 13.351361751556396, 
                                                                                   11.552864074707031, 10.682437896728516, 10.621850967407227, 7.222555875778198, 9.543353080749512, 
                                                                                   10.435199975967407, 10.084581136703491, 11.306042909622192, 12.817527055740356, 8.33662486076355, 
                                                                                   10.45873498916626, 11.3093900680542, 8.099806070327759, 10.421024084091187, 14.024849891662598, 
                                                                                   8.016197681427002, 11.732913970947266, 11.626074075698853, 12.83341908454895, 12.183241844177246, 
                                                                                   11.139633893966675, 9.099885702133179, 9.223941802978516, 11.819476842880249, 9.258319854736328, 
                                                                                   10.784610986709595, 7.936137914657593, 7.6837029457092285, 7.857706785202026, 11.53866195678711, 
                                                                                   10.071309089660645, 9.193354845046997, 8.42283320426941, 14.498060941696167, 15.400633096694946, 
                                                                                   12.906683921813965, 12.804841995239258, 9.027191638946533, 12.624915838241577, 14.10118818283081, 
                                                                                   8.544509887695312, 10.486035108566284, 13.761353015899658, 9.055784940719604, 11.697642803192139, 
                                                                                   10.145763158798218, 12.498543739318848, 11.528660297393799, 7.321206092834473, 13.920256853103638, 
                                                                                   11.943364143371582, 8.295090913772583, 11.669752836227417, 10.307748317718506, 14.031654119491577))

bra.comptime <- c(result.branin.comptime[,1], result.branin.comptime[,2], c(139.29871487617493, 162.51807618141174, 146.56250405311584, 140.66339707374573, 158.13290071487427, 
                                                                            147.68331122398376, 155.9636242389679, 145.16405200958252, 150.1841320991516, 147.00824093818665, 
                                                                            160.47584414482117, 150.96788907051086, 140.0498218536377, 147.64123916625977, 170.92029285430908, 
                                                                            147.69771003723145, 142.63090920448303, 147.36768412590027, 158.81940603256226, 155.16283583641052, 
                                                                            146.09358882904053, 156.11737513542175, 141.87225008010864, 147.26390290260315, 149.22877407073975, 
                                                                            143.82359290122986, 139.06521797180176, 151.26320695877075, 159.24053812026978, 159.7500557899475, 
                                                                            146.2730691432953, 137.30180287361145, 146.83052515983582, 138.1788101196289, 149.07607412338257, 
                                                                            152.51751899719238, 158.44029903411865, 149.99088191986084, 153.45190000534058, 146.9187879562378, 
                                                                            135.79412007331848, 137.6956911087036, 131.6539227962494, 121.37452006340027, 130.94265818595886, 
                                                                            133.0022599697113, 131.50036096572876, 130.9718210697174, 133.38457798957825, 112.36569213867188, 
                                                                            132.10127425193787, 115.22399830818176, 139.21120619773865, 130.9202868938446, 131.9804151058197, 
                                                                            141.5452208518982, 156.39409589767456, 160.71244382858276, 158.8587532043457, 157.62313866615295, 
                                                                            151.53280806541443, 158.65682983398438, 145.93891191482544, 152.2540988922119, 128.97311997413635, 
                                                                            140.41571307182312, 124.02498006820679, 122.42486214637756, 130.61541986465454, 124.27168083190918, 
                                                                            129.37533283233643, 134.17281985282898, 130.25748801231384, 131.57602286338806, 129.56485891342163, 
                                                                            113.14144897460938, 109.5228500366211, 100.79314517974854, 100.6189169883728, 97.40959692001343, 
                                                                            101.4079077243805, 96.97514390945435, 103.08406090736389, 102.71305871009827, 106.89587378501892, 
                                                                            106.34330081939697, 99.79533576965332, 97.57815408706665, 106.1473662853241, 94.20387697219849, 
                                                                            102.51707100868225, 101.41523885726929, 93.11203002929688, 93.22160816192627, 97.22728395462036, 
                                                                            89.63261795043945, 95.05528521537781, 90.62356877326965, 93.9550347328186, 95.87637901306152))

park.comptime <- c(result.park.comptime[,1], result.park.comptime[,2], c(37.59697079658508, 38.284141063690186, 48.01426720619202, 49.72091484069824, 40.468273878097534, 
                                                                         72.36816811561584, 53.14524006843567, 57.79141402244568, 57.793092012405396, 67.21538805961609, 
                                                                         63.308337926864624, 70.35309505462646, 55.83506202697754, 75.59888887405396, 56.60691595077515, 
                                                                         65.50443291664124, 72.42974328994751, 52.16663980484009, 55.941831827163696, 51.217947006225586, 
                                                                         61.65242314338684, 59.505435943603516, 69.45346474647522, 48.75571894645691, 49.88899302482605, 
                                                                         46.45164394378662, 47.82719326019287, 72.98168015480042, 52.37946820259094, 72.53865194320679, 
                                                                         66.71575880050659, 47.082996129989624, 54.51527810096741, 56.45328712463379, 54.735270977020264, 
                                                                         45.335487842559814, 64.70993304252625, 74.58084201812744, 60.62234807014465, 85.19076013565063, 
                                                                         40.25997591018677, 56.78615617752075, 53.69546890258789, 81.31187105178833, 66.47054076194763, 
                                                                         64.86935496330261, 58.13776111602783, 49.64909291267395, 57.85577392578125, 68.54154515266418, 
                                                                         39.86267280578613, 41.29041910171509, 69.75921583175659, 53.208172082901, 62.091586112976074, 
                                                                         64.73214483261108, 64.6875729560852, 53.90369987487793, 65.9427239894867, 54.63913702964783, 
                                                                         62.38015675544739, 63.685375928878784, 57.33064889907837, 51.27383279800415, 46.67505383491516, 
                                                                         72.6160261631012, 51.69126915931702, 57.06262683868408, 65.0576171875, 87.18740820884705, 
                                                                         58.21577310562134, 55.46340012550354, 42.650124073028564, 56.20947003364563, 47.498610973358154, 
                                                                         46.2004599571228, 58.55592608451843, 73.1507019996643, 44.24789786338806, 82.06606888771057, 
                                                                         60.84273290634155, 56.2906608581543, 51.66313624382019, 51.88989973068237, 44.84682822227478, 
                                                                         60.61213994026184, 53.8550820350647, 63.2651309967041, 61.434524059295654, 62.52141308784485, 
                                                                         75.24477577209473, 76.4413161277771, 53.047552824020386, 45.835747957229614, 71.38743782043457, 
                                                                         56.08518886566162, 36.63569402694702, 61.86357665061951, 61.033599853515625, 56.07898783683777))

borehole.comptime <- c(result.borehole.comptime[,1], result.borehole.comptime[,2], rep(1000,100))

currin.comptime <- c(result.currin.comptime[,1], result.currin.comptime[,2], c(38.46457505226135, 45.29613900184631, 41.67369318008423, 13.852711915969849, 66.1725401878357, 
                                                                               48.59166121482849, 42.97861981391907, 53.27743482589722, 42.89934515953064, 42.774128913879395, 
                                                                               46.519551038742065, 36.85790181159973, 41.81056785583496, 40.66720390319824, 41.6387300491333, 
                                                                               35.29694604873657, 45.70790982246399, 45.42031407356262, 56.07017803192139, 60.64559721946716, 
                                                                               26.092360973358154, 53.423269271850586, 38.765719175338745, 29.527088165283203, 38.00006413459778, 
                                                                               53.50367999076843, 45.0527720451355, 73.1283929347992, 29.322190046310425, 47.27030396461487, 
                                                                               62.09379696846008, 24.999520301818848, 49.110761880874634, 18.03965473175049, 58.86049675941467, 
                                                                               66.25800490379333, 37.40067386627197, 33.23986196517944, 44.91863489151001, 67.02728796005249, 
                                                                               63.02077388763428, 33.2865469455719, 41.235963106155396, 50.644490003585815, 44.38124108314514, 
                                                                               52.33662271499634, 25.506235122680664, 41.71459484100342, 71.1365180015564, 51.028106927871704, 
                                                                               51.71041417121887, 36.98497986793518, 34.90028500556946, 37.8846960067749, 52.581767082214355, 
                                                                               48.512667179107666, 36.47699785232544, 42.46441912651062, 49.92925405502319, 49.961191177368164, 
                                                                               35.59922981262207, 44.32139492034912, 65.25653123855591, 51.57457995414734, 31.872881174087524, 
                                                                               36.46962308883667, 44.602237939834595, 30.923030853271484, 30.89026689529419, 40.68335199356079, 
                                                                               55.80264902114868, 29.661072254180908, 47.33366298675537, 34.87639594078064, 46.38042712211609, 
                                                                               42.20359015464783, 44.61259698867798, 26.800437927246094, 66.36775994300842, 32.11775302886963, 
                                                                               42.078288316726685, 27.806385040283203, 56.80334281921387, 44.93955183029175, 36.94848990440369, 
                                                                               38.07759118080139, 31.344873905181885, 60.28279399871826, 77.21438503265381, 16.452447175979614, 
                                                                               41.06525111198425, 42.43532395362854, 53.989501953125, 49.655259132385254, 49.13196682929993, 
                                                                               34.81497597694397, 55.88296604156494, 46.69751811027527, 35.98944282531738, 32.41061568260193))

franke.comptime <- c(result.frank.comptime[,1], result.frank.comptime[,2], c(131.58816719055176, 193.98071193695068, 175.7153582572937, 196.55459213256836, 189.7441108226776, 
                                                                             180.15572500228882, 190.9534330368042, 167.68277502059937, 178.62781882286072, 159.25829005241394, 
                                                                             172.88579297065735, 212.43094873428345, 175.0363998413086, 182.27984929084778, 176.92628717422485, 
                                                                             183.5248727798462, 184.2406063079834, 165.13108801841736, 213.01585602760315, 196.10092306137085, 
                                                                             171.23573780059814, 184.95077896118164, 187.24046802520752, 184.54044795036316, 151.77827095985413, 
                                                                             167.4719729423523, 168.96429085731506, 170.85820078849792, 176.580904006958, 207.66964411735535, 
                                                                             156.30769085884094, 181.0253129005432, 169.366192817688, 142.02578592300415, 151.27464389801025, 
                                                                             139.34179711341858, 138.86328291893005, 147.94982600212097, 175.71044993400574, 155.05779314041138, 
                                                                             156.02310585975647, 175.8110179901123, 198.23550510406494, 183.66553497314453, 156.96297192573547, 
                                                                             170.78364396095276, 189.10905170440674, 184.97563695907593, 173.18184113502502, 187.47554183006287, 
                                                                             200.18465995788574, 182.57795095443726, 133.4921441078186, 140.89774918556213, 150.2748682498932, 
                                                                             130.7517387866974, 151.8335099220276, 130.1663429737091, 140.09835171699524, 137.67639589309692, 
                                                                             176.2342598438263, 140.63138008117676, 114.55709886550903, 111.80307674407959, 131.57559084892273, 
                                                                             98.31709218025208, 99.95511603355408, 107.54872965812683, 105.12071967124939, 127.22969818115234, 
                                                                             105.37112188339233, 98.97530913352966, 96.53790712356567, 99.96590375900269, 122.99387812614441, 
                                                                             105.3841028213501, 113.64740800857544, 111.2693121433258, 98.93330216407776, 99.08079481124878, 
                                                                             113.4723129272461, 113.66482305526733, 122.41466021537781, 102.70020604133606, 113.67346811294556, 
                                                                             120.50010991096497, 108.61943507194519, 111.11429405212402, 110.45885372161865, 100.26177430152893, 
                                                                             101.47745490074158, 102.3727560043335, 101.48106694221497, 118.30065703392029, 109.76288890838623, 
                                                                             98.3718409538269, 127.53809475898743, 105.30002117156982, 99.86162495613098, 98.55333828926086))


model <- c(rep("RNAmf",100),rep("CoKriging",100),rep("NARGP",100))
model <- factor(model, levels=c("RNAmf", "CoKriging", "NARGP"))



df.comptime <- data.frame(comptime=c(perd.comptime, bra.comptime, park.comptime, borehole.comptime, currin.comptime, franke.comptime), 
                          model=rep(model,6), 
                          Function=c(rep("Perdikaris",300),rep("Branin",300),rep("Park",300),rep("Borehole",300),rep("Currin",300),rep("Franke",300)), 
                          xx=c(rep("Perdikaris-RNAmf",100),rep("Perdikaris-CoKriging",100),rep("Perdikaris-NARGP",100),
                               rep("Branin-RNAmf",100),rep("Branin-CoKriging",100),rep("Branin-NARGP",100),
                               rep("Park-RNAmf",100), rep("Park-CoKriging",100), rep("Park-NARGP",100),
                               rep("Borehole-RNAmf",100),rep("Borehole-CoKriging",100),rep("Borehole-NARGP",100),
                               rep("Currin-RNAmf",100),rep("Currin-CoKriging",100),rep("Currin-NARGP",100),
                               rep("Franke-RNAmf",100),rep("Franke-CoKriging",100),rep("Franke-NARGP",100))
)

df.comptime$Function <- factor(df.comptime$Function , levels=c("Perdikaris", "Branin", "Park", "Borehole", "Currin", "Franke"))

ggplot(df.comptime, 
       aes(x=Function, y=comptime, fill=model, color=model)) + 
  geom_boxplot(alpha=0.5)  +
  theme_ipsum() +  ggtitle("") +
  theme(
    # axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        # axis.text.x = element_blank(),
        # axis.title.y = element_blank(),
        axis.title.y = element_text(size=16, margin = margin(r = 10), hjust=0.5),
        strip.text.x = element_blank(),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="", y = "Computational time (seconds)")+
  theme(legend.position="bottom",
        panel.spacing = unit(0, "lines")#,
        # panel.border = element_rect(fill = NA, color = "white")
  ) +
  geom_vline(xintercept = c(0.4,6.6), linetype="dotted")+
  scale_y_continuous(
    limits = c(0,214)) +
  facet_grid(. ~ model)





### Blade comparison ### load blade mean2 rdata

model <- c(rep("RNAmf",50),rep("CoKriging",50),rep("NARGP",50))
model <- factor(model, levels=c("RNAmf", "CoKriging", "NARGP"))


blade.l2error <- c(result.blade.rmse[,1],result.blade.rmse[,2],c(0.08296562284056933, 0.26944568322745277, 0.08139932349505961, 0.06591543174151877, 0.07539813412693266, 
                                                                 0.054232110710972156, 0.08305997653664889, 0.06354733784199235, 0.05857402773073935, 0.1358282065717591, 
                                                                 0.05184968080311192, 0.060033594336085724, 0.08171029623631187, 0.17766831135547617, 0.05291498799731878, 
                                                                 0.06584816603092461, 0.03009986302911407, 0.04810813255207206, 0.06714902793781717, 0.10171316513508435, 
                                                                 0.06671993366295795, 0.050836904522805114, 0.12107033428616551, 0.06044577038368387, 0.0677842654704768, 
                                                                 0.0712270546508961, 0.08555776795597886, 0.06678370359631516, 0.07219886395976191, 0.055628482275447574, 
                                                                 0.08072438995677511, 0.055625115918103755, 0.07098293549593825, 0.06434466047336887, 0.10888285474073596, 
                                                                 0.07422725981139436, 0.0633772920358081, 0.05796373984467013, 0.05492804850918673, 0.08209933537618808,
                                                                 0.05932549535668908, 0.05262276438485215, 0.05718896275464744, 0.06434278236699577, 0.0738492087860961, 
                                                                 0.06421319849590207, 0.11957146209403301, 0.12055272961881312, 0.07923707571466478, 0.0432476866983063))* (sqrt(sum((y.test)^2))) / sqrt(100) # 224.5119 / sqrt(100)

p13 <- ggplot(data.frame(blade.l2error, model), aes(x=model, y=blade.l2error, fill=model)) + 
  geom_boxplot(alpha=0.5)  + 
  theme_ipsum() +  ggtitle("Blade data") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.title.y = element_text(size=18, margin = margin(r = 10), hjust=0.5),
        axis.text.x = element_blank(),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="", y = "RMSE")+
  theme(legend.position="none")


apply(matrix(blade.l2error, nrow=50), 2, mean) # 1.825184 3.213079 1.739078

# blade.meanscore <- c(result.blade.meanscore[,1],result.blade.meanscore[,2],c(-1.226594881411466, -5.918071411935286, -1.9558134642277685, -1.1085329528957575, -1.488964865587706,
#                                                                              -1.0080137923033425, -1.4315922877597598, -1.7087502834053412, -1.2037667689509974, -1.972561304609182, 
#                                                                              -1.4881905697693818, -1.6943623807970218, -14.979284910707147, -2.1177263505148507, -2.569032939900389,
#                                                                              -1.4842397062145218, -0.8799405048802597, -1.090908447988584, -1.4026922180858028, -7.702982712178904,
#                                                                              -1.6542554315350002, -1.2317690422106842, -1.774163331111413, -1.3164878700815088, -2.016471303590664,
#                                                                              -1.7111590264163186, -1.7097943540423484, -1.815623748849739, -1.8387515849737113, -1.3239455410075098, 
#                                                                              -2.9060598636177826,  -1.409988196868963, -1.7499975034338784, -1.2701768252674852,-4.132453999219285, 
#                                                                              -1.6906501203035058, -1.3586512654534149, -1.4651167782278274, -1.1106552422294056, -2.0169536276943574, 
#                                                                              -1.4044022126524454, -1.1321504460605327, -1.46435129580732, -1.7986160831840154, -1.0989744158653132, 
#                                                                              -1.501738937874928, -1.5566837480191227, -2.2544837541746277, -1.6985180936658764, -1.6379162095425286))
# 
# p14 <- ggplot(data.frame(blade.meanscore, model), aes(x=model, y=blade.meanscore, fill=model)) + 
#   geom_boxplot(alpha=0.5)  + 
#   theme_ipsum() +  ggtitle("Blade data") +
#   theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
#         axis.title.y = element_text(size=18, margin = margin(r = 10), hjust=0.5),
#         axis.text.x = element_blank(),
#         text=element_text(size=16,  family="serif"),
#         plot.title = element_text(family="serif",hjust = 0.5, size=20))+
#   labs(x="", y = "Mean score")+
#   theme(legend.position="none")

blade.meancrps <- c(result.blade.meancrps[,1],result.blade.meancrps[,2],c(0.8665093009231732, 3.996773266725916, 0.9802456550518821, 0.6696942386202817, 0.8310129570466523, 
                                                                          0.6079964990070106, 0.8242488435136438, 0.7511318502757607, 0.6535410076500228, 1.103256421420058, 
                                                                          0.6275321686864954, 0.7299826779332969, 1.0389419776713766, 2.146308417291548, 0.7090765223645336, 
                                                                          0.7745347977930385, 0.4509185498132722, 0.5678350412189771, 0.7818611893871226, 1.0411421605781788,
                                                                          0.7248813654899693, 0.6415190011288245, 1.316842718739142, 0.6781357965493838, 0.8546743101841922, 
                                                                          0.8404619653981069, 0.8762054608835957, 0.855108086986757, 0.8582863606674513, 0.6732108401997121, 
                                                                          0.967623455683443, 0.6701625579856461, 0.8315109786094997, 0.6687164971636793, 1.1600744842502486, 
                                                                          0.8909211338068274, 0.7130320995639049, 0.6955417205478769, 0.6198301378860349, 0.9090131335823368,
                                                                          0.6687303161557062, 0.6335639385687759, 0.6685754180500859, 0.7814467739024099, 0.7310736180116385, 
                                                                          0.7622633520528582, 0.9875835837418133, 1.1314801482714438, 0.8567294136300516, 0.6545641497112342))

p15 <- ggplot(data.frame(blade.meancrps, model), aes(x=model, y=blade.meancrps, fill=model)) + 
  geom_boxplot(alpha=0.5)  + 
  theme_ipsum() +  ggtitle("Blade data") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.title.y = element_text(size=18, margin = margin(r = 10), hjust=0.5),
        axis.text.x = element_blank(),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="", y = "CRPS")+
  theme(legend.position="none")


apply(matrix(blade.meancrps, nrow=50), 2, mean) # 0.9119910 1.7862862 0.8894861

blade.comptime <- c(result.blade.comptime[,1],result.blade.comptime[,2],c(59.946934938430786, 52.38716220855713, 45.430317878723145, 53.316006898880005, 48.21130299568176, 
                                                                          56.10589289665222, 44.731395959854126, 47.18866419792175, 60.56380605697632, 57.43849778175354, 
                                                                          48.32838797569275, 53.80115723609924, 45.838510274887085, 46.60700488090515, 45.262287855148315, 
                                                                          40.9620361328125, 34.07575297355652, 57.633909940719604, 43.846577167510986, 44.71290373802185, 
                                                                          57.686196088790894, 50.56412315368652, 62.51956605911255, 51.00797200202942, 53.17219305038452, 
                                                                          51.76465892791748, 58.04509210586548, 109.99326395988464, 47.08927607536316, 47.45210385322571, 
                                                                          41.366633892059326, 49.39511013031006, 57.91219902038574, 45.636048316955566, 39.643718242645264, 
                                                                          50.07098174095154, 46.20361399650574, 63.02182722091675, 60.79853701591492, 44.08073091506958,
                                                                          48.222208976745605, 53.435014963150024, 41.8203980922699, 50.117005825042725, 52.72865009307861, 
                                                                          52.23375916481018, 33.32086992263794, 55.11390280723572, 46.73026132583618, 51.14283609390259))

p16 <- ggplot(data.frame(blade.comptime, model), aes(x=model, y=blade.comptime, fill=model)) + 
  geom_boxplot(alpha=0.5)  + 
  theme_ipsum() +  ggtitle("Blade data") +
  theme(axis.title.x = element_text(margin = margin(t = 10), hjust=0.5),
        axis.title.y = element_text(size=18, margin = margin(r = 10), hjust=0.5),
        axis.text.x = element_blank(),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="", y = "Computation time")+
  theme(legend.position="none")


apply(matrix(blade.comptime, nrow=50), 2, mean) # 0.27988  0.06636 51.17355

p.blade <- ggarrange(p13, p15, p16, ncol=3, nrow=1, common.legend = TRUE, legend="bottom")
p.blade

mean(blade.l2error[1:50]) # 1.825184
mean(blade.l2error[51:100]) # 3.213079
mean(blade.l2error[101:150]) # 1.739078
  
# mean(blade.meanscore[1:50]) # -5.091343
# mean(blade.meanscore[51:100]) # -5.496169
# mean(blade.meanscore[101:150]) # -2.10966

mean(blade.meancrps[1:50]) # 0.911991
mean(blade.meancrps[51:100]) # 1.786286
mean(blade.meancrps[101:150]) # 0.8894861

mean(blade.comptime[1:50]) # 0.27988
mean(blade.comptime[51:100]) # 0.06636
mean(blade.comptime[101:150]) # 51.17355



### Blade active learning ###

costc <-  c(0.00,  0.62,  1.24,  2.89,  3.51,  4.13,  4.75,
            5.37,  5.99,  6.61,  8.26,  8.88, 10.53, 12.18,
            13.83, 15.48, 17.13, 18.78, 20.43, 22.08, 23.73,
            25.38)

rmsec <- c(3.2611319, 1.9233651, 1.8363613, 1.5997293,
           1.3733761, 1.3515449, 1.3865496, 1.3740616,
           1.3680320, 1.3767876, 1.2847067, 1.3014988,
           1.1535156, 0.8852040, 0.8428560, 0.8294780,
           0.7333394, 0.7886327, 0.7166451, 0.7190991,
           0.6916721, 0.7031986)

costc2 <- c(0.00,  0.62,  1.24,  2.89,  3.51,  4.13,  5.78,
            6.40,  7.02,  7.64,  8.26,  8.88,  9.50, 11.15,
            11.77, 12.39, 13.01, 13.63, 14.25, 14.87, 15.49,
            16.11, 16.73, 17.35, 17.97, 18.59, 19.21, 19.83,
            20.45, 21.07, 21.69, 22.31, 22.93, 23.55, 24.17,
            24.79, 25.41)

rmsec2 <- c(3.2611319, 1.9315585, 1.6128726, 1.4028079,
            1.3854075, 1.4223574, 1.3481041, 1.3464402,
            1.3205793, 1.1906808, 1.0899642, 0.9940933,
            1.0188526, 1.1885785, 0.9822382, 1.0201724,
            1.0693308, 1.0616464, 0.8565424, 0.6788938,
            0.5354022, 0.4947483, 0.5140573, 0.4683984,
            0.4608091, 0.4805101, 0.4781361, 0.5405246,
            0.5066524, 0.3588916, 0.3694267, 0.3407210,
            0.3538892, 0.3643914, 0.3408677, 0.3353673,
            0.3198383)

costco <- c(0.00, 0.62, 1.24, 1.86, 2.48, 3.10, 3.72,
            4.34, 4.96, 5.58, 6.20, 6.82, 7.44, 8.06,
            8.68, 9.30, 9.92, 10.54, 11.16, 11.78, 12.40,
            13.02, 13.64, 14.26, 14.88, 15.50, 16.12, 16.74,
            17.36, 17.98, 18.60, 19.22, 19.84, 20.46, 21.08,
            21.70, 22.32, 22.94, 23.56, 24.18, 24.80, 25.42)

rmseco <- c(6.5116637, 1.6451782, 1.1920471, 1.2482818,
            1.3037054, 1.3195764, 1.3250443, 1.3514703,
            1.3316573, 1.3884059, 7.5710551, 1.2704820,
            1.2402900, 1.3467219, 1.2394631, 1.0044693,
            1.1315698, 1.0017485, 7.4709064, 7.0190235,
            6.9817714, 8.9545034, 7.4590203, 3.1407430,
            8.4408065, 6.8334551, 3.6085404, 0.5124625,
            6.8825590, 2.5887414, 7.0842701, 6.4936995,
            6.4759153, 6.4558973, 2.0079818, 0.4453532,
            6.4477626, 6.1542354, 6.6337098, 2.8434502,
            6.0050222, 1.1241586) #* (sqrt(sum((y.test)^2))) / sqrt(100) # 203.0631 / sqrt(100)

costk <- c(0.00, 0.62, 2.27, 3.92, 5.57, 7.22, 8.87,
           10.52, 12.17, 13.82, 14.44, 16.09, 17.74, 19.39,
           21.04, 22.69, 24.34, 25.99)

rmsek <- c(2.523179, 2.430208, 2.530261, 1.972223, 2.022362,
           2.317552, 2.616220, 2.541949, 2.628952, 4.241900,
           4.266809, 4.254287, 3.486376, 3.444703, 2.468400,
           2.614164, 3.131806, 3.000546) #* (sqrt(sum((y.test)^2))) / sqrt(100) # 203.0631 / sqrt(100)


parkalresult <- data.frame(x=rep(0:100,4), 
                           Strategy=factor(c(rep("Comprehensive searching",101), rep("Optimal point-first",101), rep("Cokriging",101), rep("MR-SUR",101)), 
                                           levels=c("Comprehensive searching","Optimal point-first","Cokriging","MR-SUR")),
                           Mean = c(resultmeanc, resultmeanc2, resultmeanco, resultmeank),
                           Max = c(resultmaxc, resultmaxc2, resultmaxco, resultmaxk),
                           Min = c(resultminc, resultminc2, resultminco, resultmink))

resultdfc <- data.frame(cost=costc, error=rmsec)
resultdfc2 <- data.frame(cost=costc2, error=rmsec2)
resultdfco <- data.frame(cost=costco, error=rmseco)
resultdfk <- data.frame(cost=costk, error=rmsek)

Strategy=factor(c(rep("Comprehensive searching",nrow(resultdfc)), rep("Optimal point-first",nrow(resultdfc2)), 
                  rep("Cokriging",nrow(resultdfco)), rep("MR-SUR",nrow(resultdfk))), 
                levels=c("Comprehensive searching","Optimal point-first","Cokriging","MR-SUR"))
resultdf <- data.frame(rbind(resultdfc, resultdfc2,resultdfco, resultdfk), Strategy)


ggplot(data=resultdf, aes(x=cost, y=error), color=group) +
  theme_ipsum() +  ggtitle("Active learning for Blade data") +
  theme(axis.title.x = element_text(size=14, margin = margin(t = 10), hjust=0.5),
        axis.title.y = element_text(size=14, margin = margin(r = 10), hjust=0.5),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="Costs", y = "RMSE") +
  geom_step(aes(group=Strategy, color=Strategy, linetype=Strategy)) +
  scale_color_manual(
    values = c("red", "green", "blue", "purple"),
    labels = c("Comprehensive searching", "Optimal point-first",
               "Co-kriging", "MR-SUR")
  )+
  scale_linetype_manual(values=c(1,2,6,5))+
  scale_color_manual(values=c("red", "green", "blue", "purple"))
  


### Comparison to recursive co-kriging, Qian(2008) ###

sqrt(sum(((predy-ytest)/ytest)^2)/12) # closed form 0.1904997
sqrt(sum(((pred.muficokm$mean-ytest)/ytest)^2)/12) # Cokm 0.513249
sqrt(sum(((pred.NARGP-ytest)/ytest)^2)/12) # NARGP 0.141051
sqrt(sum(((pred.qian-ytest)/ytest)^2)/12) # Qian 0.3873585

qiancomp <- data.frame(y = rep(ytest,4), Prediction = c(predy, pred.muficokm$mean, pred.NARGP, pred.qian), group = factor(c(rep("Proposed",12),rep("Recursive co-kriging",12),rep("NARGP",12),rep("BHGP",12)), 
                                                                                                              levels=c("Proposed","Recursive co-kriging","NARGP","BHGP")))


ggplot(data=qiancomp, aes(x=Prediction, y=y, group=group)) +
  theme_ipsum() +  ggtitle("Model comparison with linear cellular alloys data") +
  theme(axis.title.x = element_text(size=14, margin = margin(t = 10), hjust=0.5),
        axis.title.y = element_text(size=14, margin = margin(r = 10), hjust=0.5),
        text=element_text(size=16,  family="serif"),
        plot.title = element_text(family="serif",hjust = 0.5, size=20))+
  labs(x="Prediction", y = "True y") +
  geom_point(aes(x=Prediction, y=y, color=group, shape=group), size=2.7) +
  geom_abline(intercept = 0, slope = 1) + 
  xlim(0,42) + ylim(0,42)




  