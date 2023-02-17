library(rreliability)

blankbackground <-   theme(
  panel.background = element_rect(fill = "transparent"), # bg of the panel
  plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
  panel.grid.major = element_blank(), # get rid of major grid
  panel.grid.minor = element_blank(), # get rid of minor grid
  legend.background = element_rect(fill = "transparent"), # get rid of legend bg
  legend.box.background = element_blank()#element_rect(fill = "transparent", colour="none") # get rid of legend panel bg
)

m1 = c(1,5,10,20,50,40,50,60,70,80,90,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800,850,900,950,1000)
m2 = c(8,16,30,24,39,54,40,68,72,62,122,80,181,259,275,380,320,434,479,587,626,648,738,766,793,851,871,957,1001,960)
m3 = m2 + round(runif(30)*1)+100

data = measures(m1, m2)
(ba <- bland_altman(data, method = 2))
(ba2 <- bland_altman(data, method = 3))
(co <- correlation(data))
cor(data)
(co2 <- rankcorrelation(data))
cor(data, method = "spearman")
(i <- rreliability::icc(data))
(i <- rreliability::icc(data, method="randeff"))
irr::icc(data, type="agreement")$value

(co <- correlation(as.measures(data[1:10,])))
(co <- correlation(as.measures(data[11:20,])))
(co <- correlation(as.measures(data[21:30,])))
(co <- correlation(data))


(ic <- rreliability::icc(as.measures(data[1:10,])))
(ic <- rreliability::icc(as.measures(data[11:20,])))
(ic <- rreliability::icc(as.measures(data[21:30,])))
(ic <- rreliability::icc(data))

(re <- rreliability::regression(as.measures(data[1:10,])))
(re <- rreliability::regression(as.measures(data[11:20,])))
(re <- rreliability::regression(as.measures(data[21:30,])))
(re <- rreliability::regression(data))

plot(data)
plot(ba)
#psych::ICC(data)$results$ICC[1]  #[1] 0.9918106 0.9918228 0.9947956 0.9958884 0.9958946 0.9973910
p = plot(ba, data = TRUE, ideal=FALSE, bias = FALSE, error = FALSE, meanerror = FALSE, regline = FALSE)
p = plot(ba, data = TRUE, ideal=TRUE, bias = FALSE, error = FALSE, meanerror = FALSE, regline = FALSE)
p = plot(ba, data = TRUE, ideal=TRUE, bias = TRUE, error = FALSE, meanerror = FALSE, regline = FALSE)
p = plot(ba, data = TRUE, ideal=TRUE, bias = TRUE, error = TRUE, meanerror = FALSE, regline = FALSE)
p = plot(ba, data = TRUE, ideal=TRUE, bias = TRUE, error = TRUE, meanerror = TRUE, regline = FALSE)
p = plot(ba, data = TRUE, ideal=TRUE, bias = TRUE, error = TRUE, meanerror = TRUE, regline = TRUE)

p = plot(ba)

