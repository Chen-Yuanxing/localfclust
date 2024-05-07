
rm(list = ls())

########################################################################                                  
############  fit the LFC model to the COVID-19 dataset     ############ 
########################################################################   
library(localfclust)
library(fda)

setwd("/home/yc947/project")  # set working directory to the folder of this Rdata object
load("covid19_normalized_nytimes.RData")  # load the normalized covid-19 data

#row.names(covid19.normalized)

n <- nrow(covid19.normalized)
n0 <- ncol(covid19.normalized)
nknots <- 30
L <- nknots + 1
r <- 3
order <- r + 1
p <- nknots + order



basisobj <- create.bspline.basis(c(0,1), nbasis = p, norder = order)

Times <- seq(0, 1, length.out = n0)
times <- vector(mode = "list", length = n)
for (i in 1:n) {
  times[[i]] <- Times
}

Z_lst <- vector(mode = "list", length = n)
for (i in 1:n) {
  Z_lst[[i]] <- eval.basis(Times, basisobj)
}

Y_lst <- vector(mode = "list", length = n)
for (i in 1:n) {
  Y_lst[[i]] <- covid19.normalized[i,]
}


D <- eval.penalty(basisobj = basisobj, Lfdobj=int2Lfd(2))




lambda1.opt <- 3e-7
lambda2.opt <- 0.1
lambda3.opt <- 0.3


model.lfc <- LFC.fit(times = times, oridata = Y_lst, rangeval = c(0, 1), nknots = nknots, order = order, nu = 2,
                     tau = 3, K0 = 8, lambda1 = lambda1.opt, lambda2 = lambda2.opt, lambda3 = lambda3.opt)

model.lfc$cls_mem








###########################################################################################                                   
############        plot a graph with 6 global clusters, which corresponds to  ############ 
############        Figure 5 in Section 6 (Chen et al., 2024)                  ############ 
########################################################################################### 
library(ggplot2)

state <- row.names(covid19.normalized)
cluster1.curve <- as.numeric(Z_lst[[1]] %*% model.lfc$centers[,1])
cluster2.curve <- as.numeric(Z_lst[[1]] %*% model.lfc$centers[,2])
cluster3.curve <- as.numeric(Z_lst[[1]] %*% model.lfc$centers[,3])
cluster4.curve <- as.numeric(Z_lst[[1]] %*% model.lfc$centers[,4])
cluster5.curve <- as.numeric(Z_lst[[1]] %*% model.lfc$centers[,5])
cluster6.curve <- as.numeric(Z_lst[[1]] %*% model.lfc$centers[,6])



Y_hat_mat <- matrix(0, nrow = n, ncol = n0)
for (i in 1:n) {
  Y_hat_mat[i,] <- as.numeric(Z_lst[[i]] %*% model.lfc$Beta_ini[,i])
}



# generate the curves of Cluster 1

covid19.group1 <- data.frame(
  Date = rep(seq(from = as.Date("2020-03-08"), to = as.Date("2021-02-22") ,by = "day"), 23),  # date-time series
  State = rep(c(state[model.lfc$cls_mem == 1], "Center"), each = 352),  # classification lable via states
  value = as.vector(cbind(t(Y_hat_mat[model.lfc$cls_mem == 1,]), cluster1.curve))  # repeated measurements
)
covid19.group1$State <- factor(covid19.group1$State, levels = c(state[model.lfc$cls_mem == 1], "Center"))



p1 <- ggplot(covid19.group1, aes(x = Date, y = value, color = State, linetype = State, linewidth = State)) +
  geom_line() + scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("Washington" = "#FF0000", "California" = "#FFC0CB", "Arizona" = "#FF4500", "Massachusetts" = "#DC143C",
                                "Oregon" = "#00FF00", "New York" = "#008000", "Rhode Island" = "#32CD32",
                                "New Hampshire" = "#0000FF", "North Carolina" = "#ADD8E6", "New Jersey" = "#000080", "Tennessee" = "#4169E1",
                                "Kentucky" = "#FFFF00", "Oklahoma" = "#FFD700", "Pennsylvania" = "#FFA500", "Vermont" = "#FFCC00",
                                "Virginia" = "#800080", "Connecticut" = "#DA70D6", "Ohio" = "#9400D3", "Arkansas" = "#8A2BE2",
                                "Delaware" = "#008080", "Maine" = "#C0C0C0", "West Virginia" = "#CD7F32", "Center" = "black")) +
  scale_linetype_manual(values = c("Washington" = "dotted", "California" = "dashed", "Arizona" = "dotdash", "Massachusetts" = "longdash",
                                   "Oregon" = "twodash", "New York" = "dotted", "Rhode Island" = "dashed", 
                                   "New Hampshire" = "dotdash", "North Carolina" = "longdash", "New Jersey" = "twodash", "Tennessee" = "dotted",
                                   "Kentucky" = "dashed", "Oklahoma" = "dotdash", "Pennsylvania" = "longdash", "Vermont" = "twodash",
                                   "Virginia" = "dotted", "Connecticut" = "dashed", "Ohio" = "dotdash", "Arkansas" = "longdash",
                                   "Delaware" = "twodash", "Maine" = "dotted", "West Virginia" = "dashed", "Center" = "solid")) +
  scale_y_continuous(breaks = seq(-0.5, 5, 0.5), limits = c(-0.5, 5)) +
  ylab("Normalized counts of newly confirmed cases") + 
  scale_linewidth_manual(values = c(rep(1.5,22),3)) +
  theme(axis.title.y = element_text(size = 25), axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), axis.title.x = element_text(size = 25)) + theme(legend.text = element_text(size = 20)) +
  theme_bw(base_size = 40) + guides(color = guide_legend(ncol = 1)) + theme(legend.key.size = unit(3, "lines"))

p1

ggsave("COVID_cluster1.pdf", plot = p1, device = "pdf", width = 20, height = 16)






# generate the curves of Cluster 2
covid19.group2 <- data.frame(
  Date = rep(seq(from = as.Date("2020-03-08"), to = as.Date("2021-02-22") ,by = "day"), 11),  
  State = rep(c(state[model.lfc$cls_mem == 2], "Center"), each = 352),  
  value = as.vector(cbind(t(Y_hat_mat[model.lfc$cls_mem == 2,]), cluster2.curve))  
)
covid19.group2$State <- factor(covid19.group2$State, levels = c(state[model.lfc$cls_mem == 2], "Center"))




p2 <- ggplot(covid19.group2, aes(x = Date, y = value, color = State, linetype = State, linewidth = State)) +
  geom_line() + scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("Illinois" = "#08519C", "Utah" = "#BDD7E7", "Colorado" = "#3182BD", "Indiana" = "#A50F15",
                                "Kansas" = "#74C476", "Missouri" = "#E6550D", "Michigan" = "#9E9AC8", "New Mexico" = "#54278F",
                                "Alaska" = "#008080",     "Idaho" = "#9400D3", "Center" = "black")) +
  scale_linetype_manual(values = c("Illinois" = "dotted", "Utah" = "dashed", "Colorado" = "dotdash", "Indiana" = "longdash",
                                   "Kansas" = "twodash", "Missouri" = "dotted", "Michigan" = "dashed",
                                   "New Mexico" = "dotdash", "Alaska" = "longdash", "Idaho" = "twodash", "Center" = "solid")) +
  ylab("Normalized counts of newly confirmed cases") + 
  scale_y_continuous(breaks = seq(-0.5, 5, 0.5), limits = c(-0.5, 5)) +
  scale_linewidth_manual(values = c(rep(1.5,10),3)) +
  theme(axis.title.y = element_text(size = 25), axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), axis.title.x = element_text(size = 25)) + theme(legend.text = element_text(size = 20)) +
  theme_bw(base_size = 40) + guides(color = guide_legend(ncol = 1)) + theme(legend.key.size = unit(3, "lines"))

p2

ggsave("COVID_cluster2.pdf", plot = p2, device = "pdf", width = 20, height = 16)









# generate the curves of Cluster 3
covid19.group3 <- data.frame(
  Date = rep(seq(from = as.Date("2020-03-08"), to = as.Date("2021-02-22") ,by = "day"), 9),  
  State = rep(c(state[model.lfc$cls_mem == 3], "Center"), each = 352),  
  value = as.vector(cbind(t(Y_hat_mat[model.lfc$cls_mem == 3,]), cluster3.curve))  
)
covid19.group3$State <- factor(covid19.group3$State, levels = c(state[model.lfc$cls_mem == 3], "Center"))




p3 <- ggplot(covid19.group3, aes(x = Date, y = value, color = State, linetype = State, linewidth = State)) +
  geom_line() + scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("Wisconsin" = "#08519C", "Nebraska" = "#BDD7E7", "Minnesota" = "#3182BD", "Iowa" = "#A50F15",
                                "South Dakota" = "#74C476", "North Dakota" = "#E6550D", "Wyoming" = "#9E9AC8", "Montana" = "#54278F",
                                "Center" = "black")) +
  scale_linetype_manual(values = c("Wisconsin" = "dotted", "Nebraska" = "dashed", "Minnesota" = "dotdash", 
                                   "Iowa" = "dotted", "South Dakota" = "dashed", "North Dakota" = "dotdash",
                                   "Wyoming" = "longdash", "Montana" = "twodash", "Center" = "solid")) +
  ylab("Normalized counts of newly confirmed cases") + 
  scale_y_continuous(breaks = seq(-0.5, 5, 0.5), limits = c(-0.5, 5)) +
  scale_linewidth_manual(values = c(rep(1.5,8),3)) +
  theme(axis.title.y = element_text(size = 25), axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), axis.title.x = element_text(size = 25)) + theme(legend.text = element_text(size = 20)) +
  theme_bw(base_size = 40) + guides(color = guide_legend(ncol = 1)) + theme(legend.key.size = unit(3, "lines"))

p3

ggsave("COVID_cluster3.pdf", plot = p3, device = "pdf", width = 20, height = 16)





# generate the curves of Cluster 4
covid19.group4 <- data.frame(
  Date = rep(seq(from = as.Date("2020-03-08"), to = as.Date("2021-02-22") ,by = "day"), 9),  
  State = rep(c(state[model.lfc$cls_mem == 4], "Center"), each = 352),  
  value = as.vector(cbind(t(Y_hat_mat[model.lfc$cls_mem == 4,]), cluster4.curve))  
)
covid19.group4$State <- factor(covid19.group4$State, levels = c(state[model.lfc$cls_mem == 4], "Center"))




p4 <- ggplot(covid19.group4, aes(x = Date, y = value, color = State, linetype = State, linewidth = State)) +
  geom_line() + scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("Texas" = "#08519C", "Florida" = "#BDD7E7", "Georgia" = "#3182BD", "Nevada" = "#A50F15",
                                "South Carolina" = "#74C476", "Louisiana" = "#E6550D", "Mississippi" = "#9E9AC8", "Alabama" = "#54278F",
                                "Center" = "black")) +
  scale_linetype_manual(values = c("Texas" = "dotted", "Florida" = "dashed", "Georgia" = "dotdash", "Nevada" = "longdash",
                                   "South Carolina" = "twodash", "Louisiana" = "dotted", "Mississippi" = "dashed",
                                   "Alabama" = "dotdash", "Center" = "solid")) +
  ylab("Normalized counts of newly confirmed cases") + 
  scale_y_continuous(breaks = seq(-0.5, 5, 0.5), limits = c(-0.5, 5)) + 
  scale_linewidth_manual(values = c(rep(1.5,8),3)) +
  theme(axis.title.y = element_text(size = 25), axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), axis.title.x = element_text(size = 25)) + theme(legend.text = element_text(size = 20)) +
  theme_bw(base_size = 40) + guides(color = guide_legend(ncol = 1)) + theme(legend.key.size = unit(3, "lines"))

p4

ggsave("COVID_cluster4.pdf", plot = p4, device = "pdf", width = 20, height = 16)








# generate the curves of Cluster 5
covid19.group5 <- data.frame(
  Date = rep(seq(from = as.Date("2020-03-08"), to = as.Date("2021-02-22") ,by = "day"), 3), 
  State = rep(c(state[model.lfc$cls_mem == 5], "Center"), each = 352),  
  value = as.vector(cbind(t(Y_hat_mat[model.lfc$cls_mem == 5,]), cluster5.curve)) 
)
covid19.group5$State <- factor(covid19.group5$State, levels = c(state[model.lfc$cls_mem == 5], "Center"))



p5 <- ggplot(covid19.group5, aes(x = Date, y = value, color = State, linetype = State, linewidth = State)) +
  geom_line() + scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("Maryland" = "#08519C", "District of Columbia" = "#E6550D", "Center" = "black")) +
  scale_linetype_manual(values = c("Maryland" = "dotted", "District of Columbia" = "dashed", "Center" = "solid")) +
  ylab("Normalized counts of newly confirmed cases") + 
  scale_y_continuous(breaks = seq(-0.5, 5, 0.5), limits = c(-0.5, 5)) + 
  scale_linewidth_manual(values = c(rep(1.5,2),3)) +
  theme(axis.title.y = element_text(size = 25), axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), axis.title.x = element_text(size = 25)) + theme(legend.text = element_text(size = 20)) +
  theme_bw(base_size = 40) + guides(color = guide_legend(ncol = 1)) + theme(legend.key.size = unit(3, "lines"))

p5

ggsave("COVID_cluster5.pdf", plot = p5, device = "pdf", width = 20, height = 16)





# generate the curves of Cluster 6
covid19.group6 <- data.frame(
  Date = rep(seq(from = as.Date("2020-03-08"), to = as.Date("2021-02-22") ,by = "day"), 2),  
  State = rep(c("Hawaii", "Center"), each = 352),  
  value = as.vector(cbind(Y_hat_mat[model.lfc$cls_mem == 6,], cluster6.curve))  
)
covid19.group6$State <- factor(covid19.group6$State, levels = c("Hawaii", "Center"))


p6 <- ggplot(covid19.group6, aes(x = Date, y = value, color = State, linetype = State, linewidth = State)) +
  geom_line() + scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("Hawaii" = "red", "Center" = "black")) +
  scale_linetype_manual(values = c("Hawaii" = "dotted", "Center" = "solid")) +
  ylab("Normalized counts of newly confirmed cases") + 
  scale_y_continuous(breaks = seq(-0.5, 5, 0.5), limits = c(-0.5, 5)) + 
  scale_linewidth_manual(values = c(rep(1.5,1),3)) +
  theme(axis.title.y = element_text(size = 25), axis.text.x = element_text(size = 20), 
        axis.text.y = element_text(size = 20), axis.title.x = element_text(size = 25)) + theme(legend.text = element_text(size = 20)) +
  theme_bw(base_size = 40) + guides(color = guide_legend(ncol = 1)) + theme(legend.key.size = unit(3, "lines"))

p6

ggsave("COVID_cluster6.pdf", plot = p6, device = "pdf", width = 20, height = 16)






#######################################################################################################                                
############        plot a graph to show different local clustering structures,            ############ 
############        which corresponds to Figure 6 in Section 6 (Chen et al., 2024)         ############ 
####################################################################################################### 

#install.packages("latex2exp")
library(latex2exp)


covid19.cluster.centers <- data.frame(
  Date = rep(seq(from = as.Date("2020-03-08"), to = as.Date("2021-02-22") ,by = "day"), 5),  # 日期列
  Clusters = rep(c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"), each = 352),  # 时间序列标识列
  value = as.vector(c(cluster1.curve, cluster2.curve, cluster3.curve, cluster4.curve, cluster5.curve))  # 数值列
)
covid19.cluster.centers$Clusters <- factor(covid19.cluster.centers$Clusters, levels = c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5"))



# 创建图表
p7 <- ggplot(covid19.cluster.centers, aes(x = Date, y = value, color = Clusters, linetype = Clusters, size = Clusters)) +
  geom_line() + scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  scale_color_manual(values = c("Cluster 1" = "#3182BD", "Cluster 2" = "#E6550D", "Cluster 3" = "purple", "Cluster 4" = "#A50F15",
                                "Cluster 5" = "#74C476")) +
  scale_linetype_manual(values = c("Cluster 1" = "dashed", "Cluster 2" = "dotted", "Cluster 3" = "longdash", "Cluster 4" = "dotdash",
                                   "Cluster 5" = "twodash")) +
  scale_y_continuous(breaks = seq(0, 4, 0.5), limits = c(0, 4)) + 
  scale_size_manual(values = rep(1,5)) + 
  geom_vline(xintercept = as.Date("2020-03-08"), linetype = "dotted", color = "black") +
  geom_vline(xintercept = as.Date("2020-04-12"), linetype = "dotted", color = "black") +
  geom_vline(xintercept = as.Date("2020-06-18"), linetype = "dotted", color = "black") +
  geom_vline(xintercept = as.Date("2020-08-12"), linetype = "dotted", color = "black") +
  geom_vline(xintercept = as.Date("2020-10-10"), linetype = "dotted", color = "black") +
  geom_vline(xintercept = as.Date("2020-11-13"), linetype = "dotted", color = "black") +
  geom_vline(xintercept = as.Date("2021-01-28"), linetype = "dotted", color = "black") +
  geom_vline(xintercept = as.Date("2021-02-22"), linetype = "dotted", color = "black") +
  annotate("text", x = as.Date("2020-03-26"), y = 4, label = expression(T[1]), color = "blue", size = 10, parse = TRUE) +
  annotate("text", x = as.Date("2020-05-15"), y = 4, label = expression(T[2]), color = "blue", size = 10, parse = TRUE) +
  annotate("text", x = as.Date("2020-07-15"), y = 4, label = expression(T[3]), color = "blue", size = 10, parse = TRUE) +
  annotate("text", x = as.Date("2020-09-12"), y = 4, label = expression(T[4]), color = "blue", size = 10, parse = TRUE) +
  annotate("text", x = as.Date("2020-10-26"), y = 4, label = expression(T[5]), color = "blue", size = 10, parse = TRUE) +
  annotate("text", x = as.Date("2020-12-20"), y = 4, label = expression(T[6]), color = "blue", size = 10, parse = TRUE) +
  annotate("text", x = as.Date("2021-02-10"), y = 4, label = expression(T[7]), color = "blue", size = 10, parse = TRUE) +
  ylab("Normalized counts of newly confirmed cases") + 
  theme(axis.title.y = element_text(size = 20), axis.text.x = element_text(size = 15), 
        axis.text.y = element_text(size = 15), axis.title.x = element_text(size = 20)) + theme(legend.text = element_text(size = 15), legend.position = "top") +
  theme_bw(base_size = 25) + theme(legend.key.size = unit(4, "lines"))
# geom_text(x = as.Date("2020-05-21"), y = 4, label = expression(paste(sigma, ' = 1')), color = "blue", size = 8) 

p7

ggsave("COVID_local_clusters.pdf", plot = p7, device = "pdf", width = 15, height = 10)





#######################################################################################################                                
############      plot a graph used as a motivation example in Introduction Section        ############ 
####################################################################################################### 

Date <- seq(from = as.Date("2020-03-08"), to = as.Date("2021-02-22"), by = "day")
monthly_dates <- c(as.Date("2020-03-01"), Date[format(Date, "%d") == "01"], as.Date("2021-03-01"))
matplot(as.numeric(Date), t(Y_hat_mat[-22,]),type = "l", ylim = c(-0.5,5), xaxt = "n", yaxt = "n", xlab = "Date",
        ylab = "Normalized counts of newly confirmed cases", cex.lab = 1.5, lwd = 2)

axis(1, at = monthly_dates, labels = format(monthly_dates, "%b"), cex.axis = 1.5)
axis(2, at = seq(-0.5, 5, by = 0.5), cex.axis = 1.5)


abline(v = monthly_dates, col = "gray", lty = "solid")
abline(h = seq(-0.5, 5, by = 0.5), col = "gray", lty = "solid")


abline(v = as.Date("2020-03-08"), col = "blue", lty = "dashed", lwd = 2)
abline(v = as.Date("2020-06-15"), col = "blue", lty = "dashed", lwd = 2)
abline(v = as.Date("2020-08-15"), col = "blue", lty = "dashed", lwd = 2)
abline(v = as.Date("2020-10-15"), col = "blue", lty = "dashed", lwd = 2)
abline(v = as.Date("2021-02-22"), col = "blue", lty = "dashed", lwd = 2)










